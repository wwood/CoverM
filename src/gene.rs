use std;
use std::collections::HashMap;
use std::io::BufRead;

use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use bam_generator::*;
use coverage_takers::*;
use mosdepth_genome_coverage_estimators::*;
use nm;
use FlagFilter;
use ReadsMapped;

/// A gene (or, more generally, a GFF feature) defined as a region of a contig.
#[derive(Clone, Debug, PartialEq)]
pub struct Gene {
    /// Name used in the output. Taken from the GFF `ID` attribute where
    /// available, otherwise generated from the coordinates.
    pub name: String,
    /// Name of the contig (GFF seqid / column 1) the gene is located on.
    pub contig: String,
    /// 0-based, inclusive start of the gene on the contig.
    pub start: u64,
    /// 0-based, exclusive end of the gene on the contig.
    pub end: u64,
}

/// Read one or more GFF3 files, returning one [`Gene`] per feature line. Lines
/// beginning with `#` and blank lines are ignored, as is the FASTA section of a
/// GFF3 file (everything after a `##FASTA` directive).
pub fn read_gff_files(gff_files: &[&str]) -> Vec<Gene> {
    let mut genes = vec![];
    for gff_file in gff_files {
        debug!("Reading GFF file {gff_file}");
        let file = std::fs::File::open(gff_file)
            .unwrap_or_else(|e| panic!("Failed to open GFF file '{}': {}", gff_file, e));
        let reader = std::io::BufReader::new(file);
        let mut num_read = 0;
        for line_result in reader.lines() {
            let line = line_result
                .unwrap_or_else(|e| panic!("Error reading GFF file '{}': {}", gff_file, e));
            let trimmed = line.trim_end();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with('#') {
                // The FASTA section that can be appended to GFF3 files starts
                // with a ##FASTA directive - stop parsing this file there.
                if trimmed.starts_with("##FASTA") {
                    break;
                }
                continue;
            }
            let fields: Vec<&str> = trimmed.split('\t').collect();
            if fields.len() < 9 {
                warn!(
                    "Ignoring GFF line in '{gff_file}' with fewer than 9 tab-separated fields: {trimmed}"
                );
                continue;
            }
            let contig = fields[0].to_string();
            let start: u64 = fields[3].parse().unwrap_or_else(|_| {
                panic!(
                    "Failed to parse GFF start coordinate '{}' in {}",
                    fields[3], gff_file
                )
            });
            let end: u64 = fields[4].parse().unwrap_or_else(|_| {
                panic!(
                    "Failed to parse GFF end coordinate '{}' in {}",
                    fields[4], gff_file
                )
            });
            if start < 1 || end < start {
                warn!(
                    "Ignoring GFF line in '{gff_file}' with nonsensical coordinates ({start}, {end}): {trimmed}"
                );
                continue;
            }
            let name = parse_gff_id(fields[8]).unwrap_or_else(|| format!("{contig}_{start}_{end}"));
            // GFF is 1-based and inclusive on both ends; convert to 0-based
            // half-open [start, end).
            genes.push(Gene {
                name,
                contig,
                start: start - 1,
                end,
            });
            num_read += 1;
        }
        info!("Read {num_read} features from GFF file {gff_file}");
    }
    genes
}

/// Extract a name for a feature from the GFF3 attributes column (column 9),
/// preferring the `ID` attribute and falling back to `Name`.
fn parse_gff_id(attributes: &str) -> Option<String> {
    let mut name_attribute = None;
    for attribute in attributes.split(';') {
        let attribute = attribute.trim();
        if let Some(value) = attribute.strip_prefix("ID=") {
            if !value.is_empty() {
                return Some(value.to_string());
            }
        } else if let Some(value) = attribute.strip_prefix("Name=") {
            if !value.is_empty() && name_attribute.is_none() {
                name_attribute = Some(value.to_string());
            }
        }
    }
    name_attribute
}

/// Per-gene accumulator used while iterating over a single contig's reads.
struct GeneCalculator {
    global_index: usize,
    start: u64,
    end: u64,
    name: String,
    num_mapped_reads: u64,
    total_edit_distance: u64,
    total_indels: u64,
    sum_identity: f64,
}

pub fn gene_coverage<R: NamedBamReader, G: NamedBamReaderGenerator<R>, T: CoverageTaker>(
    bam_readers: Vec<G>,
    coverage_taker: &mut T,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    print_zero_coverage_genes: bool,
    flag_filters: &FlagFilter,
    threads: u16,
    genes: &[Gene],
) -> Vec<ReadsMapped> {
    // Contig-end exclusion is meant to avoid mapping artifacts at the ends of
    // contigs. In gene mode each gene is a sub-region of a contig, so applying
    // the exclusion to gene boundaries would incorrectly trim every gene (and
    // zero out any gene shorter than twice the exclusion). Disable it here.
    for estimator in coverage_estimators.iter_mut() {
        estimator.disable_contig_end_exclusion();
    }

    let mut reads_mapped_vector = vec![];
    for bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();
        bam_generated.set_threads(threads as usize);

        let stoit_name = &(bam_generated.name().to_string());
        coverage_taker.start_stoit(stoit_name);
        let header = bam_generated.header().clone();
        let target_names = header.target_names();

        // Map contig name -> tid, then assign each relevant gene to its tid.
        let mut contig_name_to_tid: HashMap<&[u8], u32> = HashMap::new();
        for (tid, name) in target_names.iter().enumerate() {
            contig_name_to_tid.insert(name, tid as u32);
        }

        let mut genes_by_tid: Vec<Vec<GeneCalculator>> =
            (0..target_names.len()).map(|_| vec![]).collect();
        let mut num_genes_not_in_bam = 0;
        for gene in genes.iter() {
            match contig_name_to_tid.get(gene.contig.as_bytes()) {
                Some(&tid) => {
                    let contig_length = header.target_len(tid).expect("Corrupt BAM file?");
                    if gene.start >= contig_length {
                        warn!(
                            "Gene '{}' starts ({}) beyond the end of contig '{}' (length {}); skipping",
                            gene.name, gene.start, gene.contig, contig_length
                        );
                        continue;
                    }
                    let end = std::cmp::min(gene.end, contig_length);
                    genes_by_tid[tid as usize].push(GeneCalculator {
                        global_index: 0, // assigned below
                        start: gene.start,
                        end,
                        name: gene.name.clone(),
                        num_mapped_reads: 0,
                        total_edit_distance: 0,
                        total_indels: 0,
                        sum_identity: 0.0,
                    });
                }
                None => num_genes_not_in_bam += 1,
            }
        }
        if num_genes_not_in_bam > 0 {
            warn!(
                "{num_genes_not_in_bam} gene(s) were defined on contigs that are not \
                 reference sequences in the BAM file for {stoit_name}"
            );
        }
        // Sort genes on each contig by start position and assign a stable global
        // index in (tid, start) order. This index is used as the entry order id
        // so that output is emitted in a consistent order across samples.
        let mut next_global_index = 0;
        for tid_genes in genes_by_tid.iter_mut() {
            tid_genes.sort_by_key(|g| (g.start, g.end));
            for g in tid_genes.iter_mut() {
                g.global_index = next_global_index;
                next_global_index += 1;
            }
        }

        let mut record: bam::record::Record = bam::record::Record::new();
        let mut last_tid: i32 = -2; // no such tid in a real BAM file
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let mut gene_sweep_lo: usize = 0;
        let mut num_mapped_reads_total: u64 = 0;

        loop {
            match bam_generated.read(&mut record) {
                None => break,
                Some(Ok(())) => {}
                Some(e) => panic!("Error reading BAM record: {:?}", e),
            }

            if !flag_filters.passes(&record) {
                trace!("Skipping read based on flag filtering");
                continue;
            }
            if record.is_unmapped() {
                continue;
            }
            let tid = record.tid();
            if tid != last_tid {
                if tid < last_tid {
                    error!("BAM file appears to be unsorted. Input BAM files must be sorted by reference (i.e. by samtools sort)");
                    panic!("BAM file appears to be unsorted. Input BAM files must be sorted by reference (i.e. by samtools sort)");
                }
                if last_tid != -2 {
                    process_previous_contig_genes(
                        &mut genes_by_tid[last_tid as usize],
                        &ups_and_downs,
                        coverage_estimators,
                        coverage_taker,
                        print_zero_coverage_genes,
                        &mut num_mapped_reads_total,
                    );
                }
                // Emit zero coverage for genes on contigs between last_tid and
                // tid that had no reads mapped.
                if print_zero_coverage_genes {
                    let from = if last_tid == -2 {
                        0
                    } else {
                        last_tid as u32 + 1
                    };
                    for skipped_tid in from..(tid as u32) {
                        print_zero_coverage_genes_for_contig(
                            &genes_by_tid[skipped_tid as usize],
                            coverage_estimators,
                            coverage_taker,
                        );
                    }
                }
                ups_and_downs =
                    vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                gene_sweep_lo = 0;
                last_tid = tid;
            }

            let is_primary = !record.is_supplementary() && !record.is_secondary();
            if is_primary {
                num_mapped_reads_total += 1;
            }

            // Walk the cigar, updating the coverage array and computing the
            // reference span, aligned length and indels of this read.
            let read_start: u64 = record.pos() as u64;
            let mut cursor: usize = record.pos() as usize;
            let mut aligned_len: u64 = 0;
            let mut indels: u64 = 0;
            for cig in record.cigar().iter() {
                match cig {
                    Cigar::Match(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                        ups_and_downs[cursor] += 1;
                        let final_pos = cursor + cig.len() as usize;
                        if final_pos < ups_and_downs.len() {
                            ups_and_downs[final_pos] -= 1;
                        }
                        cursor += cig.len() as usize;
                        aligned_len += cig.len() as u64;
                    }
                    Cigar::Del(_) => {
                        cursor += cig.len() as usize;
                        indels += cig.len() as u64;
                        aligned_len += cig.len() as u64;
                    }
                    Cigar::RefSkip(_) => {
                        cursor += cig.len() as usize;
                    }
                    Cigar::Ins(_) => {
                        indels += cig.len() as u64;
                        aligned_len += cig.len() as u64;
                    }
                    Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
                }
            }
            let read_end: u64 = cursor as u64;

            let edit = nm(&record);
            let identity = if aligned_len > 0 {
                (aligned_len as f64 - edit as f64) / aligned_len as f64
            } else {
                0.0
            };

            // Attribute this read to every gene on the contig that it overlaps.
            // Reads arrive sorted by position, so genes that end before this
            // read starts can never overlap a later read either.
            let tid_genes = &mut genes_by_tid[tid as usize];
            while gene_sweep_lo < tid_genes.len() && tid_genes[gene_sweep_lo].end <= read_start {
                gene_sweep_lo += 1;
            }
            for gene in tid_genes[gene_sweep_lo..].iter_mut() {
                if gene.start >= read_end {
                    break;
                }
                if gene.end > read_start {
                    // The read overlaps this gene. Edit distance and identity are
                    // derived from the whole-read NM tag (which cannot be split
                    // positionally between genes), so a read that straddles a gene
                    // boundary contributes its whole-read statistics to each gene
                    // it overlaps.
                    gene.total_edit_distance += edit;
                    gene.total_indels += indels;
                    if is_primary {
                        gene.num_mapped_reads += 1;
                        if aligned_len > 0 {
                            gene.sum_identity += identity;
                        }
                    }
                }
            }
        }

        // Process the final contig with reads, then any trailing contigs.
        if last_tid != -2 {
            process_previous_contig_genes(
                &mut genes_by_tid[last_tid as usize],
                &ups_and_downs,
                coverage_estimators,
                coverage_taker,
                print_zero_coverage_genes,
                &mut num_mapped_reads_total,
            );
        }
        if print_zero_coverage_genes {
            let from = if last_tid == -2 {
                0
            } else {
                last_tid as u32 + 1
            };
            for skipped_tid in from..(target_names.len() as u32) {
                print_zero_coverage_genes_for_contig(
                    &genes_by_tid[skipped_tid as usize],
                    coverage_estimators,
                    coverage_taker,
                );
            }
        }

        let reads_mapped = ReadsMapped {
            num_mapped_reads: num_mapped_reads_total,
            num_reads: bam_generated.num_detected_primary_alignments(),
        };
        info!(
            "In sample '{}', found {} reads mapped out of {} total ({:.*}%)",
            stoit_name,
            reads_mapped.num_mapped_reads,
            reads_mapped.num_reads,
            2,
            (reads_mapped.num_mapped_reads * 100) as f64 / reads_mapped.num_reads as f64
        );
        reads_mapped_vector.push(reads_mapped);

        if bam_generated.num_detected_primary_alignments() == 0 {
            warn!(
                "No primary alignments were observed for sample {stoit_name} \
                 - perhaps something went wrong in the mapping?"
            );
        }

        bam_generated.finish();
    }
    reads_mapped_vector
}

/// Compute and emit coverage statistics for all genes on a contig, given the
/// contig's coverage (ups_and_downs) array.
fn process_previous_contig_genes<T: CoverageTaker>(
    tid_genes: &mut [GeneCalculator],
    ups_and_downs: &[i32],
    coverage_estimators: &mut [CoverageEstimator],
    coverage_taker: &mut T,
    print_zero_coverage_genes: bool,
    _num_mapped_reads_total: &mut u64,
) {
    if tid_genes.is_empty() {
        return;
    }
    // Cumulative coverage across the contig, so that a per-gene ups_and_downs
    // array can be reconstructed for any sub-region.
    let mut coverage: Vec<i32> = vec![0; ups_and_downs.len()];
    let mut cumulative_sum: i32 = 0;
    for (i, ud) in ups_and_downs.iter().enumerate() {
        cumulative_sum += ud;
        coverage[i] = cumulative_sum;
    }

    // Genes are sorted by start, which matches their global index order.
    for gene in tid_genes.iter() {
        let start = gene.start as usize;
        let end = std::cmp::min(gene.end as usize, coverage.len());
        if start >= end {
            continue;
        }
        // Reconstruct a delta-encoded ups_and_downs array for the gene region
        // whose running sum reproduces the contig coverage within the gene.
        let gene_len = end - start;
        let mut gene_ups_and_downs: Vec<i32> = vec![0; gene_len];
        gene_ups_and_downs[0] = coverage[start];
        for i in 1..gene_len {
            gene_ups_and_downs[i] = coverage[start + i] - coverage[start + i - 1];
        }

        for estimator in coverage_estimators.iter_mut() {
            estimator.setup();
            estimator.add_contig(
                &gene_ups_and_downs,
                gene.num_mapped_reads,
                gene.total_edit_distance - gene.total_indels,
                gene.sum_identity,
            );
        }
        let coverages: Vec<f32> = coverage_estimators
            .iter_mut()
            .map(|estimator| estimator.calculate_coverage(&[0]))
            .collect();
        let has_nonzero_coverage = coverages.iter().any(|&c| c > 0.0);

        if print_zero_coverage_genes || has_nonzero_coverage {
            coverage_taker.start_entry(gene.global_index, &gene.name);
            for (coverage, estimator) in coverages.into_iter().zip(coverage_estimators.iter_mut()) {
                estimator.print_coverage(coverage, coverage_taker);
            }
            coverage_taker.finish_entry();
        }
    }
    for estimator in coverage_estimators.iter_mut() {
        estimator.setup();
    }
}

/// Emit zero coverage entries for all genes on a contig that had no reads.
fn print_zero_coverage_genes_for_contig<T: CoverageTaker>(
    tid_genes: &[GeneCalculator],
    coverage_estimators: &[CoverageEstimator],
    coverage_taker: &mut T,
) {
    for gene in tid_genes.iter() {
        coverage_taker.start_entry(gene.global_index, &gene.name);
        for estimator in coverage_estimators.iter() {
            estimator.print_zero_coverage(coverage_taker, gene.end - gene.start);
        }
        coverage_taker.finish_entry();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use OutputWriter;

    fn test_with_stream<R: NamedBamReader, G: NamedBamReaderGenerator<R>>(
        expected: &str,
        bam_readers: Vec<G>,
        coverage_estimators: &mut Vec<CoverageEstimator>,
        genes: &[Gene],
        print_zero_coverage_genes: bool,
    ) -> Vec<ReadsMapped> {
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();

        let flag_filters = FlagFilter {
            include_improper_pairs: true,
            include_secondary: false,
            include_supplementary: false,
        };
        let reads_mapped_vec;
        {
            let mut coverage_taker =
                CoverageTakerType::new_single_float_coverage_streaming_coverage_printer(
                    OutputWriter::generate(Some(t)),
                );
            reads_mapped_vec = gene_coverage(
                bam_readers,
                &mut coverage_taker,
                coverage_estimators,
                print_zero_coverage_genes,
                &flag_filters,
                1,
                genes,
            );
        }
        assert_eq!(expected, std::fs::read_to_string(tf.path()).unwrap());
        reads_mapped_vec
    }

    #[test]
    fn test_parse_gff_id() {
        assert_eq!(Some("gene1".to_string()), parse_gff_id("ID=gene1;Name=abc"));
        assert_eq!(Some("abc".to_string()), parse_gff_id("Name=abc;locus=x"));
        assert_eq!(None, parse_gff_id("locus=x"));
    }

    #[test]
    fn test_read_gff_files() {
        let genes = read_gff_files(&["tests/data/genes/2seqs.gff"]);
        assert_eq!(
            vec![
                Gene {
                    name: "gene1".to_string(),
                    contig: "seq1".to_string(),
                    start: 0,
                    end: 1000,
                },
                Gene {
                    name: "gene2".to_string(),
                    contig: "seq1".to_string(),
                    start: 199,
                    end: 400,
                },
                Gene {
                    name: "gene3".to_string(),
                    contig: "seq2".to_string(),
                    start: 0,
                    end: 500,
                },
            ],
            genes
        );
    }

    #[test]
    fn test_gene_spanning_whole_contig_matches_contig_coverage() {
        // gene1 spans the whole of seq1, so its mean coverage should equal the
        // whole-contig mean coverage of seq1 (1.2, per the contig-mode tests).
        // gene3 spans the whole of seq2, which has no reads (0).
        test_with_stream(
            "2seqs.reads_for_seq1\tgene1\t1.2\n2seqs.reads_for_seq1\tgene3\t0\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            &mut vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
            &[
                Gene {
                    name: "gene1".to_string(),
                    contig: "seq1".to_string(),
                    start: 0,
                    end: 1000,
                },
                Gene {
                    name: "gene3".to_string(),
                    contig: "seq2".to_string(),
                    start: 0,
                    end: 1000,
                },
            ],
            true,
        );
    }

    #[test]
    fn test_contig_end_exclusion_disabled_in_gene_mode() {
        // Even with a non-zero contig-end-exclusion, gene mode does not trim the
        // ends of the gene, so a gene spanning the whole of seq1 reports the same
        // value (1.2) as it does with no exclusion. Without disabling exclusion,
        // 75 bases would be trimmed from each end.
        test_with_stream(
            "2seqs.reads_for_seq1\tgene1\t1.2\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            &mut vec![CoverageEstimator::new_estimator_mean(0.0, 75, false)],
            &[Gene {
                name: "gene1".to_string(),
                contig: "seq1".to_string(),
                start: 0,
                end: 1000,
            }],
            false,
        );
    }

    #[test]
    fn test_no_zeros() {
        test_with_stream(
            "2seqs.reads_for_seq1\tgene1\t1.2\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            &mut vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
            &[
                Gene {
                    name: "gene1".to_string(),
                    contig: "seq1".to_string(),
                    start: 0,
                    end: 1000,
                },
                Gene {
                    name: "gene3".to_string(),
                    contig: "seq2".to_string(),
                    start: 0,
                    end: 1000,
                },
            ],
            false,
        );
    }

    #[test]
    fn test_count_estimator_per_gene() {
        // Whole-contig read count for seq1 equals the number of reads mapped.
        let reads_mapped = test_with_stream(
            "2seqs.reads_for_seq1\tgene1\t12\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            &mut vec![CoverageEstimator::new_estimator_read_count()],
            &[Gene {
                name: "gene1".to_string(),
                contig: "seq1".to_string(),
                start: 0,
                end: 1000,
            }],
            false,
        );
        assert_eq!(
            vec![ReadsMapped {
                num_mapped_reads: 12,
                num_reads: 12,
            }],
            reads_mapped
        );
    }
}
