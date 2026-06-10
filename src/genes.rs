use std;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};

use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use bam_generator::*;
use coverage_takers::*;
use mosdepth_genome_coverage_estimators::*;
use nm;
use FlagFilter;
use ReadsMapped;

/// A single gene (or other feature) parsed from a GFF file.
#[derive(Clone, Debug, PartialEq)]
pub struct Gene {
    /// Identifier reported for the gene in the output.
    pub id: String,
    /// Name of the contig / reference sequence the gene resides on.
    pub contig: String,
    /// 0-based, inclusive start coordinate.
    pub start: u64,
    /// 0-based, exclusive end coordinate.
    pub end: u64,
}

/// A set of gene definitions, typically parsed from a GFF file.
#[derive(Clone, Debug)]
pub struct GeneDefinitions {
    pub genes: Vec<Gene>,
}

impl GeneDefinitions {
    /// Parse gene definitions from a GFF (or GTF) file. Each non-comment line is
    /// treated as a separate feature. When `feature_type` is `Some`, only lines
    /// whose third column matches are retained.
    pub fn read_gff(path: &str, feature_type: Option<&str>) -> GeneDefinitions {
        let file = std::fs::File::open(path)
            .unwrap_or_else(|e| panic!("Failed to open GFF file {}: {}", path, e));
        let reader = BufReader::new(file);
        let mut genes = vec![];
        let mut auto_id = 0u64;

        for (line_number, line_result) in reader.lines().enumerate() {
            let line =
                line_result.unwrap_or_else(|e| panic!("Failed to read GFF file {}: {}", path, e));
            let trimmed = line.trim_end();
            // Skip comments, directives (e.g. "##gff-version 3") and blank lines.
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = trimmed.split('\t').collect();
            if fields.len() < 8 {
                warn!(
                    "Skipping malformed GFF line {} in {} (expected at least 8 \
                     tab-separated fields, found {})",
                    line_number + 1,
                    path,
                    fields.len()
                );
                continue;
            }

            if let Some(wanted) = feature_type {
                if fields[2] != wanted {
                    continue;
                }
            }

            let contig = fields[0].to_string();
            let start_1based: u64 = match fields[3].parse() {
                Ok(s) => s,
                Err(_) => {
                    warn!(
                        "Skipping GFF line {} in {}: could not parse start coordinate '{}'",
                        line_number + 1,
                        path,
                        fields[3]
                    );
                    continue;
                }
            };
            let end_1based: u64 = match fields[4].parse() {
                Ok(e) => e,
                Err(_) => {
                    warn!(
                        "Skipping GFF line {} in {}: could not parse end coordinate '{}'",
                        line_number + 1,
                        path,
                        fields[4]
                    );
                    continue;
                }
            };
            if start_1based == 0 || end_1based < start_1based {
                warn!(
                    "Skipping GFF line {} in {}: invalid coordinates {}-{}",
                    line_number + 1,
                    path,
                    start_1based,
                    end_1based
                );
                continue;
            }

            let attributes = fields.get(8).copied().unwrap_or("");
            let id = parse_gff_id(attributes).unwrap_or_else(|| {
                auto_id += 1;
                format!("{contig}_gene_{auto_id}")
            });

            genes.push(Gene {
                id,
                contig,
                // Convert from 1-based inclusive (GFF) to 0-based half-open.
                start: start_1based - 1,
                end: end_1based,
            });
        }

        info!("Read in {} gene definitions from {}", genes.len(), path);
        GeneDefinitions { genes }
    }
}

/// Extract a human-readable identifier from a GFF/GTF attributes column,
/// trying a number of common attribute keys in turn.
fn parse_gff_id(attributes: &str) -> Option<String> {
    for key in &["ID", "locus_tag", "gene_id", "Name", "gene", "Parent"] {
        if let Some(value) = parse_gff_attribute(attributes, key) {
            if !value.is_empty() {
                return Some(value);
            }
        }
    }
    None
}

/// Find the value of a single attribute, supporting both GFF3 (`key=value`) and
/// GTF (`key "value"`) styles.
fn parse_gff_attribute(attributes: &str, key: &str) -> Option<String> {
    for entry in attributes.split(';') {
        let entry = entry.trim();
        if entry.is_empty() {
            continue;
        }
        // GFF3 style: key=value
        if let Some(rest) = entry.strip_prefix(&format!("{key}=")) {
            return Some(rest.trim().to_string());
        }
        // GTF style: key "value"
        if let Some(rest) = entry.strip_prefix(&format!("{key} ")) {
            return Some(rest.trim().trim_matches('"').to_string());
        }
    }
    None
}

/// A gene resolved against a particular BAM header, ready for coverage
/// calculation. Coordinates are clamped to the length of the contig.
struct ResolvedGene {
    entry_id: usize,
    name: String,
    start: usize,
    end: usize,
}

/// Calculate per-gene coverage for each BAM, using the coordinates defined in
/// `gene_definitions`. Reads are assigned to a gene (for read-count based
/// methods) when their leftmost mapped position falls within the gene's
/// coordinates.
pub fn gene_coverage<R: NamedBamReader, G: NamedBamReaderGenerator<R>, T: CoverageTaker>(
    bam_readers: Vec<G>,
    coverage_taker: &mut T,
    coverage_estimators: &mut [CoverageEstimator],
    gene_definitions: &GeneDefinitions,
    print_zero_coverage_genes: bool,
    flag_filters: &FlagFilter,
    threads: u16,
) -> Vec<ReadsMapped> {
    let mut reads_mapped_vector = vec![];

    for bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();
        bam_generated.set_threads(threads as usize);

        let stoit_name = bam_generated.name().to_string();
        coverage_taker.start_stoit(&stoit_name);

        let header = bam_generated.header().clone();
        let genes_by_tid = resolve_genes_against_header(gene_definitions, &header);

        let mut record: bam::record::Record = bam::record::Record::new();
        let mut last_tid: i32 = -2; // no such tid in a real BAM file
        let mut ups_and_downs: Vec<i32> = Vec::new();
        // Per-contig record of mapped reads, in ascending start-position order.
        let mut read_starts: Vec<u64> = vec![];
        let mut read_is_primary: Vec<u64> = vec![];
        let mut read_mismatches: Vec<u64> = vec![];
        let mut read_identities: Vec<f64> = vec![];
        let mut num_mapped_reads_total: u64 = 0;

        loop {
            match bam_generated.read(&mut record) {
                None => break,
                Some(Ok(())) => {}
                Some(e) => panic!("Error reading BAM record: {:?}", e),
            }
            if !flag_filters.passes(&record) {
                continue;
            }
            if record.is_unmapped() {
                continue;
            }

            let tid = record.tid();
            if tid != last_tid {
                if tid < last_tid {
                    error!(
                        "BAM file appears to be unsorted. Input BAM files must be \
                         sorted by reference (i.e. by samtools sort)"
                    );
                    panic!("BAM file appears to be unsorted.");
                }
                process_previous_genes(
                    last_tid,
                    tid,
                    &ups_and_downs,
                    &read_starts,
                    &read_is_primary,
                    &read_mismatches,
                    &read_identities,
                    &genes_by_tid,
                    coverage_estimators,
                    coverage_taker,
                    print_zero_coverage_genes,
                );
                ups_and_downs =
                    vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                last_tid = tid;
                read_starts.clear();
                read_is_primary.clear();
                read_mismatches.clear();
                read_identities.clear();
            }

            let is_primary = !record.is_supplementary() && !record.is_secondary();
            if is_primary {
                num_mapped_reads_total += 1;
            }

            let mut cursor: usize = record.pos() as usize;
            let mut aligned_len: u64 = 0;
            let mut indels: u64 = 0;
            let contig_len = ups_and_downs.len();
            for cig in record.cigar().iter() {
                match cig {
                    Cigar::Match(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                        ups_and_downs[cursor] += 1;
                        let final_pos = cursor + cig.len() as usize;
                        if final_pos < contig_len {
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

            let edit = nm(&record);
            read_starts.push(record.pos() as u64);
            read_is_primary.push(if is_primary { 1 } else { 0 });
            // Number of substitutions (NM minus indel length).
            read_mismatches.push(edit.saturating_sub(indels));
            read_identities.push(if is_primary && aligned_len > 0 {
                (aligned_len as f64 - edit as f64) / aligned_len as f64
            } else {
                0.0
            });
        }

        process_previous_genes(
            last_tid,
            genes_by_tid.len() as i32,
            &ups_and_downs,
            &read_starts,
            &read_is_primary,
            &read_mismatches,
            &read_identities,
            &genes_by_tid,
            coverage_estimators,
            coverage_taker,
            print_zero_coverage_genes,
        );

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

/// Map each gene onto a tid in the BAM header, dropping genes whose contig is
/// absent. Returns a per-tid list of genes, each sorted by start coordinate,
/// with a globally-consistent entry id assigned in (tid, start) order.
fn resolve_genes_against_header(
    gene_definitions: &GeneDefinitions,
    header: &bam::HeaderView,
) -> Vec<Vec<ResolvedGene>> {
    let target_names = header.target_names();
    let mut name_to_tid: HashMap<&[u8], usize> = HashMap::new();
    for (tid, name) in target_names.iter().enumerate() {
        name_to_tid.insert(name, tid);
    }

    let mut genes_by_tid: Vec<Vec<ResolvedGene>> =
        (0..target_names.len()).map(|_| vec![]).collect();
    let mut num_skipped = 0;
    for gene in &gene_definitions.genes {
        match name_to_tid.get(gene.contig.as_bytes()) {
            Some(&tid) => {
                let contig_len = header.target_len(tid as u32).expect("Corrupt BAM file?");
                let start = gene.start.min(contig_len);
                let end = gene.end.min(contig_len);
                if start >= end {
                    num_skipped += 1;
                    continue;
                }
                genes_by_tid[tid].push(ResolvedGene {
                    entry_id: 0, // assigned below
                    name: gene.id.clone(),
                    start: start as usize,
                    end: end as usize,
                });
            }
            None => num_skipped += 1,
        }
    }

    if num_skipped > 0 {
        warn!(
            "{num_skipped} gene(s) were ignored because their contig was not \
             present in the reference, or they had invalid coordinates"
        );
    }

    let mut next_entry_id = 0;
    for genes in genes_by_tid.iter_mut() {
        genes.sort_by_key(|g| g.start);
        for gene in genes.iter_mut() {
            gene.entry_id = next_entry_id;
            next_entry_id += 1;
        }
    }
    genes_by_tid
}

#[allow(clippy::too_many_arguments)]
fn process_previous_genes<T: CoverageTaker>(
    last_tid: i32,
    current_tid: i32,
    ups_and_downs: &[i32],
    read_starts: &[u64],
    read_is_primary: &[u64],
    read_mismatches: &[u64],
    read_identities: &[f64],
    genes_by_tid: &[Vec<ResolvedGene>],
    coverage_estimators: &mut [CoverageEstimator],
    coverage_taker: &mut T,
    print_zero_coverage_genes: bool,
) {
    // Emit coverage for the contig that was just finished.
    if last_tid != -2 {
        emit_genes_for_contig(
            &genes_by_tid[last_tid as usize],
            ups_and_downs,
            read_starts,
            read_is_primary,
            read_mismatches,
            read_identities,
            coverage_estimators,
            coverage_taker,
            print_zero_coverage_genes,
        );
    }

    // Emit zero coverage for genes on contigs that had no mapped reads at all.
    if print_zero_coverage_genes {
        let mut my_tid = if last_tid == -2 { 0 } else { last_tid + 1 };
        while my_tid < current_tid {
            emit_zero_coverage_genes(
                &genes_by_tid[my_tid as usize],
                coverage_estimators,
                coverage_taker,
            );
            my_tid += 1;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn emit_genes_for_contig<T: CoverageTaker>(
    genes: &[ResolvedGene],
    ups_and_downs: &[i32],
    read_starts: &[u64],
    read_is_primary: &[u64],
    read_mismatches: &[u64],
    read_identities: &[f64],
    coverage_estimators: &mut [CoverageEstimator],
    coverage_taker: &mut T,
    print_zero_coverage_genes: bool,
) {
    if genes.is_empty() {
        return;
    }

    // Cumulative coverage at each base of the contig.
    let contig_len = ups_and_downs.len();
    let mut coverage_at_base = vec![0i32; contig_len];
    let mut running = 0i32;
    for (i, ud) in ups_and_downs.iter().enumerate() {
        running += ud;
        coverage_at_base[i] = running;
    }

    // Prefix sums over the per-read contributions, so that the values for an
    // arbitrary span of read start positions can be computed in O(log n).
    let n = read_starts.len();
    let mut prefix_primary = vec![0u64; n + 1];
    let mut prefix_mismatches = vec![0u64; n + 1];
    let mut prefix_identity = vec![0f64; n + 1];
    for i in 0..n {
        prefix_primary[i + 1] = prefix_primary[i] + read_is_primary[i];
        prefix_mismatches[i + 1] = prefix_mismatches[i] + read_mismatches[i];
        prefix_identity[i + 1] = prefix_identity[i] + read_identities[i];
    }

    for gene in genes {
        let start = gene.start;
        let end = gene.end.min(contig_len);
        if start >= end {
            continue;
        }
        let len = end - start;

        // Build a delta-encoded coverage array for the gene region whose
        // cumulative sums reproduce the true coverage within the gene.
        let mut gene_ups_and_downs = vec![0i32; len];
        gene_ups_and_downs[0] = coverage_at_base[start];
        if len > 1 {
            gene_ups_and_downs[1..len].copy_from_slice(&ups_and_downs[(start + 1)..end]);
        }

        // Reads assigned to this gene by leftmost mapped position.
        let lo = read_starts.partition_point(|&x| (x as usize) < start);
        let hi = read_starts.partition_point(|&x| (x as usize) < end);
        let num_mapped_reads = prefix_primary[hi] - prefix_primary[lo];
        let mismatches = prefix_mismatches[hi] - prefix_mismatches[lo];
        let sum_identity = prefix_identity[hi] - prefix_identity[lo];

        for estimator in coverage_estimators.iter_mut() {
            estimator.add_contig(
                &gene_ups_and_downs,
                num_mapped_reads,
                mismatches,
                sum_identity,
            );
        }
        let coverages: Vec<f32> = coverage_estimators
            .iter_mut()
            .map(|estimator| estimator.calculate_coverage(&[0]))
            .collect();
        let has_nonzero_coverage = coverages.iter().any(|&c| c > 0.0);

        if print_zero_coverage_genes || has_nonzero_coverage {
            coverage_taker.start_entry(gene.entry_id, &gene.name);
            for (coverage, estimator) in coverages.into_iter().zip(coverage_estimators.iter_mut()) {
                estimator.print_coverage(coverage, coverage_taker);
            }
            coverage_taker.finish_entry();
        }

        for estimator in coverage_estimators.iter_mut() {
            estimator.setup();
        }
    }
}

fn emit_zero_coverage_genes<T: CoverageTaker>(
    genes: &[ResolvedGene],
    coverage_estimators: &[CoverageEstimator],
    coverage_taker: &mut T,
) {
    for gene in genes {
        coverage_taker.start_entry(gene.entry_id, &gene.name);
        let gene_len = (gene.end - gene.start) as u64;
        for estimator in coverage_estimators.iter() {
            estimator.print_zero_coverage(coverage_taker, gene_len);
        }
        coverage_taker.finish_entry();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use OutputWriter;

    fn run_genes(
        gene_definitions: &GeneDefinitions,
        bam_files: Vec<&str>,
        coverage_estimators: &mut [CoverageEstimator],
        print_zeros: bool,
    ) -> String {
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let flag_filters = FlagFilter {
            include_improper_pairs: true,
            include_secondary: false,
            include_supplementary: false,
        };
        {
            let mut coverage_taker =
                CoverageTakerType::new_single_float_coverage_streaming_coverage_printer(
                    OutputWriter::generate(Some(tf.path().to_str().unwrap())),
                );
            gene_coverage(
                generate_named_bam_readers_from_bam_files(bam_files),
                &mut coverage_taker,
                coverage_estimators,
                gene_definitions,
                print_zeros,
                &flag_filters,
                1,
            );
        }
        std::fs::read_to_string(tf.path()).unwrap()
    }

    #[test]
    fn test_gff_parsing_gff3() {
        let defs = GeneDefinitions::read_gff("tests/data/2seqs.gff", None);
        assert_eq!(
            defs.genes,
            vec![
                Gene {
                    id: "gene1".to_string(),
                    contig: "seq1".to_string(),
                    start: 0,
                    end: 1000,
                },
                Gene {
                    id: "gene2".to_string(),
                    contig: "seq1".to_string(),
                    start: 99,
                    end: 200,
                },
                Gene {
                    id: "gene3".to_string(),
                    contig: "seq2".to_string(),
                    start: 0,
                    end: 1000,
                },
            ]
        );
    }

    #[test]
    fn test_whole_contig_gene_matches_contig_coverage() {
        // A gene spanning an entire contig must give the same mean coverage as
        // running coverm in contig mode (with no contig-end-exclusion).
        let defs = GeneDefinitions {
            genes: vec![
                Gene {
                    id: "gene_seq1".to_string(),
                    contig: "seq1".to_string(),
                    start: 0,
                    end: 1000,
                },
                Gene {
                    id: "gene_seq2".to_string(),
                    contig: "seq2".to_string(),
                    start: 0,
                    end: 1000,
                },
            ],
        };
        let observed = run_genes(
            &defs,
            vec!["tests/data/2seqs.reads_for_seq1.bam"],
            &mut [CoverageEstimator::new_estimator_mean(0.0, 0, false)],
            true,
        );
        assert_eq!(
            "2seqs.reads_for_seq1\tgene_seq1\t1.2\n\
             2seqs.reads_for_seq1\tgene_seq2\t0\n",
            observed
        );
    }

    #[test]
    fn test_no_zeros_omits_uncovered_genes() {
        let defs = GeneDefinitions {
            genes: vec![
                Gene {
                    id: "gene_seq1".to_string(),
                    contig: "seq1".to_string(),
                    start: 0,
                    end: 1000,
                },
                Gene {
                    id: "gene_seq2".to_string(),
                    contig: "seq2".to_string(),
                    start: 0,
                    end: 1000,
                },
            ],
        };
        let observed = run_genes(
            &defs,
            vec!["tests/data/2seqs.reads_for_seq1.bam"],
            &mut [CoverageEstimator::new_estimator_mean(0.0, 0, false)],
            false,
        );
        assert_eq!("2seqs.reads_for_seq1\tgene_seq1\t1.2\n", observed);
    }

    #[test]
    fn test_count_method_per_gene() {
        // The whole-contig gene should accumulate every read mapped to seq1.
        let defs = GeneDefinitions {
            genes: vec![Gene {
                id: "gene_seq1".to_string(),
                contig: "seq1".to_string(),
                start: 0,
                end: 1000,
            }],
        };
        let observed = run_genes(
            &defs,
            vec!["tests/data/2seqs.reads_for_seq1.bam"],
            &mut [CoverageEstimator::new_estimator_read_count()],
            false,
        );
        assert_eq!("2seqs.reads_for_seq1\tgene_seq1\t12\n", observed);
    }
}
