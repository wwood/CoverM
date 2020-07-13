use std;
use std::process;

use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use bam_generator::*;
use coverage_takers::*;
use mosdepth_genome_coverage_estimators::*;
use FlagFilter;
use ReadsMapped;

pub fn contig_coverage<R: NamedBamReader, G: NamedBamReaderGenerator<R>, T: CoverageTaker>(
    bam_readers: Vec<G>,
    coverage_taker: &mut T,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    print_zero_coverage_contigs: bool,
    flag_filters: FlagFilter,
    threads: usize,
) -> Vec<ReadsMapped> {
    let mut reads_mapped_vector = vec![];
    for bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();
        bam_generated.set_threads(threads);

        let stoit_name = &(bam_generated.name().to_string());
        coverage_taker.start_stoit(stoit_name);
        let mut record: bam::record::Record = bam::record::Record::new();
        let mut last_tid: i32 = -2; // no such tid in a real BAM file
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let header = bam_generated.header().clone();
        let target_names = header.target_names();

        let mut num_mapped_reads_total: u64 = 0;
        let mut num_mapped_reads_in_current_contig: u64 = 0;
        let mut total_indels_in_current_contig: u64 = 0;
        let mut total_edit_distance_in_current_contig: u64 = 0;

        let mut process_previous_contigs =
            |last_tid,
             tid,
             coverage_estimators: &mut Vec<CoverageEstimator>,
             ups_and_downs,
             num_mapped_reads_in_current_contig,
             total_edit_distance_in_current_contig,
             total_indels_in_current_contig,
             num_mapped_reads_total: &mut u64| {
                if last_tid != -2 {
                    debug!(
                        "Found {} reads mapped to tid {}, with total edit \
                        distance {} and {} indels",
                        num_mapped_reads_in_current_contig,
                        last_tid,
                        total_edit_distance_in_current_contig,
                        total_indels_in_current_contig
                    );
                    for estimator in coverage_estimators.iter_mut() {
                        estimator.add_contig(
                            &ups_and_downs,
                            num_mapped_reads_in_current_contig,
                            total_edit_distance_in_current_contig - total_indels_in_current_contig,
                        )
                    }
                    let coverages: Vec<f32> = coverage_estimators
                        .iter_mut()
                        .map(|estimator| estimator.calculate_coverage(&vec![0]))
                        .collect();
                    let has_nonzero_coverage = coverages.iter().any(|&coverage| coverage > 0.0);
                    debug!("Found nonzero coverage?: {}", has_nonzero_coverage);
                    if has_nonzero_coverage {
                        *num_mapped_reads_total += num_mapped_reads_in_current_contig;
                    }
                    if print_zero_coverage_contigs || has_nonzero_coverage {
                        coverage_taker.start_entry(
                            last_tid as usize,
                            std::str::from_utf8(target_names[last_tid as usize]).unwrap(),
                        );
                        for (coverage, estimator) in
                            coverages.iter().zip(coverage_estimators.iter_mut())
                        {
                            estimator.print_coverage(&coverage, coverage_taker);
                            estimator.setup();
                        }
                        coverage_taker.finish_entry();
                    }
                }
                if print_zero_coverage_contigs {
                    print_previous_zero_coverage_contigs(
                        match last_tid {
                            -2 => -1,
                            _ => last_tid,
                        },
                        tid,
                        coverage_estimators,
                        &target_names,
                        coverage_taker,
                        &header,
                    );
                }
            };

        // for record in records
        while bam_generated
            .read(&mut record)
            .expect("Error while reading BAM record")
            == true
        {
            trace!("Starting with a new read.. {:?}", record);
            if (!flag_filters.include_supplementary && record.is_supplementary())
                || (!flag_filters.include_secondary && record.is_secondary())
                || (!flag_filters.include_improper_pairs && !record.is_proper_pair())
            {
                trace!("Skipping read based on flag filtering");
                continue;
            }
            // if reference has changed, print the last record
            let tid = record.tid();
            if !record.is_unmapped() {
                // if mapped
                // if reference has changed, print the last record
                if tid != last_tid {
                    if tid < last_tid {
                        error!("BAM file appears to be unsorted. Input BAM files must be sorted by reference (i.e. by samtools sort)");
                        panic!("BAM file appears to be unsorted. Input BAM files must be sorted by reference (i.e. by samtools sort)");
                    }
                    process_previous_contigs(
                        last_tid,
                        tid,
                        coverage_estimators,
                        ups_and_downs,
                        num_mapped_reads_in_current_contig,
                        total_edit_distance_in_current_contig,
                        total_indels_in_current_contig,
                        &mut num_mapped_reads_total,
                    );
                    ups_and_downs =
                        vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    debug!(
                        "Working on new reference {}",
                        std::str::from_utf8(target_names[tid as usize]).unwrap()
                    );
                    last_tid = tid;
                    num_mapped_reads_in_current_contig = 0;
                    total_edit_distance_in_current_contig = 0;
                    total_indels_in_current_contig = 0;
                }

                if !record.is_supplementary() && !record.is_secondary() {
                    num_mapped_reads_in_current_contig += 1;
                }

                // for each chunk of the cigar string
                trace!(
                    "read name {:?}",
                    std::str::from_utf8(record.qname()).unwrap()
                );
                let mut cursor: usize = record.pos() as usize;
                for cig in record.cigar().iter() {
                    trace!("Found cigar {:} from {}", cig, cursor);
                    match cig {
                        Cigar::Match(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                            // if M, X, or = increment start and decrement end index
                            trace!(
                                "Adding M, X, or = at {} and {}",
                                cursor,
                                cursor + cig.len() as usize
                            );
                            ups_and_downs[cursor] += 1;
                            let final_pos = cursor + cig.len() as usize;
                            if final_pos < ups_and_downs.len() {
                                // True unless the read hits the contig end.
                                ups_and_downs[final_pos] -= 1;
                            }
                            cursor += cig.len() as usize;
                        }
                        Cigar::Del(_) => {
                            cursor += cig.len() as usize;
                            total_indels_in_current_contig += cig.len() as u64;
                        }
                        Cigar::RefSkip(_) => {
                            // if D or N, move the cursor
                            cursor += cig.len() as usize;
                        }
                        Cigar::Ins(_) => {
                            total_indels_in_current_contig += cig.len() as u64;
                        }
                        Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
                    }
                }

                // Determine the number of mismatching bases in this read by
                // looking at the NM tag.
                total_edit_distance_in_current_contig += match record.aux("NM".as_bytes()) {
                    Some(aux) => aux.integer() as u64,
                    None => {
                        error!(
                            "Mapping record encountered that does not have an 'NM' \
                                    auxiliary tag in the SAM/BAM format. This is required \
                                    to work out some coverage statistics"
                        );
                        process::exit(1);
                    }
                };

                trace!("At end of loop")
            }
        }

        process_previous_contigs(
            last_tid,
            target_names.len() as i32,
            coverage_estimators,
            ups_and_downs,
            num_mapped_reads_in_current_contig,
            total_edit_distance_in_current_contig,
            total_indels_in_current_contig,
            &mut num_mapped_reads_total,
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
                "No primary alignments were observed for sample {} \
                   - perhaps something went wrong in the mapping?",
                stoit_name
            );
        }

        bam_generated.finish();
    }
    return reads_mapped_vector;
}

fn print_previous_zero_coverage_contigs<T: CoverageTaker>(
    last_tid: i32,
    current_tid: i32,
    coverage_estimators: &Vec<CoverageEstimator>,
    target_names: &Vec<&[u8]>,
    coverage_taker: &mut T,
    header: &bam::HeaderView,
) {
    let mut my_tid = last_tid + 1;
    while my_tid < current_tid {
        debug!("printing zero coverage for tid {}", my_tid);
        coverage_taker.start_entry(
            my_tid as usize,
            std::str::from_utf8(target_names[my_tid as usize]).unwrap(),
        );
        for ref coverage_estimator in coverage_estimators.iter() {
            coverage_estimator
                .print_zero_coverage(coverage_taker, header.target_len(my_tid as u32).unwrap());
        }
        coverage_taker.finish_entry();
        my_tid += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use genome_exclusion::*;
    use mapping_parameters::*;
    use shard_bam_reader::*;
    use std::io::Cursor;
    use std::str;

    fn test_with_stream<R: NamedBamReader, G: NamedBamReaderGenerator<R>>(
        expected: &str,
        bam_readers: Vec<G>,
        coverage_estimators: &mut Vec<CoverageEstimator>,
        print_zero_coverage_contigs: bool,
        proper_pairs_only: bool,
    ) -> Vec<ReadsMapped> {
        let mut stream = Cursor::new(Vec::new());
        let flag_filters = FlagFilter {
            include_improper_pairs: !proper_pairs_only,
            include_secondary: false,
            include_supplementary: false,
        };
        let reads_mapped_vec;
        {
            let mut coverage_taker =
                CoverageTakerType::new_single_float_coverage_streaming_coverage_printer(
                    &mut stream,
                );
            reads_mapped_vec = contig_coverage(
                bam_readers,
                &mut coverage_taker,
                coverage_estimators,
                print_zero_coverage_contigs,
                flag_filters,
                1,
            );
        }
        assert_eq!(expected, str::from_utf8(stream.get_ref()).unwrap());
        return reads_mapped_vec;
    }

    #[test]
    fn test_one_genome_two_contigs_first_covered_no_zeros() {
        test_with_stream(
            "7seqs.reads_for_seq1_and_seq2\tgenome2~seq1\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome5~seq2\t1.2\n",
            generate_named_bam_readers_from_bam_files(
                vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            &mut vec!(CoverageEstimator::new_estimator_mean(0.0,0,false)),
            false,
            false);
    }

    #[test]
    fn test_one_genome_two_contigs_first_covered() {
        test_with_stream(
            "7seqs.reads_for_seq1_and_seq2\tgenome1~random_sequence_length_11000\t0\n7seqs.reads_for_seq1_and_seq2\tgenome1~random_sequence_length_11010\t0\n7seqs.reads_for_seq1_and_seq2\tgenome2~seq1\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome3~random_sequence_length_11001\t0\n7seqs.reads_for_seq1_and_seq2\tgenome4~random_sequence_length_11002\t0\n7seqs.reads_for_seq1_and_seq2\tgenome5~seq2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6~random_sequence_length_11003\t0\n",
            generate_named_bam_readers_from_bam_files(
                vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            &mut vec!(CoverageEstimator::new_estimator_mean(0.0,0,false)),
            true,
            false);
    }

    #[test]
    fn test_sharded_bams_some_not_covered() {
        test_with_stream(
            "shard1|shard2\tgenome3~random_sequence_length_11001\t0.10908099\
            \nshard1|shard2\tgenome4~random_sequence_length_11002\t0.109071076\
            \nshard1|shard2\tgenome5~seq2\t0\
            \nshard1|shard2\tgenome6~random_sequence_length_11003\t0.10906117\
            \nshard1|shard2\tgenome1~random_sequence_length_11000\t0.10909091\
            \nshard1|shard2\tgenome1~random_sequence_length_11010\t0.108991824\
            \nshard1|shard2\tgenome2~seq1\t0\n",
            generate_sharded_bam_reader_from_bam_files(
                vec!["tests/data/shard1.bam", "tests/data/shard2.bam"],
                4,
                &NoExclusionGenomeFilter {},
            ),
            &mut vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
            true,
            false,
        );
    }

    #[test]
    fn test_proper_pairs_only() {
        test_with_stream(
            "",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/1.bam"]),
            &mut vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
            false,
            true,
        );
    }

    #[test]
    fn test_one_contig_variance() {
        test_with_stream(
            &("2seqs.reads_for_seq1\tseq1\t0.9489489\n".to_owned()
                + "2seqs.reads_for_seq1\tseq2\t0\n"),
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            &mut vec![CoverageEstimator::new_estimator_variance(0.0, 0)],
            true,
            false,
        );
    }

    #[test]
    fn test_streaming_bam_file() {
        test_with_stream(
            "7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome2~seq1\t1.2\n\
             7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome5~seq2\t1.2\n",
            vec![generate_named_bam_readers_from_reads(
                MappingProgram::BWA_MEM,
                "tests/data/7seqs.fna",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                Some("tests/data/reads_for_seq1_and_seq2.2.fq.gz"),
                ReadFormat::Coupled,
                4,
                None,
                false,
                None,
                true,
            )],
            &mut vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
            false,
            false,
        );
    }

    #[test]
    fn test_multiple_coverage_methods() {
        test_with_stream(
            &("2seqs.reads_for_seq1\tseq1\t1.2\t0.9489489\n".to_owned()
                + "2seqs.reads_for_seq1\tseq2\t0\t0\n"),
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            &mut vec![
                CoverageEstimator::new_estimator_mean(0.0, 0, false),
                CoverageEstimator::new_estimator_variance(0.0, 0),
            ],
            true,
            false,
        );
    }

    #[test]
    fn test_julian_error() {
        test_with_stream(
            "2seqs.reads_for_seq1.with_unmapped\tseq1\t1.497\n\
             2seqs.reads_for_seq1.with_unmapped\tseq2\t1.5\n",
            // has unmapped reads, which caused problems with --proper-pairs-only.
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/2seqs.reads_for_seq1.with_unmapped.bam",
            ]),
            &mut vec![CoverageEstimator::new_estimator_mean(0.0, 0, true)],
            true,
            false,
        );
    }

    #[test]
    fn test_trimmed_mean_bug() {
        test_with_stream(
            &("2seqs.reads_for_seq1\tseq1\t0\n".to_owned() + "2seqs.reads_for_seq1\tseq2\t0\n"),
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            &mut vec![CoverageEstimator::new_estimator_trimmed_mean(
                0.0, 0.05, 0.0, 0,
            )],
            true,
            false,
        );
    }

    #[test]
    fn test_multiple_outputs_one_zero_no_print_zeroes() {
        test_with_stream(
            "2seqs.reads_for_seq1\tseq1\t1.2\t0\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            &mut vec![
                CoverageEstimator::new_estimator_mean(0.0, 0, false),
                // covered fraction is 0.727, so go lower so trimmed mean is 0,
                // mean > 0.
                CoverageEstimator::new_estimator_trimmed_mean(0.0, 0.05, 0.0, 0),
            ],
            false,
            false,
        );
    }

    #[test]
    fn test_multiple_outputs_one_zero_no_print_zeroes_reverse_order() {
        test_with_stream(
            "2seqs.reads_for_seq1\tseq1\t0\t1.2\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            &mut vec![
                // covered fraction is 0.727, so go lower so trimmed mean is 0,
                // mean > 0.
                CoverageEstimator::new_estimator_trimmed_mean(0.0, 0.05, 0.0, 0),
                CoverageEstimator::new_estimator_mean(0.0, 0, false),
            ],
            false,
            false,
        );
    }

    #[test]
    fn test_contig_end_exclusion() {
        // From https://bitbucket.org/berkeleylab/metabat/issues/48/jgi_summarize_bam_contig_depths-coverage
        test_with_stream(
            &("7seqs.reads_for_seq1_and_seq2\tgenome2~seq1\t1.4117647\t1.3049262\n".to_owned()
                + "7seqs.reads_for_seq1_and_seq2\tgenome5~seq2\t1.2435294\t0.6862065\n"),
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
            ]),
            &mut vec![
                CoverageEstimator::new_estimator_mean(0.0, 75, false),
                // covered fraction is 0.727, so go lower so trimmed mean is 0,
                // mean > 0.
                CoverageEstimator::new_estimator_variance(0.0, 75),
            ],
            false,
            false,
        );
    }

    #[test]
    fn test_one_read_of_pair_mapped() {
        // In the second read, tid is != -1, because the first in the pair is assigned somewhere.
        test_with_stream(
            &("1read_of_pair_mapped\t73.20100900_E1D.16_contig_9606\t0.011293635\n".to_owned()),
            generate_named_bam_readers_from_bam_files(vec!["tests/data/1read_of_pair_mapped.bam"]),
            &mut vec![CoverageEstimator::new_estimator_mean(0.0, 75, true)],
            false,
            false,
        );
    }

    #[test]
    fn test_variance_estimator_all_bases_covered() {
        // In the past this threw up a underflow error
        test_with_stream(
            &("k141_2005182	k141_2005182	5.107387\n".to_owned()),
            generate_named_bam_readers_from_bam_files(vec!["tests/data/k141_2005182.bam"]),
            &mut vec![CoverageEstimator::new_estimator_variance(0.0, 75)],
            false,
            false,
        );
    }

    #[test]
    fn test_reads_not_counting_when_sufficient_min_covered() {
        // In the past this threw up a underflow error
        let reads_mapped = test_with_stream(
            &("2seqs.reads_for_seq1_and_seq2	seq1	1.3049262\n\
               2seqs.reads_for_seq1_and_seq2	seq2	0.6862065\n"
                .to_owned()),
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/2seqs.reads_for_seq1_and_seq2.bam",
            ]),
            &mut vec![CoverageEstimator::new_estimator_variance(0.0, 75)],
            false,
            false,
        );
        assert_eq!(
            vec!(ReadsMapped {
                num_mapped_reads: 24,
                num_reads: 24
            }),
            reads_mapped
        );
    }

    #[test]
    fn test_reads_not_counting_when_insufficient_min_covered() {
        // In the past this threw up a underflow error
        let reads_mapped = test_with_stream(
            &("".to_owned()),
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/2seqs.reads_for_seq1_and_seq2.bam",
            ]),
            &mut vec![CoverageEstimator::new_estimator_variance(0.99, 75)],
            false,
            false,
        );
        assert_eq!(
            vec!(ReadsMapped {
                num_mapped_reads: 0,
                num_reads: 24
            }),
            reads_mapped
        );
    }

    #[test]
    fn test_reads_per_base_estimator() {
        test_with_stream(
            "7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11000\t0\n\
             7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11010\t0\n\
             7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome2~seq1\t0.012\n\
             7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome3~random_sequence_length_11001\t0\n\
             7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome4~random_sequence_length_11002\t0\n\
             7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome5~seq2\t0.012\n\
             7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome6~random_sequence_length_11003\t0\n",
            vec![generate_named_bam_readers_from_reads(
                MappingProgram::BWA_MEM,
                "tests/data/7seqs.fna",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                Some("tests/data/reads_for_seq1_and_seq2.2.fq.gz"),
                ReadFormat::Coupled,
                4,
                None,
                false,
                None,
                true,
            )],
            &mut vec![CoverageEstimator::new_estimator_reads_per_base()],
            true,
            false,
        );
    }
}
