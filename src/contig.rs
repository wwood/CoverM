
use std;
use std::result::Result;
use std::any::Any;
use rust_htslib::bam;
use std::process;
use mosdepth_genome_coverage_estimators::*;
use bam_generator::*;
use get_trimmed_mean_estimator;

pub fn method_match(method: &str, min_fraction_covered: f32, m: &clap::ArgMatches)-> CoverageEstimatorMethods{
    let mut coverage_method;
    match method {
        "mean" =>
            coverage_method = CoverageEstimatorMethods::MeanGenomeCoverageEstimator(
                MeanGenomeCoverageEstimator::new(min_fraction_covered)),
        "coverage_histogram" =>
            coverage_method = CoverageEstimatorMethods::PileupCountsGenomeCoverageEstimator(
                PileupCountsGenomeCoverageEstimator::new(min_fraction_covered)),
        "trimmed_mean" =>
            coverage_method = get_trimmed_mean_estimator(m, min_fraction_covered),
        "covered_fraction" =>
            coverage_method = CoverageEstimatorMethods::CoverageFractionGenomeCoverageEstimator(
                CoverageFractionGenomeCoverageEstimator::new(min_fraction_covered)),
        "variance" =>
            coverage_method = CoverageEstimatorMethods::VarianceGenomeCoverageEstimator(
                VarianceGenomeCoverageEstimator::new(min_fraction_covered)),
        _ => panic!("programming error")
    }
    return coverage_method;
}

pub fn contig_coverage<T: MosdepthGenomeCoverageEstimator<T>,
                       R: NamedBamReader,
                       G: NamedBamReaderGenerator<R>>(
    bam_readers: Vec<G>,
    print_stream: &mut std::io::Write,
    min_fraction_covered: f32,
    methods: Vec<&str>,
    m: &clap::ArgMatches,
    print_zero_coverage_contigs: bool,
    flag_filtering: bool) {
//    let mut coverage_estimator_box: Box<Any>;
    for mut bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();

        let stoit_name = &(bam_generated.name().to_string());
        let mut record: bam::record::Record = bam::record::Record::new();
        let mut last_tid: i32 = -1; // no such tid in a real BAM file
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let header = bam_generated.header().clone();
        let target_names = header.target_names();

        // for record in records
        while bam_generated.read(&mut record).is_ok() {
            debug!("Starting with a new read.. {:?}", record);
            if flag_filtering &&
                (record.is_secondary() ||
                    record.is_supplementary() ||
                    !record.is_proper_pair()) {
                continue;
            }
            // if reference has changed, print the last record
            let tid = record.tid();
            print_contig(stoit_name, std::str::from_utf8(target_names[last_tid as usize]).unwrap(), print_stream);
            for method in methods {
//                match method {
//                    "mean" =>
//                        coverage_estimator_box = Box::new(MeanGenomeCoverageEstimator::new(
//                            min_fraction_covered)),
//                    "coverage_histogram" =>
//                        coverage_estimator_box = Box::new(PileupCountsGenomeCoverageEstimator::new(
//                            min_fraction_covered)),
//                    "trimmed_mean" =>
//                        coverage_estimator_box = Box::new(get_trimmed_mean_estimator(m, min_fraction_covered)),
//                    "covered_fraction" =>
//                        coverage_estimator_box = Box::new(CoverageFractionGenomeCoverageEstimator::new(
//                            min_fraction_covered)),
//                    "variance" =>
//                        coverage_estimator_box= Box::new(VarianceGenomeCoverageEstimator::new(
//                            min_fraction_covered)),
//                    _ => panic!("programming error")
//                }

                let mut coverage_method = method_match(method, min_fraction_covered, m);
                if let CoverageEstimatorMethods::MeanGenomeCoverageEstimator(MeanGenomeCoverageEstimator) = coverage_method {
                    let mut coverage_estimator = MeanGenomeCoverageEstimator;
                    contig_print_reference(coverage_estimator,
                                           tid,
                                           last_tid,
                                           &target_names,
                                           header.target_len(tid as u32).expect("Corrupt BAM file?") as usize,
                                           ups_and_downs.clone(),
                                           print_zero_coverage_contigs,
                                           print_stream);
                }
                if let CoverageEstimatorMethods::PileupCountsGenomeCoverageEstimator(PileupCountsGenomeCoverageEstimator) = coverage_method {
                    let mut coverage_estimator = PileupCountsGenomeCoverageEstimator;
                    contig_print_reference(coverage_estimator,
                                           tid,
                                           last_tid,
                                           &target_names,
                                           header.target_len(tid as u32).expect("Corrupt BAM file?") as usize,
                                           ups_and_downs.clone(),
                                           print_zero_coverage_contigs,
                                           print_stream);
                }
                if let CoverageEstimatorMethods::CoverageFractionGenomeCoverageEstimator(CoverageFractionGenomeCoverageEstimator) = coverage_method {
                    let mut coverage_estimator = CoverageFractionGenomeCoverageEstimator;
                    contig_print_reference(coverage_estimator,
                                           tid,
                                           last_tid,
                                           &target_names,
                                           header.target_len(tid as u32).expect("Corrupt BAM file?") as usize,
                                           ups_and_downs.clone(),
                                           print_zero_coverage_contigs,
                                           print_stream);
                }
                if let CoverageEstimatorMethods::TrimmedMeanGenomeCoverageEstimator(TrimmedMeanGenomeCoverageEstimator) = coverage_method {
                    let mut coverage_estimator = TrimmedMeanGenomeCoverageEstimator;
                    contig_print_reference(coverage_estimator,
                                           tid,
                                           last_tid,
                                           &target_names,
                                           header.target_len(tid as u32).expect("Corrupt BAM file?") as usize,
                                           ups_and_downs.clone(),
                                           print_zero_coverage_contigs,
                                           print_stream);
                }
                if let CoverageEstimatorMethods::VarianceGenomeCoverageEstimator(VarianceGenomeCoverageEstimator) = coverage_method {
                    let mut coverage_estimator = VarianceGenomeCoverageEstimator;
                    contig_print_reference(coverage_estimator,
                                           tid,
                                           last_tid,
                                           &target_names,
                                           header.target_len(tid as u32).expect("Corrupt BAM file?") as usize,
                                           ups_and_downs.clone(),
                                           print_zero_coverage_contigs,
                                           print_stream);
                }
//                if tid != last_tid {
//                    if last_tid != -1 {
//                        coverage_estimator.add_contig(&ups_and_downs);
//                        let coverage = coverage_estimator.calculate_coverage(0);
//
//                        if coverage > 0.0 {
//                            coverage_estimator.print_coverage(
//                                &coverage,
//                                print_stream);
//                        } else if print_zero_coverage_contigs {
//                            coverage_estimator.print_zero_coverage(
//                                print_stream);
//                        }
//                    }
//                    // reset for next time
//                    coverage_estimator.setup();
//                    if print_zero_coverage_contigs {
//                        print_previous_zero_coverage_contigs(last_tid, tid, stoit_name, coverage_estimator, &target_names, print_stream);
//                    }
//                    ups_and_downs = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
//                    debug!("Working on new reference {}",
//                           std::str::from_utf8(target_names[tid as usize]).unwrap());
//                    last_tid = tid;
//                }
//
                // for each chunk of the cigar string
                debug!("read name {:?}", std::str::from_utf8(record.qname()).unwrap());
                let mut cursor: usize = record.pos() as usize;
                for cig in record.cigar().iter() {
                    debug!("Found cigar {:} from {}", cig, cursor);
                    match cig.char() {
                        'M' => {
                            // if M, increment start and decrement end index
                            debug!("Adding M at {} and {}", cursor, cursor + cig.len() as usize);
                            ups_and_downs[cursor] += 1;
                            let final_pos = cursor + cig.len() as usize;
                            if final_pos < ups_and_downs.len() { // True unless the read hits the contig end.
                                ups_and_downs[final_pos] -= 1;
                            }
                            cursor += cig.len() as usize;
                        },
                        'D' => {
                            // if D, move the cursor
                            cursor += cig.len() as usize;
                        },
                        '=' => panic!("CIGAR '=' detected, but this case is not correctly handled for now"),
                        _ => {}
                    }
                }
                debug!("At end of loop")
            }
        }
        for method in methods {
//            match method {
//                "mean" =>
//                    coverage_estimator_box = Box::new(MeanGenomeCoverageEstimator::new(
//                        min_fraction_covered)),
//                "coverage_histogram" =>
//                    coverage_estimator_box = Box::new(PileupCountsGenomeCoverageEstimator::new(
//                        min_fraction_covered)),
//                "trimmed_mean" =>
//                    coverage_estimator_box = Box::new(get_trimmed_mean_estimator(m, min_fraction_covered)),
//                "covered_fraction" =>
//                    coverage_estimator_box = Box::new(CoverageFractionGenomeCoverageEstimator::new(
//                        min_fraction_covered)),
//                "variance" =>
//                    coverage_estimator_box= Box::new(VarianceGenomeCoverageEstimator::new(
//                        min_fraction_covered)),
//                _ => panic!("programming error")
//            }

            let mut coverage_method = method_match(method, min_fraction_covered, m);
            if let CoverageEstimatorMethods::MeanGenomeCoverageEstimator(MeanGenomeCoverageEstimator) = coverage_method {
                let mut coverage_estimator = MeanGenomeCoverageEstimator;
                contig_print_last_reference(coverage_estimator,
                                            target_names,
                                            last_tid,
                                            ups_and_downs,
                                            print_zero_coverage_contigs,
                                            print_stream);
            }
            if let CoverageEstimatorMethods::PileupCountsGenomeCoverageEstimator(PileupCountsGenomeCoverageEstimator) = coverage_method {
                let mut coverage_estimator = PileupCountsGenomeCoverageEstimator;
                contig_print_last_reference(coverage_estimator,
                                            target_names,
                                            last_tid,
                                            ups_and_downs,
                                            print_zero_coverage_contigs,
                                            print_stream);
            }
            if let CoverageEstimatorMethods::CoverageFractionGenomeCoverageEstimator(CoverageFractionGenomeCoverageEstimator) = coverage_method {
                let mut coverage_estimator = CoverageFractionGenomeCoverageEstimator;
                contig_print_last_reference(coverage_estimator,
                                            target_names,
                                            last_tid,
                                            ups_and_downs,
                                            print_zero_coverage_contigs,
                                            print_stream);
            }
            if let CoverageEstimatorMethods::TrimmedMeanGenomeCoverageEstimator(TrimmedMeanGenomeCoverageEstimator) = coverage_method {
                let mut coverage_estimator = TrimmedMeanGenomeCoverageEstimator;
                contig_print_last_reference(coverage_estimator,
                                            target_names,
                                            last_tid,
                                            ups_and_downs,
                                            print_zero_coverage_contigs,
                                            print_stream);
            }
            if let CoverageEstimatorMethods::VarianceGenomeCoverageEstimator(VarianceGenomeCoverageEstimator) = coverage_method {
                let mut coverage_estimator = VarianceGenomeCoverageEstimator;
                contig_print_last_reference(coverage_estimator,
                                            target_names,
                                            last_tid,
                                            ups_and_downs,
                                            print_zero_coverage_contigs,
                                            print_stream);
            }
//            if last_tid != -1 {
//                coverage_estimator.add_contig(&ups_and_downs);
//                let coverage = coverage_estimator.calculate_coverage(0);
//                coverage_estimator.print_coverage(
//                    &coverage,
//                    print_stream);
//            }
//            // print zero coverage contigs at the end
//            if print_zero_coverage_contigs {
//                print_previous_zero_coverage_contigs(last_tid, target_names.len() as i32, stoit_name, coverage_estimator, &target_names, print_stream);
//            }
        }
        debug!("Outside loop");
        // print the last ref, unless there was no alignments
        print_contig(stoit_name, std::str::from_utf8(target_names[last_tid as usize]).unwrap(), print_stream);

        bam_generated.finish();
    }
}

fn contig_print_reference<T: MosdepthGenomeCoverageEstimator<T>>(
    mut coverage_estimator: T,
    mut tid: i32,
    mut last_tid: i32,
    target_names: &Vec<&[u8]>,
    header: usize,
    mut ups_and_downs: Vec<i32>,
    print_zero_coverage_contigs: bool,
    print_stream: &mut std::io::Write
    ){
    if tid != last_tid {
        if last_tid != -1 {
            coverage_estimator.add_contig(&ups_and_downs);
            let coverage = coverage_estimator.calculate_coverage(0);

            if coverage > 0.0 {
                coverage_estimator.print_coverage(
                    &coverage,
                    print_stream);
            } else if print_zero_coverage_contigs {
                coverage_estimator.print_zero_coverage(
                    print_stream);
            }
        }
        // reset for next time
        coverage_estimator.setup();
        if print_zero_coverage_contigs {
            print_previous_zero_coverage_contigs(last_tid, tid,  &coverage_estimator, print_stream);
        }
        ups_and_downs = vec![0; header];
        debug!("Working on new reference {}",
               std::str::from_utf8(target_names[tid as usize]).unwrap());
        last_tid = tid;
    }
}

fn contig_print_last_reference<T: MosdepthGenomeCoverageEstimator<T>>(
    mut coverage_estimator: T,
    target_names: Vec<&[u8]>,
    mut last_tid: i32,
    mut ups_and_downs: Vec<i32>,
    print_zero_coverage_contigs: bool,
    print_stream: &mut std::io::Write){
    if last_tid != -1 {
        coverage_estimator.add_contig(&ups_and_downs);
        let coverage = coverage_estimator.calculate_coverage(0);
        coverage_estimator.print_coverage(
            &coverage,
            print_stream);
    }
    // print zero coverage contigs at the end
    if print_zero_coverage_contigs {
        print_previous_zero_coverage_contigs(last_tid, target_names.len() as i32, &coverage_estimator, print_stream);
    }
}

fn print_previous_zero_coverage_contigs<T>(
    last_tid: i32,
    current_tid: i32,
    coverage_estimator: &MosdepthGenomeCoverageEstimator<T>,
    print_stream: &mut std::io::Write) {
    let mut my_tid = last_tid + 1;
    while my_tid < current_tid {
        coverage_estimator.print_zero_coverage(print_stream);
        my_tid += 1;
    };
}

fn print_contig<'a >(stoit_name: &str,
                     contig: &str,
                     print_stream: &'a mut std::io::Write) -> &'a mut std::io::Write {
    write!(print_stream, "\n{}\t{}",
           stoit_name,
           contig).unwrap();
    return print_stream;
}


//
//#[cfg(test)]
//mod tests {
//    use super::*;
//    use std::io::Cursor;
//    use std::str;
//
//    #[test]
//    fn test_one_genome_two_contigs_first_covered_no_zeros(){
//        let mut stream = Cursor::new(Vec::new());
//        contig_coverage(
//            generate_named_bam_readers_from_bam_files(
//                vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
//            &mut stream,
//            &mut MeanGenomeCoverageEstimator::new(0.0),
//            false,
//            false);
//        assert_eq!(
//            "7seqs.reads_for_seq1_and_seq2\tgenome2~seq1\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome5~seq2\t1.2\n",
//            str::from_utf8(stream.get_ref()).unwrap())
//    }
//
//    #[test]
//    fn test_one_genome_two_contigs_first_covered(){
//        let mut stream = Cursor::new(Vec::new());
//        contig_coverage(
//            generate_named_bam_readers_from_bam_files(
//                vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
//            &mut stream,
//            &mut MeanGenomeCoverageEstimator::new(0.0),
//            true,
//            false);
//        assert_eq!(
//            "7seqs.reads_for_seq1_and_seq2\tgenome1~random_sequence_length_11000\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome1~random_sequence_length_11010\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome2~seq1\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome3~random_sequence_length_11001\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome4~random_sequence_length_11002\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome5~seq2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6~random_sequence_length_11003\t0.0\n",
//            str::from_utf8(stream.get_ref()).unwrap())
//    }
//
//    #[test]
//    fn test_flag_filtering(){
//        let mut stream = Cursor::new(Vec::new());
//        contig_coverage(
//            generate_named_bam_readers_from_bam_files(
//                vec!["tests/data/1.bam"]),
//            &mut stream,
//            &mut MeanGenomeCoverageEstimator::new(0.0),
//            false,
//            true);
//        assert_eq!(
//            "",
//            str::from_utf8(stream.get_ref()).unwrap())
//    }
//
//    #[test]
//    fn test_one_contig_variance(){
//        let mut stream = Cursor::new(Vec::new());
//        contig_coverage(
//            generate_named_bam_readers_from_bam_files(
//                vec!["tests/data/2seqs.reads_for_seq1.bam"]),
//            &mut stream,
//            &mut VarianceGenomeCoverageEstimator::new(0.0),
//            true,
//            false);
//        assert_eq!(
//            "2seqs.reads_for_seq1\tseq1\t0.9489489\n".to_owned()+
//                "2seqs.reads_for_seq1\tseq2\t0.0\n",
//            str::from_utf8(stream.get_ref()).unwrap())
//    }
//
//    #[test]
//    fn test_streaming_bam_file(){
//        let mut stream = Cursor::new(Vec::new());
//        contig_coverage(
//            vec![
//                generate_named_bam_readers_from_read_couple(
//                    "tests/data/7seqs.fna",
//                    "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
//                    "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
//                    4)],
//            &mut stream,
//            &mut MeanGenomeCoverageEstimator::new(0.0),
//            false,
//            false);
//        assert_eq!(
//            "7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome2~seq1\t1.2\n7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome5~seq2\t1.2\n",
//            str::from_utf8(stream.get_ref()).unwrap())
//    }
//
//}
