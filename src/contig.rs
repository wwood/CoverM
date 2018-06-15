use std;
use rust_htslib::bam;
use rust_htslib::bam::Read;

use mosdepth_genome_coverage_estimators::*;


pub fn contig_coverage<T: MosdepthGenomeCoverageEstimator<T>>(
    bam_files: &Vec<&str>,
    print_stream: &mut std::io::Write,
    coverage_estimator: &mut T,
    print_zero_coverage_contigs: bool,
    flag_filtering: bool) {

    for bam_file in bam_files {
        debug!("Working on BAM file {}", bam_file);
        let mut bam = bam::Reader::from_path(bam_file).expect(
            &format!("Unable to find BAM file {}", bam_file));

        let mut record: bam::record::Record = bam::record::Record::new();
        let mut last_tid: i32 = -1; // no such tid in a real BAM file
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let header = bam.header().clone();
        let stoit_name = std::path::Path::new(bam_file).file_stem().unwrap().to_str().expect(
            "failure to convert bam file name to stoit name - UTF8 error maybe?");
        let target_names = header.target_names();

        // for record in records
        while bam.read(&mut record).is_ok() {
            if flag_filtering &&
                (record.is_secondary() ||
                 record.is_supplementary() ||
                 !record.is_proper_pair()) {
                    continue;
                }
            // if reference has changed, print the last record
            let tid = record.tid();
            if tid != last_tid {
                if last_tid != -1 {
                    coverage_estimator.add_contig(&ups_and_downs);
                    let coverage = coverage_estimator.calculate_coverage(0);

                    if coverage > 0.0 {
                        coverage_estimator.print_genome(
                            stoit_name,
                            std::str::from_utf8(target_names[last_tid as usize]).unwrap(),
                            &coverage,
                            print_stream);
                    } else if print_zero_coverage_contigs {
                        coverage_estimator.print_zero_coverage(
                            stoit_name,
                            std::str::from_utf8(target_names[last_tid as usize]).unwrap(),
                            print_stream);
                    }
                }
                // reset for next time
                coverage_estimator.setup();
                if print_zero_coverage_contigs {
                    print_previous_zero_coverage_contigs(last_tid, tid, stoit_name, coverage_estimator, &target_names, print_stream);
                }
                ups_and_downs = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                debug!("Working on new reference {}",
                       std::str::from_utf8(target_names[tid as usize]).unwrap());
                last_tid = tid;
            }

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
                        if final_pos < ups_and_downs.len(){ // True unless the read hits the contig end.
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
        }
        // print the last ref, unless there was no alignments
        if last_tid != -1 {
            coverage_estimator.add_contig(&ups_and_downs);
            let coverage = coverage_estimator.calculate_coverage(0);
            coverage_estimator.print_genome(
                stoit_name,
                std::str::from_utf8(target_names[last_tid as usize]).unwrap(),
                &coverage,
                print_stream);
        }
        // print zero coverage contigs at the end
        if print_zero_coverage_contigs {
            print_previous_zero_coverage_contigs(last_tid, target_names.len() as i32, stoit_name, coverage_estimator, &target_names, print_stream);
        }
    }
}


fn print_previous_zero_coverage_contigs<T>(
    last_tid: i32,
    current_tid: i32,
    stoit_name: &str,
    coverage_estimator: &MosdepthGenomeCoverageEstimator<T>,
    target_names: &Vec<&[u8]>,
    print_stream: &mut std::io::Write) {
    let mut my_tid = last_tid + 1;
    while my_tid < current_tid {
        coverage_estimator.print_zero_coverage(
            stoit_name, std::str::from_utf8(target_names[my_tid as usize]).unwrap(), print_stream);
        my_tid += 1;
    };
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use std::str;

    #[test]
    fn test_one_genome_two_contigs_first_covered_no_zeros(){
        let mut stream = Cursor::new(Vec::new());
        contig_coverage(
            &vec!["test/data/7seqs.reads_for_seq1_and_seq2.bam"],
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            false,
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome2~seq1\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome5~seq2\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_two_contigs_first_covered(){
        let mut stream = Cursor::new(Vec::new());
        contig_coverage(
            &vec!["test/data/7seqs.reads_for_seq1_and_seq2.bam"],
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true,
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1~random_sequence_length_11000\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome1~random_sequence_length_11010\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome2~seq1\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome3~random_sequence_length_11001\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome4~random_sequence_length_11002\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome5~seq2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6~random_sequence_length_11003\t0.0\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_flag_filtering(){
        let mut stream = Cursor::new(Vec::new());
        contig_coverage(
            &vec!["test/data/1.bam"],
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            false,
            true);
        assert_eq!(
            "",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_contig_variance(){
        let mut stream = Cursor::new(Vec::new());
        contig_coverage(
            &vec!["test/data/2seqs.reads_for_seq1.bam"],
            &mut stream,
            &mut VarianceGenomeCoverageEstimator::new(0.0),
            true,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1\tseq1\t0.9489489\n".to_owned()+
                "2seqs.reads_for_seq1\tseq2\t0.0\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }
}
