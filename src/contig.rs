use std;
use rust_htslib::bam;

use mosdepth_genome_coverage_estimators::*;
use bam_generator::*;


pub fn contig_coverage<T: MosdepthGenomeCoverageEstimator<T>,
                       R: NamedBamReader,
                       G: NamedBamReaderGenerator<R>>(
    bam_readers: Vec<G>,
    coverage_estimator: &mut T,
    print_zero_coverage_contigs: bool,
    flag_filtering: bool) -> Vec<OutputStream>{

    let mut output_vec: Vec<OutputStream> = Vec::new();
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
            if tid != last_tid {
                if last_tid != -1 {
                    coverage_estimator.add_contig(&ups_and_downs);
                    let coverage = coverage_estimator.calculate_coverage(0);
                    if coverage > 0.0 {
                        let mut output = OutputStream::new(
                            stoit_name.to_string(),
                            std::str::from_utf8(target_names[last_tid as usize]).unwrap().to_string(),
                            coverage
                        );
                        if coverage_estimator.is_histogram(){
                            let mut hist_vec = coverage_estimator.add_to_stream(output);
                            for hvec in hist_vec{
                                output_vec.push(hvec);
                            }
                        }else{
                            output_vec.push(output);
                        }
                    } else if print_zero_coverage_contigs {
                        let mut output = OutputStream::new(
                            stoit_name.to_string(),
                            std::str::from_utf8(target_names[last_tid as usize]).unwrap().to_string(),
                            0.0
                        );
                        if coverage_estimator.is_histogram(){
                            let mut hist_vec = coverage_estimator.add_to_stream(output);
                            for hvec in hist_vec{
                                output_vec.push(hvec);
                            }
                        }else{
                            output_vec.push(output);
                        }
                        // coverage_estimator.print_genome(output);
                    }
                }
                // reset for next time
                coverage_estimator.setup();
                if print_zero_coverage_contigs {
                    let zero_vec = print_previous_zero_coverage_contigs(last_tid, tid, stoit_name, coverage_estimator, &target_names);
                    for out in zero_vec{
                        output_vec.push(out);
                    }                
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
            debug!("At end of loop")
        }
        debug!("Outside loop");
        // print the last ref, unless there was no alignments
        if last_tid != -1 {
            coverage_estimator.add_contig(&ups_and_downs);
            let coverage = coverage_estimator.calculate_coverage(0);
            let mut output = OutputStream::new(
                stoit_name.to_string(),
                std::str::from_utf8(target_names[last_tid as usize]).unwrap().to_string(),
                coverage
            );
            if coverage_estimator.is_histogram(){
                let mut hist_vec = coverage_estimator.add_to_stream(output);
                for hvec in hist_vec{
                    output_vec.push(hvec);
                }
            }else{
                output_vec.push(output);
            }
        }
        // print zero coverage contigs at the end
        if print_zero_coverage_contigs {
            let zero_vec = print_previous_zero_coverage_contigs(last_tid, target_names.len() as i32, stoit_name, coverage_estimator, &target_names);
            for out in zero_vec{
                output_vec.push(out);
            }
        }

    bam_generated.finish();
    }
    return output_vec;
}


fn print_previous_zero_coverage_contigs<'a, T>(
    last_tid: i32,
    current_tid: i32,
    stoit_name: &str,
    coverage_estimator: &'a mut MosdepthGenomeCoverageEstimator<T>,
    target_names: &Vec<&[u8]>)
    -> Vec<OutputStream> {
    
    let mut output_vec = Vec::new();
    let mut my_tid = last_tid + 1;
    while my_tid < current_tid {
        let mut output = OutputStream::new(
            stoit_name.to_string(), 
            std::str::from_utf8(target_names[my_tid as usize]).unwrap().to_string(),
            0.0);
        output_vec.push(output);
        
        my_tid += 1;
    }    
    return output_vec;
}

#[cfg(test)]
mod tests {

    use super::*;

    pub fn print_output_stream(output_stream: Vec<OutputStream>) -> String {
        let mut out_st = Vec::new();
        for mut out in output_stream{
            out_st.push(print_output(out));
        }
        return out_st.join("")
    }
    pub fn print_output(out: OutputStream) -> String {
        let mut out_string = format!("{}\t{}", out.filename, out.genome);
        for c in out.methods.iter(){
            out_string=format!("{}{}", out_string,format!("\t{}", c));
        }
        out_string=format!("{}{}", out_string,"\n");
        return out_string
    }
    #[test]
    fn test_one_genome_two_contigs_first_covered_no_zeros(){
        let stream = contig_coverage(
                                    generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
                                    &mut MeanGenomeCoverageEstimator::new(0.0),
                                    false,
                                    false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome2~seq1\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome5~seq2\t1.2\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_one_genome_two_contigs_first_covered(){
       let stream = contig_coverage(generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
                                                                    &mut MeanGenomeCoverageEstimator::new(0.0),
                                                                    true,
                                                                    false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1~random_sequence_length_11000\t0\n7seqs.reads_for_seq1_and_seq2\tgenome1~random_sequence_length_11010\t0\n7seqs.reads_for_seq1_and_seq2\tgenome2~seq1\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome3~random_sequence_length_11001\t0\n7seqs.reads_for_seq1_and_seq2\tgenome4~random_sequence_length_11002\t0\n7seqs.reads_for_seq1_and_seq2\tgenome5~seq2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6~random_sequence_length_11003\t0\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_flag_filtering(){
        let stream = contig_coverage(
                            generate_named_bam_readers_from_bam_files(vec!["tests/data/1.bam"]),
                            &mut MeanGenomeCoverageEstimator::new(0.0),
                            false,
                            true);
        assert_eq!(
            "",
            print_output_stream(stream))
    }

    #[test]
    fn test_one_contig_variance(){
        let stream = contig_coverage(
                            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
                            &mut VarianceGenomeCoverageEstimator::new(0.0),
                            true,
                            false);
        assert_eq!(
            "2seqs.reads_for_seq1\tseq1\t0.9489489\n".to_owned()+
                "2seqs.reads_for_seq1\tseq2\t0\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_streaming_bam_file(){
        let stream = contig_coverage(
                            vec![
                                generate_named_bam_readers_from_read_couple(
                                    "tests/data/7seqs.fna",
                                    "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                                    "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                                    4)],
                            &mut MeanGenomeCoverageEstimator::new(0.0),
                            false,
                            false);
        assert_eq!(
            "7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome2~seq1\t1.2\n7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome5~seq2\t1.2\n",
            print_output_stream(stream))
    }

}
