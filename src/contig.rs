use std;
use rust_htslib::bam;
use rust_htslib::bam::Read;


pub fn contig_coverage(
    bam_files: &Vec<&str>,
    print_stream: &mut std::io::Write){

    let process_last_reference =
        |header: &bam::HeaderView,
         target_names: &Vec<&[u8]>,
         tid,
         ups_and_downs: &Vec<i32>,
         print_stream: & mut std::io::Write,
         stoit_name: &str| {
        let mut total: u32 = 0;
        let mut cumulative_sum: i32 = 0;
        let len = header.target_len(tid).expect("bam header parsing issue");
        //debug!("ups and downs for {:}: {:?}", std::str::from_utf8(header.target_names()[tid as usize]).unwrap(), ups_and_downs);
        for i in 0..len {
            let current = ups_and_downs[i as usize];
            cumulative_sum += current;
            assert!(cumulative_sum >= 0,
                    format!("At 0-based position {} on {}, the cumulative sum became < 0, it was {}",
                            i,
                            std::str::from_utf8(target_names[tid as usize]).unwrap(),
                            cumulative_sum));
            total += cumulative_sum as u32;
        }
        debug!("at end, sum was {}, total was {}", cumulative_sum, total);
        let mean: f32 = total as f32 / len as f32;
        writeln!(print_stream, "{}\t{}\t{}",
                 stoit_name,
                 std::str::from_utf8(target_names[tid as usize]).unwrap(),
                 mean).unwrap();
    };

    for bam_file in bam_files {
        debug!("Working on BAM file {}", bam_file);
        let mut bam = bam::Reader::from_path(bam_file).expect(
            &format!("Unable to find BAM file {}", bam_file));

        // for record in records
        let mut record: bam::record::Record = bam::record::Record::new();
        let mut last_tid: i32 = -1; // no such tid in a real BAM file
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let header = bam.header().clone();
        let stoit_name = std::path::Path::new(bam_file).file_stem().unwrap().to_str().expect(
            "failure to convert bam file name to stoit name - UTF8 error maybe?");
        let target_names = header.target_names();
        while bam.read(&mut record).is_ok() {
            // if reference has changed, print the last record
            let tid = record.tid();
            if tid != last_tid {
                if last_tid != -1 {
                    process_last_reference(&header, &target_names, last_tid as u32, &ups_and_downs, print_stream, &stoit_name);
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
                        // if D, decrement start and increment end
                        // debug!("Adding D at {} and {}", cursor+1, cursor + cig.len() as usize);
                        // ups_and_downs[cursor + 1] -= 1;
                        // ups_and_downs[cursor + cig.len() as usize] += 1;
                        cursor += cig.len() as usize;
                    },
                    '=' => panic!("CIGAR '=' detected, but this case is not correctly handled for now"),
                    _ => {}
                }
            }
        }
        // print the last ref, unless there was no alignments
        if last_tid != -1 {
            process_last_reference(&header, &target_names, last_tid as u32, &ups_and_downs, print_stream, &stoit_name);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use std::str;

    #[test]
    fn test_one_genome_two_contigs_first_covered(){
        let mut stream = Cursor::new(Vec::new());
        contig_coverage(
            &vec!["test/data/7seqs.reads_for_seq1_and_seq2.bam"],
            &mut stream);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome2~seq1\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome5~seq2\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }
}
