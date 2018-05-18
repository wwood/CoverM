use std;

use rust_htslib::bam;
use rust_htslib::bam::Read as bamRead;

use std::str;

use mosdepth_genome_coverage_estimators::*;






pub fn mosdepth_genome_coverage<T: MosdepthGenomeCoverageEstimator>(
    bam_files: &Vec<&str>,
    split_char: u8,
    print_stream: &mut std::io::Write,
    coverage_estimator: &mut T,
    print_zero_coverage_genomes: bool) {

    for bam_file in bam_files {
        debug!("Working on BAM file {}", bam_file);
        let mut bam = bam::Reader::from_path(bam_file).expect(
            &format!("Unable to find BAM file {}", bam_file));
        let header = bam.header().clone();
        let target_names = header.target_names();

        let fill_genome_length_forwards = |current_tid, target_genome| {
            // Iterating reads skips over contigs with no mapped reads, but the
            // length of these contigs is required to calculate the average
            // across all contigs. This closure returns the number of bases in
            // contigs with tid > current_tid that are part of the current
            // genome.
            let mut extra: u32 = 0;
            let total_refs = header.target_count();
            let mut my_tid = current_tid + 1;
            while my_tid < total_refs {
                let my_genome = extract_genome(my_tid, &target_names, split_char);
                if my_genome == target_genome {
                    extra += header.target_len(my_tid).expect("Malformed bam header or programming error encountered");
                    my_tid += 1;
                } else {
                    break;
                }
            }
            return extra
        };
        let fill_genome_length_backwards = |current_tid, target_genome| {
            if current_tid == 0 {return 0}
            let mut extra: u32 = 0;
            let mut my_tid = current_tid - 1;
            loop {
                let my_genome = extract_genome(my_tid, &target_names, split_char);
                if my_genome == target_genome {
                    extra += header.target_len(my_tid).expect("Malformed bam header or programming error encountered");
                    if my_tid == 0 {
                        break
                    } else {
                        my_tid -= 1;
                    }
                } else {
                    break;
                }
            }
            return extra
        };
        let fill_genome_length_backwards_to_last = |current_tid, last_tid, target_genome| {
            if current_tid == 0 {return 0};
            let mut extra: u32 = 0;
            let mut my_tid = current_tid - 1;
            while my_tid > last_tid {
                let my_genome = extract_genome(my_tid, &target_names, split_char);
                if my_genome == target_genome {
                    extra += header.target_len(my_tid).expect("Malformed bam header or programming error encountered");
                    my_tid -= 1;
                } else {
                    break;
                }
            }
            return extra
        };


        let mut last_tid: u32 = 0;
        let mut doing_first = true;
        let mut last_genome: &[u8] = "error genome".as_bytes();
        let mut unobserved_contig_length: u32 = 0;
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let stoit_name = std::path::Path::new(bam_file).file_stem().unwrap().to_str().expect(
            "failure to convert bam file name to stoit name - UTF8 error maybe?");
        let mut record: bam::record::Record = bam::record::Record::new();
        while bam.read(&mut record).is_ok() {
            // if reference has changed, finish a genome or not
            let tid = record.tid() as u32;
            let current_genome = extract_genome(tid as u32, &target_names, split_char);
            if tid != last_tid || doing_first {
                if doing_first == true {
                    coverage_estimator.setup();
                    last_genome = current_genome;
                    unobserved_contig_length = fill_genome_length_backwards(tid, current_genome);
                    doing_first = false;
                    if print_zero_coverage_genomes {
                        print_previous_zero_coverage_genomes2(
                            stoit_name, b"", current_genome, tid, coverage_estimator,
                            &target_names, split_char, print_stream);
                    }

                } else if current_genome == last_genome {
                    coverage_estimator.add_contig(&ups_and_downs);
                    // Collect the length of reference sequences from this
                    // genome that had no hits that were just skipped over.
                    unobserved_contig_length += fill_genome_length_backwards_to_last(
                        tid, last_tid as u32, current_genome);

                } else {
                    coverage_estimator.add_contig(&ups_and_downs);
                    // Collect the length of refs from the end of the last genome that had no hits
                    unobserved_contig_length += fill_genome_length_backwards_to_last(
                        tid, last_tid as u32, current_genome);
                    // Determine coverage of previous genome
                    let coverage = coverage_estimator.calculate_coverage(unobserved_contig_length);

                    // Print coverage of previous genome
                    if coverage > 0.0 {
                        coverage_estimator.print_genome(
                            &stoit_name,
                            &str::from_utf8(last_genome).unwrap(),
                            &coverage,
                            print_stream);
                    }
                    coverage_estimator.setup();
                    if print_zero_coverage_genomes {
                        print_previous_zero_coverage_genomes2(
                            stoit_name, last_genome, current_genome, tid, coverage_estimator,
                            &target_names, split_char, print_stream);
                    }
                    last_genome = current_genome;
                    unobserved_contig_length = fill_genome_length_backwards(tid, current_genome);
                }

                ups_and_downs = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                last_tid = tid;
            }


            // Add coverage info for the current record
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
        // Print the last genome
        coverage_estimator.add_contig(&ups_and_downs);
        // Collect the length of refs from the end of the last genome that had no hits
        unobserved_contig_length += fill_genome_length_forwards(last_tid, last_genome);
        // Determine coverage of previous genome
        let coverage = coverage_estimator.calculate_coverage(unobserved_contig_length);

        // Print coverage of previous genome
        if coverage > 0.0 {
            coverage_estimator.print_genome(
                &stoit_name,
                &str::from_utf8(last_genome).unwrap(),
                &coverage,
                print_stream);
        }
        if print_zero_coverage_genomes {
            print_previous_zero_coverage_genomes2(
                stoit_name, last_genome, b"", header.target_count()-1, coverage_estimator,
                &target_names, split_char, print_stream);
        }
    }
}




fn extract_genome<'a>(tid: u32, target_names: &'a Vec<&[u8]>, split_char: u8) -> &'a [u8] {
    let target_name = target_names[tid as usize];
    debug!("target name {:?}, separator {:?}", target_name, split_char);
    let offset = find_first(target_name, split_char).expect(
        &format!("Contig name {} does not contain split symbol, so cannot determine which genome it belongs to",
                 str::from_utf8(target_name).unwrap()));
    return &target_name[(0..offset)];
}


// Print zero coverage for genomes that have no reads mapped. Genomes are
// detected from the header, counting backwards from the current tid until the
// last seen genome is encountered, or we reach the beginning of the tid array.
fn print_previous_zero_coverage_genomes2<'a>(
    stoit_name: &str,
    last_genome: &[u8],
    current_genome: &[u8],
    current_tid: u32,
    pileup_coverage_estimator: &'a MosdepthGenomeCoverageEstimator,
    target_names: &Vec<&[u8]>,
    split_char: u8,
    print_stream: &mut std::io::Write)
    -> &'a MosdepthGenomeCoverageEstimator {

    let mut my_current_genome = current_genome;
    let mut tid = current_tid;
    let mut genomes_to_print: Vec<&[u8]> = vec![];
    while tid > 0 {
        let genome = extract_genome(tid, &target_names, split_char);
        if genome == last_genome { break; }
        else if genome != my_current_genome {
            // In-between genome encountered for the first time.
            genomes_to_print.push(genome);
            my_current_genome = genome;
        }
        tid = tid - 1;
    };
    for i in (0..genomes_to_print.len()).rev() {
        pileup_coverage_estimator.print_zero_coverage(
            &stoit_name, &str::from_utf8(genomes_to_print[i]).unwrap(), print_stream);
    }
    return pileup_coverage_estimator;
}

/// Finds the first occurence of element in a slice
fn find_first<T>(slice: &[T], element: T) -> Result<usize, &'static str>
    where T: std::cmp::PartialEq<T> {

    let mut index: usize = 0;
    for el in slice {
        if *el == element {
            return Ok(index)
            //let res: Result<usize, None> = Ok(index)
        }
        index += 1;
    }
    return Err("Element not found in slice")
}





#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use env_logger;

    #[test]
    fn initialize_logger() {
        env_logger::init().unwrap();
    }

    #[test]
    fn test_one_genome_two_contigs_first_covered(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            &vec!["test/data/2seqs.reads_for_seq1.bam"],
            'q' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true);
        assert_eq!(
            "2seqs.reads_for_seq1\tse\t0.6\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_two_contigs_second_covered(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            &vec!["test/data/2seqs.reads_for_seq2.bam"],
            'q' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true);
        assert_eq!(
            "2seqs.reads_for_seq2\tse\t0.6\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_two_contigs_both_covered(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            'e' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_min_fraction_covered_under_min(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            'e' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.76),
            true);
        assert_eq!(
            "",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_min_fraction_covered_just_ok(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            'e' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.759),
            true);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_two_contigs_trimmed_mean(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            'e' as u8,
            &mut stream,
            &mut TrimmedMeanGenomeCoverageEstimator::new(0.1, 0.9, 0.0),
            true);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.08875\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_two_contigs_pileup_counts_estimator(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            'e' as u8,
            &mut stream,
            &mut PileupCountsGenomeCoverageEstimator::new(0.0),
            true);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t0\t482\n2seqs.reads_for_seq1_and_seq2\ts\t1\t922\n2seqs.reads_for_seq1_and_seq2\ts\t2\t371\n2seqs.reads_for_seq1_and_seq2\ts\t3\t164\n2seqs.reads_for_seq1_and_seq2\ts\t4\t61\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_zero_coverage_genomes(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            &vec!["test/data/7seqs.reads_for_seq1_and_seq2.bam"],
            '~' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.1),
            true);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome3\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome4\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6\t0.0\n",
            str::from_utf8(stream.get_ref()).unwrap());

        stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            &vec!["test/data/7seqs.reads_for_seq1_and_seq2.bam"],
            '~' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.1),
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }
}
