use std;

use rust_htslib::bam;
use rust_htslib::bam::Read as bamRead;

use std::str;
use std::collections::BTreeSet;

use mosdepth_genome_coverage_estimators::*;
use genomes_and_contigs::GenomesAndContigs;

pub fn mosdepth_genome_coverage_with_contig_names<T: 'static +  'static +  MosdepthGenomeCoverageEstimator<T> + std::fmt::Debug>(
    bam_files: &Vec<&str>,
    contigs_and_genomes: &GenomesAndContigs,
    limit_stream: bool,
    coverage_estimators: &mut T,
    print_zero_coverage_genomes: bool,
    flag_filtering: bool) -> Vec<OutputStream>{

    let mut output_vec: Vec<OutputStream> = Vec::new();
    for bam_file in bam_files {
        debug!("Working on BAM file {}", bam_file);
        let mut bam = bam::Reader::from_path(bam_file).expect(
            &format!("Unable to find BAM file {}", bam_file));
        let header = bam.header().clone();
        let target_names = header.target_names();

        // Collect reference numbers for each genome's contigs
        let mut reference_number_to_genome_index: Vec<usize> = vec![];
        for name in target_names {
            let genome_index = contigs_and_genomes.genome_index_of_contig(
                &String::from(std::str::from_utf8(name)
                              .expect("UTF8 encoding error in BAM header file")))
                .expect(
                    &format!("Found reference that was not part of any genome: {:?}", name));
            reference_number_to_genome_index.push(genome_index);
        }
        debug!("Read in {} reference IDs mapped to genomes {:?}",
               reference_number_to_genome_index.len(), reference_number_to_genome_index);


            let mut per_genome_coverage_estimators: Vec<T> = vec!();
            for _ in contigs_and_genomes.genomes.iter() {
                let cov_clone = coverage_estimators.copy();
                per_genome_coverage_estimators.push(cov_clone);
            }


            // Iterate through bam records
            let mut last_tid: u32 = 0;
            let mut doing_first = true;
            let mut ups_and_downs: Vec<i32> = Vec::new();
            let stoit_name = std::path::Path::new(bam_file).file_stem().unwrap().to_str().expect(
                "failure to convert bam file name to stoit name - UTF8 error maybe?");
            let mut record: bam::record::Record = bam::record::Record::new();
            let mut seen_ref_ids = BTreeSet::new();
            while bam.read(&mut record).is_ok() {
                if flag_filtering &&
                    (record.is_secondary() ||
                     record.is_supplementary() ||
                     !record.is_proper_pair()) {
                        continue;
                    }
                let tid = record.tid() as u32;
                if tid != last_tid || doing_first {
                    debug!("Came across a new tid {}", tid);
                    if doing_first == true {
                        doing_first = false;
                    } else {
                        per_genome_coverage_estimators[reference_number_to_genome_index[last_tid as usize]]
                            .add_contig(&ups_and_downs);
                    }

                    ups_and_downs = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                    last_tid = tid;
                    seen_ref_ids.insert(tid);
                }

                // TODO: move below into a function for code-reuse purposes.
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

            // Record the last contig
            per_genome_coverage_estimators[reference_number_to_genome_index[last_tid as usize]]
                .add_contig(&ups_and_downs);

            // Print the coverages of each genome
            // Calculate the unobserved length of each genome
            let mut unobserved_lengths: Vec<u32> = vec!();
            for _ in 0..contigs_and_genomes.genomes.len() {
                unobserved_lengths.push(0)
            }
            debug!("estimators: {:?}", per_genome_coverage_estimators);
            for (ref_id, genome_id) in reference_number_to_genome_index.iter().enumerate() {
                let ref_id_u32: u32 = ref_id as u32;
                debug!("Seen {:?}", seen_ref_ids);
                if !seen_ref_ids.contains(&ref_id_u32) {
                    debug!("Getting target #{} from header names", ref_id_u32);
                    unobserved_lengths[*genome_id] += header.target_len(ref_id_u32).unwrap()
                }
            }
            // print the genomes out
            for (i, genome) in contigs_and_genomes.genomes.iter().enumerate() {
                let coverage = per_genome_coverage_estimators[i]
                    .calculate_coverage(unobserved_lengths[i]);

                // Print coverage of previous genome
                debug!("Found coverage {} for genome {}", coverage, genome);

                // let mut output = per_genome_coverage_estimators[i].create_output();
                if coverage > 0.0 {
                    let mut output = OutputStream::new(
                        stoit_name.to_string(),
                        genome.to_string(),
                        coverage
                    );
                    if limit_stream{
                        per_genome_coverage_estimators[i].print_genome(output.clone());
                    }
                    output_vec.push(output);
                } else if print_zero_coverage_genomes {
                    let mut output = OutputStream::new(
                        stoit_name.to_string(),
                        genome.to_string(),
                        0.0
                    );
                    if limit_stream{
                        per_genome_coverage_estimators[i].print_genome(output.clone());
                    }
                    output_vec.push(output);
                }

            }
        }
        return output_vec;
    }




pub fn mosdepth_genome_coverage<T: MosdepthGenomeCoverageEstimator<T>>(
    bam_files: &Vec<&str>,
    split_char: u8,
    limit_stream: bool,
    coverage_estimator: &mut T,
    print_zero_coverage_genomes: bool,
    flag_filtering: bool,
    single_genome: bool) -> Vec<OutputStream>{

    let mut output_vec: Vec<OutputStream> = Vec::new();
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
                if single_genome ||
                    extract_genome(my_tid, &target_names, split_char) == target_genome {

                    extra += header.target_len(my_tid)
                        .expect("Malformed bam header or programming error encountered");
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
                if single_genome ||
                    extract_genome(my_tid, &target_names, split_char) == target_genome {

                    extra += header.target_len(my_tid)
                        .expect("Malformed bam header or programming error encountered");
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
            let mut my_tid = last_tid + 1;
            while my_tid < current_tid {
                if single_genome ||
                    extract_genome(my_tid, &target_names, split_char) == target_genome {

                    extra += header.target_len(my_tid)
                        .expect("Malformed bam header or programming error encountered");
                    my_tid += 1;
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
            if flag_filtering &&
                (record.is_secondary() ||
                 record.is_supplementary() ||
                 !record.is_proper_pair()) {
                    debug!("skipping record {:?} as it filters out based on flags", record);
                    continue;
                }
            // if reference has changed, finish a genome or not
            let tid = record.tid() as u32;
            let current_genome: &[u8] = match single_genome {
                true => "".as_bytes(),
                false => extract_genome(tid as u32, &target_names, split_char)
            };

            if tid != last_tid || doing_first {
                if doing_first == true {
                    coverage_estimator.setup();
                    last_genome = current_genome;
                    unobserved_contig_length = fill_genome_length_backwards(tid, current_genome);
                    doing_first = false;
                    if print_zero_coverage_genomes && !single_genome {
                        print_previous_zero_coverage_genomes2(
                            stoit_name, b"", current_genome, tid, coverage_estimator,
                            &target_names, split_char);
                    }

                } else if current_genome == last_genome {
                    coverage_estimator.add_contig(&ups_and_downs);
                    // Collect the length of reference sequences from this
                    // genome that had no hits that were just skipped over.
                    debug!("Filling unobserved from {} to {}", last_tid, tid);
                    unobserved_contig_length += fill_genome_length_backwards_to_last(
                        tid, last_tid as u32, current_genome);

                } else {
                    coverage_estimator.add_contig(&ups_and_downs);
                    // Collect the length of refs from the end of the last genome that had no hits
                    debug!("Filling unobserved from {} to {} for {}", last_tid, tid, &str::from_utf8(last_genome).unwrap());
                    unobserved_contig_length += fill_genome_length_backwards_to_last(
                        tid, last_tid as u32, last_genome);
                    debug!("unobserved_contig_length now {}", unobserved_contig_length);
                    // Determine coverage of previous genome
                    let coverage = coverage_estimator.calculate_coverage(unobserved_contig_length);

                    // Print coverage of previous genome
                    if coverage > 0.0 {
                        let mut output = OutputStream::new(
                            stoit_name.to_string(),
                            str::from_utf8(last_genome).unwrap().to_string(),
                            coverage
                        );
                        if limit_stream{
                            coverage_estimator.print_genome(output.clone());
                        }
                        output_vec.push(output);

                    } else if print_zero_coverage_genomes {
                        let mut output = OutputStream::new(
                            stoit_name.to_string(),
                            str::from_utf8(last_genome).unwrap().to_string(),
                            0.0
                        );
                        if limit_stream{
                            coverage_estimator.print_genome(output.clone());
                        }
                        output_vec.push(output);
                    }
                    coverage_estimator.setup();
                    if print_zero_coverage_genomes {
                        print_previous_zero_coverage_genomes2(
                            stoit_name, last_genome, current_genome, tid, coverage_estimator,
                            &target_names, split_char);
                    }
                    last_genome = current_genome;
                    unobserved_contig_length = fill_genome_length_backwards(tid, current_genome);
                    debug!("Setting unobserved contig length to be {}", unobserved_contig_length);
                }

                ups_and_downs = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                last_tid = tid;
            }


            // Add coverage info for the current record
            // for each chunk of the cigar string
            debug!("read name {:?}", std::str::from_utf8(record.qname()).unwrap());
            let mut cursor: usize = record.pos() as usize;
            for cig in record.cigar().iter() {
                //debug!("Found cigar {:} from {}", cig, cursor);
                match cig.char() {
                    'M' => {
                        // if M, increment start and decrement end index
                        //debug!("Adding M at {} and {}", cursor, cursor + cig.len() as usize);
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
        debug!("At end, unobserved_contig_length now {}", unobserved_contig_length);
        // Determine coverage of previous genome
        let coverage = coverage_estimator.calculate_coverage(unobserved_contig_length);

        // Give the single genome a dummy name
        if single_genome {
            last_genome = "genome1".as_bytes()
        }

        // Print coverage of previous genome

        if coverage > 0.0 {
            let mut output = OutputStream::new(
                stoit_name.to_string(),
                str::from_utf8(last_genome).unwrap().to_string(),
                coverage
            );
            if limit_stream{
                coverage_estimator.print_genome(output.clone());
            }
            output_vec.push(output);
        } else if print_zero_coverage_genomes {
            let mut output = OutputStream::new(
                stoit_name.to_string(),
                str::from_utf8(last_genome).unwrap().to_string(),
                0.0
            );
            if limit_stream{
                coverage_estimator.print_genome(output.clone());
            }
            output_vec.push(output);
        }
        if print_zero_coverage_genomes && !single_genome {
            print_previous_zero_coverage_genomes2(
                stoit_name, last_genome, b"", header.target_count()-1, coverage_estimator,
                &target_names, split_char);
        }
    }
    return output_vec;
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
fn print_previous_zero_coverage_genomes2<'a, T>(
    stoit_name: &str,
    last_genome: &[u8],
    current_genome: &[u8],
    current_tid: u32,
    pileup_coverage_estimator: &'a MosdepthGenomeCoverageEstimator<T>,
    target_names: &Vec<&[u8]>,
    split_char: u8)
    -> &'a MosdepthGenomeCoverageEstimator<T> {

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
        // pileup_coverage_estimator.print_zero_coverage(
        //     &stoit_name, &str::from_utf8(genomes_to_print[i]).unwrap(), print_stream);
        let mut output = OutputStream::new(
            stoit_name.to_string(),
            str::from_utf8(genomes_to_print[i]).unwrap().to_string(),
            0.0
        );
        pileup_coverage_estimator.print_genome(output.clone());

            // output_vec.push(output);
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
            true,
            false,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1\tse\t0.6\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_two_contigs_first_covered_contig_names(){
        let mut stream = Cursor::new(Vec::new());
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("se".to_string());
        geco.insert("seq1".to_string(),genome1);
        geco.insert("seq2".to_string(),genome1);
        mosdepth_genome_coverage_with_contig_names(
            &vec!["test/data/2seqs.reads_for_seq1.bam"],
            &geco,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true,
            false);
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
            true,
            false,
            false);
        assert_eq!(
            "2seqs.reads_for_seq2\tse\t0.6\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_two_contigs_second_covered_contig_names(){
        let mut stream = Cursor::new(Vec::new());
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("se".to_string());
        geco.insert("seq1".to_string(),genome1);
        geco.insert("seq2".to_string(),genome1);
        mosdepth_genome_coverage_with_contig_names(
            &vec!["test/data/2seqs.reads_for_seq2.bam"],
            &geco,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true,
            false);
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
            true,
            false,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_two_contigs_both_covered_contig_names(){
        let mut stream = Cursor::new(Vec::new());
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(),genome1);
        geco.insert("seq2".to_string(),genome1);
        mosdepth_genome_coverage_with_contig_names(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            &geco,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true,
            false);
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
            true,
            false,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t0.0\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_min_fraction_covered_under_min_contig_names(){
        let mut stream = Cursor::new(Vec::new());
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(),genome1);
        geco.insert("seq2".to_string(),genome1);
        mosdepth_genome_coverage_with_contig_names(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            &geco,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.76),
            false,
            false);
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
            true,
            false,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }


    #[test]
    fn test_one_genome_min_fraction_covered_just_ok_contig_names(){
        let mut stream = Cursor::new(Vec::new());
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(),genome1);
        geco.insert("seq2".to_string(),genome1);
        mosdepth_genome_coverage_with_contig_names(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            &geco,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.759),
            true,
            false);
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
            true,
            false,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.08875\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_two_contigs_trimmed_mean_contig_names(){
        let mut stream = Cursor::new(Vec::new());
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(),genome1);
        geco.insert("seq2".to_string(),genome1);
        mosdepth_genome_coverage_with_contig_names(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            &geco,
            &mut stream,
            &mut TrimmedMeanGenomeCoverageEstimator::new(0.1, 0.9, 0.0),
            true,
            false);
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
            true,
            false,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t0\t482\n2seqs.reads_for_seq1_and_seq2\ts\t1\t922\n2seqs.reads_for_seq1_and_seq2\ts\t2\t371\n2seqs.reads_for_seq1_and_seq2\ts\t3\t164\n2seqs.reads_for_seq1_and_seq2\ts\t4\t61\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_two_contigs_pileup_counts_estimator_contig_names(){
        let mut stream = Cursor::new(Vec::new());
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(),genome1);
        geco.insert("seq2".to_string(),genome1);
        mosdepth_genome_coverage_with_contig_names(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            &geco,
            &mut stream,
            &mut PileupCountsGenomeCoverageEstimator::new(0.0),
            true,
            false);
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
            true,
            false,
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome3\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome4\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6\t0.0\n",
            str::from_utf8(stream.get_ref()).unwrap());

        stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            &vec!["test/data/7seqs.reads_for_seq1_and_seq2.bam"],
            '~' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.1),
            false,
            false,
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_zero_coverage_genomes_after_min_fraction(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            &vec!["test/data/7seqs.reads_for_seq1_and_seq2.bam"],
            '~' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.76),
            true,
            false,
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome2\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome3\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome4\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6\t0.0\n",
            str::from_utf8(stream.get_ref()).unwrap());
    }

    #[test]
    fn test_single_genome(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            &vec!["test/data/7seqs.reads_for_seq1_and_seq2.bam"],
            '~' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true,
            false,
            true);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0.04209345\n",
            str::from_utf8(stream.get_ref()).unwrap());
    }

    #[test]
    fn test_zero_coverage_genomes_contig_names(){
        let mut stream = Cursor::new(Vec::new());
        let mut geco = GenomesAndContigs::new();
        // >genome1~random_sequence_length_11000
        //     >genome1~random_sequence_length_11010
        //     >genome2~seq1
        //     >genome3~random_sequence_length_11001
        //     >genome4~random_sequence_length_11002
        //     >genome5~seq2
        //     >genome6~random_sequence_length_11003
        let genome1 = geco.establish_genome("genome1".to_string());
        let genome2 = geco.establish_genome("genome2".to_string());
        let genome3 = geco.establish_genome("genome3".to_string());
        let genome4 = geco.establish_genome("genome4".to_string());
        let genome5 = geco.establish_genome("genome5".to_string());
        let genome6 = geco.establish_genome("genome6".to_string());
        geco.insert("genome1~random_sequence_length_11000".to_string(),genome1);
        geco.insert("genome1~random_sequence_length_11010".to_string(),genome1);
        geco.insert("genome2~seq1".to_string(),genome2);
        geco.insert("genome3~random_sequence_length_11001".to_string(),genome3);
        geco.insert("genome4~random_sequence_length_11002".to_string(),genome4);
        geco.insert("genome5~seq2".to_string(),genome5);
        geco.insert("genome6~random_sequence_length_11003".to_string(),genome6);
        mosdepth_genome_coverage_with_contig_names(
            &vec!["test/data/7seqs.reads_for_seq1_and_seq2.bam"],
            &geco,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.1),
            true,
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome3\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome4\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6\t0.0\n",
            str::from_utf8(stream.get_ref()).unwrap());

        stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage_with_contig_names(
            &vec!["test/data/7seqs.reads_for_seq1_and_seq2.bam"],
            &geco,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.1),
            false,
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }
}
