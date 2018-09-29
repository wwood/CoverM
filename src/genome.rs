use std;

use rust_htslib::bam;

use std::str;
use std::collections::BTreeSet;

use mosdepth_genome_coverage_estimators::*;
use genomes_and_contigs::GenomesAndContigs;
use bam_generator::*;

pub fn mosdepth_genome_coverage_with_contig_names<T: MosdepthGenomeCoverageEstimator<T> + std::fmt::Debug,
                                                  R: NamedBamReader,
                                                  G: NamedBamReaderGenerator<R>>(
    bam_readers: Vec<G>,
    contigs_and_genomes: &GenomesAndContigs,
    coverage_estimators: &mut T,
    print_zero_coverage_genomes: bool,
    flag_filtering: bool) -> Vec<OutputStream>{
    let mut output_vec: Vec<OutputStream> = Vec::new();
    for mut bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();

        let stoit_name = &(bam_generated.name().to_string());
        debug!("Working on stoit {}", stoit_name);
        let header = bam_generated.header().clone();
        let target_names = header.target_names();

        // Collect reference numbers for each genome's contigs
        let mut reference_number_to_genome_index: Vec<Option<usize>> = vec![];
        let mut num_refs_in_genomes: u32 = 0;
        let mut num_refs_not_in_genomes: u32 = 0;
        for name in target_names {
            let genome_index = contigs_and_genomes.genome_index_of_contig(
                &String::from(std::str::from_utf8(name)
                              .expect("UTF8 encoding error in BAM header file")));
            match genome_index {
                Some(i) => {
                    reference_number_to_genome_index.push(Some(i));
                    num_refs_in_genomes += 1;
                },
                None => {
                    reference_number_to_genome_index.push(None);
                    num_refs_not_in_genomes += 1;
                }
            }
        }
        info!("Of {} reference IDs, {} were assigned to a genome and {} were not",
              num_refs_in_genomes + num_refs_not_in_genomes,
              num_refs_in_genomes, num_refs_not_in_genomes);
        debug!("Reference number to genoems: {:?}", reference_number_to_genome_index);
        if num_refs_in_genomes == 0 {
            eprintln!("Error: There are no found reference sequences that are a part of a genome");
            std::process::exit(2);
        }


        let mut per_genome_coverage_estimators: Vec<T> = vec!();
        for _ in contigs_and_genomes.genomes.iter() {
            let cov_clone = coverage_estimators.copy();
            per_genome_coverage_estimators.push(cov_clone);
        }
        // Iterate through bam records
        let mut last_tid: u32 = 0;
        let mut doing_first = true;
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let mut record: bam::record::Record = bam::record::Record::new();
        let mut seen_ref_ids = BTreeSet::new();
        while bam_generated.read(&mut record).is_ok() {
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
                    match reference_number_to_genome_index[last_tid as usize] {
                        Some(genome_index) => {
                            per_genome_coverage_estimators[genome_index]
                                .add_contig(&ups_and_downs);
                        },
                        None => {}
                    }
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

        if doing_first {
            warn!("No reads were observed - perhaps something went wrong in the mapping?");
        } else {
            // Record the last contig
            match reference_number_to_genome_index[last_tid as usize] {
                Some(genome_index) => {
                    per_genome_coverage_estimators[genome_index]
                        .add_contig(&ups_and_downs);
                },
                None => {}
            }

            // Print the coverages of each genome
            // Calculate the unobserved length of each genome
            let mut unobserved_lengths: Vec<u32> = vec!();
            for _ in 0..contigs_and_genomes.genomes.len() {
                unobserved_lengths.push(0)
            }
            debug!("estimators: {:?}", per_genome_coverage_estimators);
            for (ref_id, genome_id_option) in reference_number_to_genome_index.iter().enumerate() {
                let ref_id_u32: u32 = ref_id as u32;
                debug!("Seen {:?}", seen_ref_ids);
                match genome_id_option {
                    Some(genome_id) => {
                        if !seen_ref_ids.contains(&ref_id_u32) {
                            debug!("Getting target #{} from header names", ref_id_u32);
                            unobserved_lengths[*genome_id] += header.target_len(ref_id_u32).unwrap()
                        }
                    },
                    None => {}
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
                    if per_genome_coverage_estimators[i].is_histogram() {
                        let mut hist_vec = per_genome_coverage_estimators[i].add_to_stream(output);
                        for hvec in hist_vec{
                            output_vec.push(hvec);
                        }                        
                    }else{
                        output_vec.push(output);
                    }
                } else if print_zero_coverage_genomes {
                    let mut output = OutputStream::new(
                        stoit_name.to_string(),
                        genome.to_string(),
                        0.0
                    );
                    if per_genome_coverage_estimators[i].is_histogram() {
                        let mut hist_vec = per_genome_coverage_estimators[i].add_to_stream(output);
                        for hvec in hist_vec{
                            output_vec.push(hvec);
                        }                        
                    }else{
                        output_vec.push(output);
                    }
                }
            }
        }
        bam_generated.finish();
        }
    return output_vec;
}




pub fn mosdepth_genome_coverage<T: MosdepthGenomeCoverageEstimator<T> + std::fmt::Debug,
                                R: NamedBamReader,
                                G: NamedBamReaderGenerator<R>>(
    bam_readers: Vec<G>,
    split_char: u8,
    coverage_estimator: &mut T,
    print_zero_coverage_genomes: bool,
    flag_filtering: bool,
    single_genome: bool) -> Vec<OutputStream>{

    let mut output_vec: Vec<OutputStream> = Vec::new();
    for bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();

        let stoit_name = &(bam_generated.name().to_string());
        debug!("Working on stoit {}", stoit_name);
        let header = bam_generated.header().clone();
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

        // Iterate through the bam
        let mut last_tid: u32 = 0;
        let mut doing_first = true;
        let mut last_genome: &[u8] = "error genome".as_bytes();
        let mut unobserved_contig_length: u32 = 0;
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let mut record: bam::record::Record = bam::record::Record::new();
        while bam_generated.read(&mut record).is_ok() {
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
                        let zero_output = print_previous_zero_coverage_genomes2(
                                stoit_name, b"", current_genome, tid, coverage_estimator,
                                &target_names, split_char);
                        for z in zero_output{
                            output_vec.push(z);
                        }
                        
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
                        if coverage_estimator.is_histogram() {
                            let mut hist_vec = coverage_estimator.add_to_stream(output);
                            for hvec in hist_vec{
                                output_vec.push(hvec);
                            }                        
                        }else{
                            output_vec.push(output);
                        }
                    } else if print_zero_coverage_genomes {
                        let mut output = OutputStream::new(
                            stoit_name.to_string(),
                            str::from_utf8(last_genome).unwrap().to_string(),
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
                    }
                    coverage_estimator.setup();
                    if print_zero_coverage_genomes {
                        let zero_output = print_previous_zero_coverage_genomes2(
                            stoit_name, last_genome, current_genome, tid, coverage_estimator,
                            &target_names, split_char);
                        for z in zero_output{
                            output_vec.push(z);
                        }
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

        if doing_first {
            warn!("No reads were observed - perhaps something went wrong in the mapping?");
        } else {
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
            if coverage_estimator.is_histogram(){
                let mut hist_vec = coverage_estimator.add_to_stream(output);
                for hvec in hist_vec{
                    output_vec.push(hvec);
                }  
            }else{
                output_vec.push(output);
            }
        } else if print_zero_coverage_genomes {
            let mut output = OutputStream::new(
                stoit_name.to_string(),
                str::from_utf8(last_genome).unwrap().to_string(),
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
        }
        if print_zero_coverage_genomes && !single_genome {
            let zero_output = print_previous_zero_coverage_genomes2(
                    stoit_name, last_genome, b"", header.target_count()-1, coverage_estimator,
                    &target_names, split_char);
            for z in zero_output{
                output_vec.push(z);
            }
        }
        bam_generated.finish();
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
    pileup_coverage_estimator: &'a mut MosdepthGenomeCoverageEstimator<T>,
    target_names: &Vec<&[u8]>,
    split_char: u8)
    -> Vec<OutputStream> {
    
    let mut output_vec = Vec::new();
    let mut my_current_genome = current_genome;
    let mut tid = current_tid;
    let mut genomes_to_print: Vec<&[u8]> = vec![];

    while tid > 0 {
        let genome = extract_genome(tid, &target_names, split_char);
        if genome == last_genome { break; }
        else if genome != my_current_genome {
            // In-between genome encountered for the first time.
            debug!("genome {:?} current genome {:?}", genome, my_current_genome);
            genomes_to_print.push(genome);
            my_current_genome = genome;
        }
        tid = tid - 1;
    };
    for i in (0..genomes_to_print.len()).rev() {
        let mut output = OutputStream::new(
            stoit_name.to_string(),
            str::from_utf8(genomes_to_print[i]).unwrap().to_string(),
            0.0
        );
        output_vec.push(output);
    }
    return output_vec;
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
    use env_logger;
    use std::io::BufRead;

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
    pub fn get_output(stdin: std::io::Stdin)-> String {
        let mut out_string = Vec::new();
        let lock = stdin.lock();
        for line in lock.lines(){
            out_string.push(line.unwrap());
        }
        return out_string.join("\n")
    }

    #[test]
    fn initialize_logger() {
        env_logger::init().unwrap();
    }

    #[test]
    fn test_one_genome_two_contigs_first_covered(){
        let stream =mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            'q' as u8,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true,
            false,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1\tse\t0.6\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_one_genome_two_contigs_first_covered_contig_names(){
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("se".to_string());
        geco.insert("seq1".to_string(),genome1);
        geco.insert("seq2".to_string(),genome1);
        let stream = mosdepth_genome_coverage_with_contig_names(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            &geco,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1\tse\t0.6\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_one_genome_two_contigs_second_covered(){
        let stream = mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq2.bam"]),
            'q' as u8,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true,
            false,
            false);
        assert_eq!(
            "2seqs.reads_for_seq2\tse\t0.6\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_one_genome_two_contigs_second_covered_contig_names(){
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("se".to_string());
        geco.insert("seq1".to_string(),genome1);
        geco.insert("seq2".to_string(),genome1);
        let stream = mosdepth_genome_coverage_with_contig_names(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq2.bam"]),
            &geco,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true,
            false);
        assert_eq!(
            "2seqs.reads_for_seq2\tse\t0.6\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_one_genome_two_contigs_both_covered(){
        let stream = mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            'e' as u8,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true,
            false,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_one_genome_two_contigs_both_covered_contig_names(){
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(),genome1);
        geco.insert("seq2".to_string(),genome1);
        let stream = mosdepth_genome_coverage_with_contig_names(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_one_genome_min_fraction_covered_under_min(){
        let stream = mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            'e' as u8,
            &mut MeanGenomeCoverageEstimator::new(0.76),
            true,
            false,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t0\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_one_genome_min_fraction_covered_under_min_contig_names(){
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(),genome1);
        geco.insert("seq2".to_string(),genome1);
        let stream = mosdepth_genome_coverage_with_contig_names(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            &mut MeanGenomeCoverageEstimator::new(0.76),
            false,
            false);
        assert_eq!(
            "",
            print_output_stream(stream))
    }

    #[test]
    fn test_one_genome_min_fraction_covered_just_ok(){
        let stream = mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            'e' as u8,
            &mut MeanGenomeCoverageEstimator::new(0.759),
            true,
            false,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            print_output_stream(stream))
    }


    #[test]
    fn test_one_genome_min_fraction_covered_just_ok_contig_names(){
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(),genome1);
        geco.insert("seq2".to_string(),genome1);
        let stream = mosdepth_genome_coverage_with_contig_names(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            &mut MeanGenomeCoverageEstimator::new(0.759),
            true,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_two_contigs_trimmed_mean(){
        let stream = mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            'e' as u8,
            &mut TrimmedMeanGenomeCoverageEstimator::new(0.1, 0.9, 0.0),
            true,
            false,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.08875\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_two_contigs_trimmed_mean_contig_names(){
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(),genome1);
        geco.insert("seq2".to_string(),genome1);
        let stream = mosdepth_genome_coverage_with_contig_names(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            &mut TrimmedMeanGenomeCoverageEstimator::new(0.1, 0.9, 0.0),
            true,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.08875\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_two_contigs_pileup_counts_estimator(){
        let stream = mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            'e' as u8,
            &mut PileupCountsGenomeCoverageEstimator::new(0.0),
            true,
            false,
            false);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t0\t482\n2seqs.reads_for_seq1_and_seq2\ts\t1\t922\n2seqs.reads_for_seq1_and_seq2\ts\t2\t371\n2seqs.reads_for_seq1_and_seq2\ts\t3\t164\n2seqs.reads_for_seq1_and_seq2\ts\t4\t61\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_two_contigs_pileup_counts_estimator_contig_names(){

        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(),genome1);
        geco.insert("seq2".to_string(),genome1);
        let stream = mosdepth_genome_coverage_with_contig_names(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            &mut PileupCountsGenomeCoverageEstimator::new(0.0),
            true,
            false);
        
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t0\t482\n2seqs.reads_for_seq1_and_seq2\ts\t1\t922\n2seqs.reads_for_seq1_and_seq2\ts\t2\t371\n2seqs.reads_for_seq1_and_seq2\ts\t3\t164\n2seqs.reads_for_seq1_and_seq2\ts\t4\t61\n",
            print_output_stream(stream))
    }

    #[test]
    fn test_zero_coverage_genomes(){
        let stream1 = mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            '~' as u8,
            &mut MeanGenomeCoverageEstimator::new(0.1),
            true,
            false,
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0\n7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome3\t0\n7seqs.reads_for_seq1_and_seq2\tgenome4\t0\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6\t0\n",
            print_output_stream(stream1));

        let stream2 = mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            '~' as u8,
            &mut MeanGenomeCoverageEstimator::new(0.1),
            false,
            false,
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n",
            print_output_stream(stream2))
    }

    #[test]
    fn test_zero_coverage_genomes_after_min_fraction(){
        let stream = mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            '~' as u8,
            &mut MeanGenomeCoverageEstimator::new(0.76),
            true,
            false,
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0\n7seqs.reads_for_seq1_and_seq2\tgenome2\t0\n7seqs.reads_for_seq1_and_seq2\tgenome3\t0\n7seqs.reads_for_seq1_and_seq2\tgenome4\t0\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6\t0\n",
            print_output_stream(stream));
    }

    #[test]
    fn test_single_genome(){
        let stream = mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            '~' as u8,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true,
            false,
            true);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0.04209345\n",
            print_output_stream(stream));
    }

    #[test]
    fn test_zero_coverage_genomes_contig_names(){
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
        let stream1 = mosdepth_genome_coverage_with_contig_names(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            &mut MeanGenomeCoverageEstimator::new(0.1),
            true,
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0\n7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome3\t0\n7seqs.reads_for_seq1_and_seq2\tgenome4\t0\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6\t0\n",
            print_output_stream(stream1));

        let stream2 = mosdepth_genome_coverage_with_contig_names(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            &mut MeanGenomeCoverageEstimator::new(0.1),
            false,
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n",
            print_output_stream(stream2))
    }
}
