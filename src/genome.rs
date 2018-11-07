
use std;
use std::any::Any;
use std::result::Result;
use rust_htslib::bam;

use std::process;
use std::str;
use std::collections::BTreeSet;

use mosdepth_genome_coverage_estimators::*;
use genomes_and_contigs::GenomesAndContigs;
use bam_generator::*;
use get_trimmed_mean_estimator;


pub fn mosdepth_genome_coverage_with_contig_names<R: NamedBamReader,
                                                  G: NamedBamReaderGenerator<R>>(
    bam_readers: Vec<G>,
    contigs_and_genomes: &GenomesAndContigs,
    min_fraction_covered: f32,
    print_stream: &mut std::io::Write,
    print_zero_coverage_genomes: bool,
    flag_filtering: bool,
    m: Vec<f32>,
    methods: Vec<&str>) {
    let mut coverage_estimator_box: Vec<CoverageEstimator> = Vec::new();
    let mut coverage_estimator_vec: Vec<CoverageEstimator>;
    for method in methods {
        match method {
            "mean" =>{
                coverage_estimator_box.push(CoverageEstimator::new("mean", 0.0, 0.0,
                                                                   min_fraction_covered));
            },
            "coverage_histogram" =>{
                coverage_estimator_box.push(CoverageEstimator::new("coverage_histogram", 0.0, 0.0,
                                                                   min_fraction_covered));
            },
            "trimmed_mean" =>{
                coverage_estimator_box.push(get_trimmed_mean_estimator(m.clone(), min_fraction_covered));
            },
            "covered_fraction" =>{
                coverage_estimator_box.push(CoverageEstimator::new("covered_fraction", 0.0, 0.0,
                                                                   min_fraction_covered));
            },
            "variance" =>{
                coverage_estimator_box.push(CoverageEstimator::new("variance", 0.0, 0.0,
                                                                   min_fraction_covered));
            },
            _ => panic!("programming error")
        }
    }
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
        let mut per_genome_coverage_estimators = vec!();
        for _ in contigs_and_genomes.genomes.iter() {
            let cov_clone = coverage_estimator_box.clone();
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
                            coverage_estimator_vec = Vec::new();
                            for mut coverage_estimator in per_genome_coverage_estimators[genome_index].clone() {
                                coverage_estimator = coverage_estimator.add_contig(&ups_and_downs);
                                coverage_estimator_vec.push(coverage_estimator);
                            };
                            per_genome_coverage_estimators[genome_index] = coverage_estimator_vec;
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
                    'M' | 'X' | '=' => {
                        // if M, X, or =, increment start and decrement end index
                        debug!("Adding M, X, or =, at {} and {}", cursor, cursor + cig.len() as usize);
                        ups_and_downs[cursor] += 1;
                        let final_pos = cursor + cig.len() as usize;
                        if final_pos < ups_and_downs.len() { // True unless the read hits the contig end.
                            ups_and_downs[final_pos] -= 1;
                        }
                        cursor += cig.len() as usize;
                    },
                    'D' | 'N' => {
                        // if D, move the cursor
                        cursor += cig.len() as usize;
                    },
                    'I' | 'S' | 'H' | 'P' => {},
                    _ => panic!("Unknown CIGAR string match")
                }
            }
        }

        if doing_first {
            warn!("No reads were observed - perhaps something went wrong in the mapping?");
        } else {
            // Record the last contig
            match reference_number_to_genome_index[last_tid as usize] {
                Some(genome_index) => {
                    coverage_estimator_vec = Vec::new();
                    for mut coverage_estimator in per_genome_coverage_estimators[genome_index].clone() {
                        coverage_estimator = coverage_estimator.add_contig(&ups_and_downs);
                        coverage_estimator_vec.push(coverage_estimator);
                    };
                    per_genome_coverage_estimators[genome_index] = coverage_estimator_vec;
                },
                None => {}
            }

            // Print the coverages of each genome
            // Calculate the unobserved length of each genome
            let mut unobserved_lengths: Vec<u32> = vec!();
            for _ in 0..contigs_and_genomes.genomes.len() {
                unobserved_lengths.push(0)
            }
//                debug!("estimators: {:?}", per_genome_coverage_estimators);
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

                for (j, mut coverage_estimator) in per_genome_coverage_estimators[i].clone().iter().enumerate() {
                    let coverage = coverage_estimator.clone()
                        .calculate_coverage(unobserved_lengths[i]);

                    // Print coverage of previous genome
                    debug!("Found coverage {} for genome {}", coverage, genome);
                    if coverage > 0.0 {
                        if j == 0 as usize{
                            print_genome(&stoit_name, &genome, print_stream);
                        }
                        coverage_estimator.print_coverage(
                            &stoit_name,
                            &genome,
                            &coverage,
                            print_stream);
                        if j+1 == per_genome_coverage_estimators[i].clone().iter().len(){
                            write!(print_stream, "\n");
                        }
                    } else if print_zero_coverage_genomes {
                        if j == 0 as usize{
                            print_genome(&stoit_name, &genome, print_stream);
                        }
                        coverage_estimator.print_zero_coverage(
                            print_stream);
                        if j+1 == per_genome_coverage_estimators[i].clone().iter().len(){
                            write!(print_stream, "\n");
                        }
                    }
                }

            }
        }

//        write!(print_stream, "\n").unwrap();
        bam_generated.finish();
    }
}


pub fn mosdepth_genome_coverage<R: NamedBamReader,
                                G: NamedBamReaderGenerator<R>>(
    bam_readers: Vec<G>,
    split_char: u8,
    print_stream: &mut std::io::Write,
    min_fraction_covered: f32,
    print_zero_coverage_genomes: bool,
    methods: Vec<&str>,
    m: Vec<f32>,
    flag_filtering: bool,
    single_genome: bool) {
    let mut coverage_estimator_box: Vec<CoverageEstimator> = Vec::new();
    let mut coverage_estimator_vec: Vec<CoverageEstimator>;
    for method in methods {
        match method {
            "mean" =>{
                coverage_estimator_box.push(CoverageEstimator::new("mean", 0.0, 0.0,
                                                                   min_fraction_covered));
            },
            "coverage_histogram" =>{
                coverage_estimator_box.push(CoverageEstimator::new("coverage_histogram", 0.0, 0.0,
                                                                   min_fraction_covered));
            },
            "trimmed_mean" =>{
                coverage_estimator_box.push(get_trimmed_mean_estimator(m.clone(), min_fraction_covered));
            },
            "covered_fraction" =>{
                coverage_estimator_box.push(CoverageEstimator::new("covered_fraction", 0.0, 0.0,
                                                                   min_fraction_covered));
            },
            "variance" =>{
                coverage_estimator_box.push(CoverageEstimator::new("variance", 0.0, 0.0,
                                                                   min_fraction_covered));
            },
            _ => panic!("programming error")
        }
    }
    for mut bam_generator in bam_readers {
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
//            let mut coverage_estimator;
            if tid != last_tid || doing_first {
                if doing_first == true {
                    coverage_estimator_vec = Vec::new();
                    for (i, mut coverage_estimator) in coverage_estimator_box.clone().iter().enumerate() {
//                        coverage_estimator = coverage_estimator_b;
                        let mut coverage_estimator = coverage_estimator.clone().setup();
                        coverage_estimator_vec.push(coverage_estimator.clone());
                    }
                    last_genome = current_genome;
                    unobserved_contig_length = fill_genome_length_backwards(tid, current_genome);
                    coverage_estimator_box = coverage_estimator_vec;
                    doing_first = false;

                    if print_zero_coverage_genomes && !single_genome {
                        print_previous_zero_coverage_genomes2(
                            stoit_name, b"", current_genome, tid, &coverage_estimator_box,
                            &target_names, split_char, print_stream);
                    }

                } else if current_genome == last_genome {
                    coverage_estimator_vec = Vec::new();
                    for mut coverage_estimator in coverage_estimator_box.clone() {
                        coverage_estimator = coverage_estimator.add_contig(&ups_and_downs);
                        coverage_estimator_vec.push(coverage_estimator.clone());
                        // Collect the length of reference sequences from this
                        // genome that had no hits that were just skipped over.
                        debug!("Filling unobserved from {} to {}", last_tid, tid);
                        unobserved_contig_length += fill_genome_length_backwards_to_last(
                            tid, last_tid as u32, current_genome);
                    }
                    coverage_estimator_box = coverage_estimator_vec;
                } else {
//                    print_genome(stoit_name, &str::from_utf8(current_genome).unwrap(), print_stream);
                    coverage_estimator_vec = Vec::new();
                    for (i, mut coverage_estimator) in coverage_estimator_box.clone().iter().enumerate(){
                        let mut coverage_estimator = coverage_estimator.clone().add_contig(&ups_and_downs);
                        coverage_estimator_vec.push(coverage_estimator.clone());
                        // Collect the length of refs from the end of the last genome that had no hits
                        debug!("Filling unobserved from {} to {} for {}", last_tid, tid, &str::from_utf8(last_genome).unwrap());
                        unobserved_contig_length += fill_genome_length_backwards_to_last(
                            tid, last_tid as u32, last_genome);
                        debug!("unobserved_contig_length now {}", unobserved_contig_length);
                        // Determine coverage of previous genome
                        let coverage = coverage_estimator.calculate_coverage(unobserved_contig_length);

                        // Print coverage of previous genome
                        if coverage > 0.0 {
                            if i == 0{
                                print_genome(stoit_name, &str::from_utf8(last_genome).unwrap(), print_stream);
                            }
                            coverage_estimator.print_coverage(
                                &stoit_name,
                                &str::from_utf8(last_genome).unwrap(),
                                &coverage,
                                print_stream);
                            if i+1 == coverage_estimator_box.clone().iter().len(){
                                write!(print_stream, "\n");
                            }
                        } else if print_zero_coverage_genomes {
                            if i == 0{
                                print_genome(stoit_name, &str::from_utf8(last_genome).unwrap(), print_stream);
                            }
                            coverage_estimator.print_zero_coverage(
                                print_stream);
                            if i+1 == coverage_estimator_box.clone().iter().len(){
                                write!(print_stream, "\n");
                            }
                        }
                        coverage_estimator.setup();

                    }
                    if print_zero_coverage_genomes {
                        print_previous_zero_coverage_genomes2(
                            stoit_name, last_genome, current_genome, tid, &coverage_estimator_vec,
                            &target_names, split_char, print_stream);
                    }
                    coverage_estimator_box = coverage_estimator_vec;
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
                    'M' | 'X' | '=' => {
                        // if M, X, or =, increment start and decrement end index
                        debug!("Adding M, X, or =, at {} and {}", cursor, cursor + cig.len() as usize);
                        ups_and_downs[cursor] += 1;
                        let final_pos = cursor + cig.len() as usize;
                        if final_pos < ups_and_downs.len() { // True unless the read hits the contig end.
                            ups_and_downs[final_pos] -= 1;
                        }
                        cursor += cig.len() as usize;
                    },
                    'D' | 'N' => {
                        // if D or N, move the cursor
                        cursor += cig.len() as usize;
                    },
                    'I' | 'S' | 'H' | 'P' => {},
                    _ => panic!("Unknown CIGAR string match")
                }
            }
        }

        if doing_first {
            warn!("No reads were observed - perhaps something went wrong in the mapping?");
        } else {
            // Print the last genome
            // Give the single genome a dummy name
            if single_genome {
                last_genome = "genome1".as_bytes()
            }

            // Collect the length of refs from the end of the last genome that had no hits
            unobserved_contig_length += fill_genome_length_forwards(last_tid, last_genome);
            debug!("At end, unobserved_contig_length now {}", unobserved_contig_length);
            // Determine coverage of previous genome
            coverage_estimator_vec = Vec::new();
            for (i, mut coverage_estimator)in coverage_estimator_box.clone().iter().enumerate(){
                let mut coverage_estimator = coverage_estimator.clone().add_contig(&ups_and_downs);
                coverage_estimator_vec.push(coverage_estimator.clone());
                let coverage = coverage_estimator.calculate_coverage(unobserved_contig_length);

                // Print coverage of previous genome
                if coverage > 0.0 {
                    if i == 0{
                        print_genome(stoit_name, &str::from_utf8(last_genome).unwrap(), print_stream);
                    }
                    coverage_estimator.print_coverage(
                        &stoit_name,
                        &str::from_utf8(last_genome).unwrap(),
                        &coverage,
                        print_stream);
                    if i+1 == coverage_estimator_box.clone().iter().len(){
                        write!(print_stream, "\n");
                    }
                } else if print_zero_coverage_genomes {
                    if i == 0{
                        print_genome(stoit_name, &str::from_utf8(last_genome).unwrap(), print_stream);
                    }

                    coverage_estimator.print_zero_coverage(
                        print_stream);

                    if i+1 == coverage_estimator_box.clone().iter().len(){
                        write!(print_stream, "\n");
                    }
                }

            }
            if print_zero_coverage_genomes && !single_genome {
                print_previous_zero_coverage_genomes2(
                    stoit_name, last_genome, b"", header.target_count() - 1, &coverage_estimator_vec,
                    &target_names, split_char, print_stream);
            }
            coverage_estimator_box = coverage_estimator_vec;
        }

        bam_generated.finish();
    }
}

fn print_genome<'a >(stoit_name: &str,
                     genome: &str,
                     print_stream: &'a mut std::io::Write) -> &'a mut std::io::Write {
    write!(print_stream, "{}\t{}",
           stoit_name,
           genome).unwrap();
    return print_stream;
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
    pileup_coverage_estimator: &'a Vec<CoverageEstimator>,
    target_names: &Vec<&[u8]>,
    split_char: u8,
    print_stream: &mut std::io::Write)
    -> &'a Vec<CoverageEstimator>{

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
        print_genome(stoit_name, &str::from_utf8(genomes_to_print[i]).unwrap(), print_stream);

        for mut coverage_estimator in pileup_coverage_estimator.clone() {
            coverage_estimator.print_zero_coverage(
                print_stream);
        }
        write!(print_stream, "\n");
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
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            'q' as u8,
            &mut stream,
            0.0,
            true,
            vec!("mean"),
            vec!(0.1, 0.9),
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
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            &geco,
            0.0,
            &mut stream,
            true,
            false,
            vec!(0.1, 0.9),
            vec!("mean"));
        assert_eq!(
            "2seqs.reads_for_seq1\tse\t0.6\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_two_contigs_second_covered(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq2.bam"]),
            'q' as u8,
            &mut stream,
            0.0,
            true,
            vec!("mean"),
            vec!(0.1, 0.9),
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
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq2.bam"]),
            &geco,
            0.0,
            &mut stream,
            true,
            false,
            vec!(0.1, 0.9),
            vec!("mean"));
        assert_eq!(
            "2seqs.reads_for_seq2\tse\t0.6\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_two_contigs_both_covered(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            'e' as u8,
            &mut stream,
            0.0,
            true,
            vec!("mean"),
            vec!(0.1, 0.9),
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
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            0.0,
            &mut stream,
            true,
            false,
            vec!(0.1, 0.9),
            vec!("mean"));
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_min_fraction_covered_under_min(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            'e' as u8,
            &mut stream,
            0.76,
            true,
            vec!("mean"),
            vec!(0.1, 0.9),
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
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            0.76,
            &mut stream,
            false,
            false,
            vec!(0.1, 0.9),
            vec!("mean"));
        assert_eq!(
            "",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_min_fraction_covered_just_ok(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            'e' as u8,
            &mut stream,
            0.759,
            true,
            vec!("mean"),
            vec!(0.1, 0.9),
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
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            0.759,
            &mut stream,
            true,
            false,
            vec!(0.1, 0.9),
            vec!("mean"));
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_two_contigs_trimmed_mean(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            'e' as u8,
            &mut stream,
            0.759,
            true,
            vec!("trimmed_mean"),
            vec!(0.1, 0.9),
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
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            0.0,
            &mut stream,
            true,
            false,
            vec!(0.1, 0.9),
            vec!("trimmed_mean"));
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.08875\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_two_contigs_pileup_counts_estimator(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            'e' as u8,
            &mut stream,
            0.0,
            true,
            vec!("coverage_histogram"),
            vec!(0.1, 0.9),
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
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            0.0,
            &mut stream,
            true,
            false,
            vec!(0.1, 0.9),
            vec!("coverage_histogram"));
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t0\t482\n2seqs.reads_for_seq1_and_seq2\ts\t1\t922\n2seqs.reads_for_seq1_and_seq2\ts\t2\t371\n2seqs.reads_for_seq1_and_seq2\ts\t3\t164\n2seqs.reads_for_seq1_and_seq2\ts\t4\t61\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_zero_coverage_genomes(){
        let mut stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            '~' as u8,
            &mut stream,
            0.1,
            true,
            vec!("mean"),
            vec!(0.1, 0.9),
            false,
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome3\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome4\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6\t0.0\n",
            str::from_utf8(stream.get_ref()).unwrap());

        stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            '~' as u8,
            &mut stream,
            0.1,
            false,
            vec!("mean"),
            vec!(0.1, 0.9),
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
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            '~' as u8,
            &mut stream,
            0.759,
            true,
            vec!("mean"),
            vec!(0.1, 0.9),
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
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            '~' as u8,
            &mut stream,
            0.0,
            true,
            vec!("mean"),
            vec!(0.1, 0.9),
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
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            0.1,
            &mut stream,
            true,
            false,
            vec!(0.1, 0.9),
            vec!("mean"));
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0.0\n\
            7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n\
            7seqs.reads_for_seq1_and_seq2\tgenome3\t0.0\n\
            7seqs.reads_for_seq1_and_seq2\tgenome4\t0.0\n\
            7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n\
            7seqs.reads_for_seq1_and_seq2\tgenome6\t0.0\n",
            str::from_utf8(stream.get_ref()).unwrap());

        stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage_with_contig_names(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            0.1,
            &mut stream,
            false,
            false,
            vec!(0.1, 0.9),
            vec!("mean"));
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_zero_coverage_genomes_contig_names_with_multiple_methods(){
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
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            0.1,
            &mut stream,
            true,
            false,
            vec!(0.1, 0.9),
            vec!("mean", "variance"));
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0.0\t0.0\n\
            7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\t1.3633634\n\
            7seqs.reads_for_seq1_and_seq2\tgenome3\t0.0\t0.0\n\
            7seqs.reads_for_seq1_and_seq2\tgenome4\t0.0\t0.0\n\
            7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\t0.6166166\n\
            7seqs.reads_for_seq1_and_seq2\tgenome6\t0.0\t0.0\n",
            str::from_utf8(stream.get_ref()).unwrap());

        stream = Cursor::new(Vec::new());
        mosdepth_genome_coverage_with_contig_names(
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            0.1,
            &mut stream,
            false,
            false,
            vec!(0.1, 0.9),
            vec!("mean","variance"));
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\t1.3633634\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\t0.6166166\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }
}
