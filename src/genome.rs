use nm;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use std;
use std::process;
use FlagFilter;

use std::collections::BTreeSet;
use std::str;

use bam_generator::*;
use coverage_takers::*;
use genomes_and_contigs::find_first;
use genomes_and_contigs::GenomesAndContigs;
use mosdepth_genome_coverage_estimators::*;
use ReadsMapped;

pub fn mosdepth_genome_coverage_with_contig_names<
    R: NamedBamReader,
    G: NamedBamReaderGenerator<R>,
    T: CoverageTaker,
>(
    bam_readers: Vec<G>,
    contigs_and_genomes: &GenomesAndContigs,
    coverage_taker: &mut T,
    print_zero_coverage_genomes: bool,
    flag_filters: &FlagFilter,
    coverage_estimators: &mut [CoverageEstimator],
    threads: u16,
) -> Vec<ReadsMapped> {
    let mut reads_mapped_vector = vec![];
    let mut is_first_bam = true;
    for bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();
        bam_generated.set_threads(threads as usize);

        let stoit_name = &(bam_generated.name().to_string());
        debug!("Working on stoit {}", stoit_name);
        coverage_taker.start_stoit(stoit_name);
        let header = bam_generated.header().clone();
        let target_names = header.target_names();

        // Collect reference numbers for each genome's contigs
        let mut reference_number_to_genome_index: Vec<Option<usize>> = vec![];
        let mut num_refs_in_genomes: u32 = 0;
        let mut num_refs_not_in_genomes: u32 = 0;
        // Collect reference numbers for each genome
        let mut genome_index_to_references: Vec<Vec<u32>> =
            vec![vec!(); contigs_and_genomes.genomes.len()];
        // Reads mapped are only counted when the genome has non-zero coverage.
        let mut reads_mapped_in_each_genome: Vec<u64> = vec![0; contigs_and_genomes.genomes.len()];
        for (tid, name) in target_names.iter().enumerate() {
            let genome_index = contigs_and_genomes.genome_index_of_contig(&String::from(
                std::str::from_utf8(name).expect("UTF8 encoding error in BAM header file"),
            ));

            match genome_index {
                Some(i) => {
                    reference_number_to_genome_index.push(Some(i));
                    num_refs_in_genomes += 1;
                    genome_index_to_references[i].push(tid as u32);
                }
                None => {
                    reference_number_to_genome_index.push(None);
                    num_refs_not_in_genomes += 1;
                }
            }
        }
        if is_first_bam {
            is_first_bam = false;
            info!(
                "Of {} reference IDs, {} were assigned to a genome and {} were not",
                num_refs_in_genomes + num_refs_not_in_genomes,
                num_refs_in_genomes,
                num_refs_not_in_genomes
            );
        }
        trace!(
            "Reference number to genomes: {:?}",
            reference_number_to_genome_index
        );
        if num_refs_in_genomes == 0 {
            error!("Error: There are no found reference sequences that are a part of a genome");
            process::exit(1);
        }
        {
            let num_unreferenced =
                contigs_and_genomes.contig_to_genome.len() as u32 - num_refs_in_genomes;
            if num_unreferenced > 0 {
                warn!(
                    "Found {} contig(s) that were defined as being part of a genome, \
                       but were not reference sequences in BAM files.",
                    num_unreferenced
                )
            }
        }
        let mut per_genome_coverage_estimators = vec![];
        for _ in contigs_and_genomes.genomes.iter() {
            let mut cov_clone = vec![coverage_estimators[0].clone(); coverage_estimators.len()];
            cov_clone.clone_from_slice(coverage_estimators);
            per_genome_coverage_estimators.push(cov_clone);
        }

        // Iterate through bam records
        let mut last_tid: u32 = 0;
        let mut doing_first = true;
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let mut record: bam::record::Record = bam::record::Record::new();
        let mut seen_ref_ids = BTreeSet::new();
        let mut num_mapped_reads_in_current_contig: u64 = 0;
        let mut total_edit_distance_in_current_contig: u64 = 0;
        let mut total_indels_in_current_contig: u64 = 0;
        loop {
            match bam_generated.read(&mut record) {
                None => {
                    break;
                }
                Some(Ok(())) => {}
                Some(e) => {
                    panic!("Error reading BAM record: {:?}", e)
                }
            }

            if !flag_filters.passes(&record) {
                trace!("Skipping read based on flag filtering");
                continue;
            }
            let original_tid = record.tid();
            if !record.is_unmapped() {
                // if mapped
                let tid = original_tid as u32;
                if tid != last_tid || doing_first {
                    debug!("Came across a new tid {}", tid);
                    if doing_first {
                        doing_first = false;
                    } else {
                        if tid < last_tid {
                            error!("BAM file appears to be unsorted. Input BAM files must be sorted by reference (i.e. by samtools sort)");
                            panic!("BAM file appears to be unsorted. Input BAM files must be sorted by reference (i.e. by samtools sort)");
                        }
                        if let Some(genome_index) =
                            reference_number_to_genome_index[last_tid as usize]
                        {
                            debug!(
                                "Found {} reads mapped to tid {}",
                                num_mapped_reads_in_current_contig, last_tid
                            );
                            for ref mut coverage_estimator in
                                per_genome_coverage_estimators[genome_index].iter_mut()
                            {
                                coverage_estimator.add_contig(
                                    &ups_and_downs,
                                    num_mapped_reads_in_current_contig,
                                    total_edit_distance_in_current_contig
                                        - total_indels_in_current_contig,
                                );
                            }
                        }
                    }

                    ups_and_downs =
                        vec![0; header.target_len(tid).expect("Corrupt BAM file?") as usize];
                    num_mapped_reads_in_current_contig = 0;
                    total_edit_distance_in_current_contig = 0;
                    total_indels_in_current_contig = 0;
                    last_tid = tid;
                    seen_ref_ids.insert(tid);
                }

                // TODO: move below into a function for code-reuse purposes.
                // Add coverage info for the current record
                // for each chunk of the cigar string
                match reference_number_to_genome_index[tid as usize] {
                    None => {}
                    Some(genome_index) => {
                        reads_mapped_in_each_genome[genome_index] += 1;
                        num_mapped_reads_in_current_contig += 1;
                        trace!(
                            "read name {:?}",
                            std::str::from_utf8(record.qname()).unwrap()
                        );
                        let mut cursor: usize = record.pos() as usize;
                        for cig in record.cigar().iter() {
                            trace!("Found cigar {:} from {}", cig, cursor);
                            match cig {
                                Cigar::Match(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                                    // if M, X, or =, increment start and decrement end index
                                    trace!(
                                        "Adding M, X, or =, at {} and {}",
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
                        total_edit_distance_in_current_contig += nm(&record);
                    }
                }
            }
        }

        let mut num_mapped_reads_total: u64 = 0;
        if doing_first && bam_generated.num_detected_primary_alignments() == 0 {
            warn!(
                "No primary alignments were observed for sample {} \
                   - perhaps something went wrong in the mapping?",
                stoit_name
            );
        } else {
            // Record the last contig
            if let Some(genome_index) = reference_number_to_genome_index[last_tid as usize] {
                for ref mut coverage_estimator in
                    per_genome_coverage_estimators[genome_index].iter_mut()
                {
                    coverage_estimator.add_contig(
                        &ups_and_downs,
                        num_mapped_reads_in_current_contig,
                        total_edit_distance_in_current_contig - total_indels_in_current_contig,
                    )
                }
            }

            // Print the coverages of each genome
            // Calculate the unobserved lengths of each genome's contigs
            let mut unobserved_lengths: Vec<Vec<u64>> = vec![];
            for _ in 0..contigs_and_genomes.genomes.len() {
                unobserved_lengths.push(vec![])
            }
            for (ref_id, genome_id_option) in reference_number_to_genome_index.iter().enumerate() {
                let ref_id_u32: u32 = ref_id as u32;
                trace!("Seen {:?}", seen_ref_ids);
                match genome_id_option {
                    Some(genome_id) => {
                        if !seen_ref_ids.contains(&ref_id_u32) {
                            debug!("Getting target #{} from header names", ref_id_u32);
                            unobserved_lengths[*genome_id]
                                .push(header.target_len(ref_id_u32).unwrap())
                        }
                    }
                    None => {}
                }
            }
            // print the genomes out
            for (i, genome) in contigs_and_genomes.genomes.iter().enumerate() {
                // Determine if any coverages are non-zero.
                let coverages: Vec<f32> = per_genome_coverage_estimators[i]
                    .iter_mut()
                    .map(|coverage_estimator| {
                        coverage_estimator.calculate_coverage(&unobserved_lengths[i])
                    })
                    .collect();
                let any_nonzero_coverage = coverages.iter().any(|c| *c > 0.0);
                if any_nonzero_coverage {
                    num_mapped_reads_total += reads_mapped_in_each_genome[i];
                }
                if print_zero_coverage_genomes || any_nonzero_coverage {
                    coverage_taker.start_entry(i, genome);
                    for (j, ref mut coverage_estimator) in
                        per_genome_coverage_estimators[i].iter_mut().enumerate()
                    {
                        let coverage = coverages[j];

                        // Print coverage of previous genome
                        debug!("Found coverage {} for genome {}", coverage, genome);
                        if coverage > 0.0 {
                            coverage_estimator.print_coverage(&coverage, coverage_taker);
                        } else {
                            coverage_estimator.print_zero_coverage(
                                coverage_taker,
                                genome_index_to_references[i]
                                    .iter()
                                    .map(|tid| header.target_len(*tid).unwrap())
                                    .sum(),
                            );
                        }
                    }
                    coverage_taker.finish_entry();
                }
            }
        }

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

        bam_generated.finish();
    }
    reads_mapped_vector
}

#[derive(PartialEq, Debug)]
struct UnobservedLengthAndFirstTid {
    unobserved_contig_lengths: Vec<u64>,
    first_tid: usize,
}

#[allow(clippy::too_many_arguments)]
fn print_last_genomes<T: CoverageTaker>(
    num_mapped_reads_in_current_contig: u64,
    last_genome: Option<&[u8]>,
    unobserved_contig_length_and_first_tid: &mut UnobservedLengthAndFirstTid,
    ups_and_downs: &[i32],
    total_edit_distance_in_current_contig: u64,
    total_indels_in_current_contig: u64,
    current_genome: &[u8],
    coverage_estimators: &mut Vec<CoverageEstimator>,
    coverage_taker: &mut T,
    print_zero_coverage_genomes: bool,
    single_genome: bool,
    target_names: &[&[u8]],
    split_char: u8,
    header: &rust_htslib::bam::HeaderView,
    tid_to_print_zeros_to: u32,
) -> bool {
    //    debug!("ups_and_downs {:?}", &ups_and_downs);
    for coverage_estimator in coverage_estimators.iter_mut() {
        coverage_estimator.add_contig(
            ups_and_downs,
            num_mapped_reads_in_current_contig,
            total_edit_distance_in_current_contig - total_indels_in_current_contig,
        );
    }

    // Determine coverage of previous genome
    let coverages: Vec<f32> = coverage_estimators
        .iter_mut()
        .map(|coverage_estimator| {
            coverage_estimator.calculate_coverage(
                &unobserved_contig_length_and_first_tid.unobserved_contig_lengths,
            )
        })
        .collect();
    let positive_coverage = coverages.iter().any(|c| *c > 0.0);

    if print_zero_coverage_genomes || positive_coverage {
        if let Some(last_genome_name) = last_genome {
            coverage_taker.start_entry(
                unobserved_contig_length_and_first_tid.first_tid,
                str::from_utf8(last_genome.unwrap()).unwrap(),
            );
            for (i, ref mut coverage_estimator) in coverage_estimators.iter_mut().enumerate() {
                let coverage = coverages[i];
                debug!(
                    "Found coverage {} for genome {} i.e. coverage estimator {}",
                    coverage,
                    &str::from_utf8(last_genome_name).unwrap(),
                    i
                );

                // Print coverage of previous genome
                if coverage > 0.0 {
                    coverage_estimator.print_coverage(&coverage, coverage_taker);
                } else {
                    coverage_estimator.print_zero_coverage(
                        coverage_taker,
                        // Length coverage is always >0, so this 0 is never used
                        9,
                    );
                }
            }
            coverage_taker.finish_entry();
        }
    }
    // Reset all estimators
    for coverage_estimator in coverage_estimators.iter_mut() {
        coverage_estimator.setup();
    }
    if print_zero_coverage_genomes && !single_genome {
        print_previous_zero_coverage_genomes2(
            last_genome,
            current_genome,
            tid_to_print_zeros_to,
            coverage_estimators,
            target_names,
            split_char,
            coverage_taker,
            header,
        );
    }
    positive_coverage
}

#[allow(clippy::too_many_arguments)]
pub fn mosdepth_genome_coverage<
    R: NamedBamReader,
    G: NamedBamReaderGenerator<R>,
    T: CoverageTaker,
>(
    bam_readers: Vec<G>,
    split_char: u8,
    coverage_taker: &mut T,
    print_zero_coverage_genomes: bool,
    coverage_estimators: &mut Vec<CoverageEstimator>,
    flag_filters: &FlagFilter,
    single_genome: bool,
    threads: u16,
) -> Vec<ReadsMapped> {
    let mut reads_mapped_vector = vec![];
    debug!(
        "Calculating coverage using a split character {}",
        str::from_utf8(&[split_char]).unwrap()
    );
    for bam_generator in bam_readers {
        let mut bam_generated = bam_generator.start();
        bam_generated.set_threads(threads as usize);

        let stoit_name = &(bam_generated.name().to_string());
        debug!("Working on stoit {}", stoit_name);
        coverage_taker.start_stoit(stoit_name);
        let header = bam_generated.header().clone();
        let target_names = header.target_names();

        let fill_genome_length_forwards = |current_tid, target_genome: Option<&[u8]>| -> Vec<u64> {
            // Iterating reads skips over contigs with no mapped reads, but the
            // length of these contigs is required to calculate the average
            // across all contigs. This closure returns the number of bases in
            // contigs with tid > current_tid that are part of the current
            // genome.
            if target_genome.is_none() {
                return vec![];
            }
            let mut extras: Vec<u64> = vec![];
            let total_refs = header.target_count();
            let mut my_tid = current_tid + 1;
            while my_tid < total_refs {
                if single_genome
                    || extract_genome(my_tid, &target_names, split_char) == target_genome.unwrap()
                {
                    extras.push(
                        header
                            .target_len(my_tid)
                            .expect("Malformed bam header or programming error encountered"),
                    );
                    my_tid += 1;
                } else {
                    break;
                }
            }
            extras
        };

        let fill_genome_length_backwards_to_last =
            |current_tid, last_tid, target_genome| -> Vec<u64> {
                if current_tid == 0 {
                    return vec![];
                };
                let mut extras: Vec<u64> = vec![];
                let mut my_tid = last_tid + 1;
                while my_tid < current_tid {
                    if single_genome
                        || extract_genome(my_tid, &target_names, split_char) == target_genome
                    {
                        extras.push(
                            header
                                .target_len(my_tid)
                                .expect("Malformed bam header or programming error encountered"),
                        );
                        my_tid += 1;
                    } else {
                        break;
                    }
                }
                extras
            };

        let mut last_tid: u32 = 0;
        let mut doing_first = true;
        let mut last_genome: Option<&[u8]> = None;
        let mut unobserved_contig_length_and_first_tid = UnobservedLengthAndFirstTid {
            unobserved_contig_lengths: vec![],
            first_tid: 0,
        };
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let mut record: bam::record::Record = bam::record::Record::new();
        let mut num_mapped_reads_total: u64 = 0;
        let mut num_mapped_reads_in_current_contig: u64 = 0;
        let mut num_mapped_reads_in_current_genome: u64 = 0;
        let mut total_edit_distance_in_current_contig: u64 = 0;
        let mut total_indels_in_current_contig: u64 = 0;
        loop {
            match bam_generated.read(&mut record) {
                None => {
                    break;
                }
                Some(Ok(())) => {}
                Some(e) => {
                    panic!("Error reading BAM record: {:?}", e)
                }
            }

            if !flag_filters.passes(&record) {
                trace!("Skipping read based on flag filtering");
                continue;
            }
            let original_tid = record.tid();
            if !record.is_unmapped() {
                // if reference has changed, finish a genome or not
                let tid = original_tid as u32;
                let current_genome: &[u8] = match single_genome {
                    true => "".as_bytes(),
                    false => extract_genome(tid, &target_names, split_char),
                };
                if tid != last_tid || doing_first {
                    debug!(
                        "Processing a change in tid, from {} to {} (first is {}). Current \
                           unobserved_and_first: unobserved {:?}, first {}",
                        last_tid,
                        tid,
                        doing_first,
                        unobserved_contig_length_and_first_tid.unobserved_contig_lengths,
                        unobserved_contig_length_and_first_tid.first_tid
                    );
                    if !doing_first && tid < last_tid {
                        error!("BAM file appears to be unsorted. Input BAM files must be sorted by reference (i.e. by samtools sort)");
                        panic!("BAM file appears to be unsorted. Input BAM files must be sorted by reference (i.e. by samtools sort)");
                    }
                    if doing_first {
                        for ref mut coverage_estimator in coverage_estimators.iter_mut() {
                            coverage_estimator.setup()
                        }
                        unobserved_contig_length_and_first_tid = fill_genome_length_backwards(
                            tid,
                            current_genome,
                            single_genome,
                            &target_names,
                            split_char,
                            &header,
                        );
                        last_genome = Some(current_genome);
                        debug!("doing first..");
                        doing_first = false;

                        if print_zero_coverage_genomes && !single_genome {
                            print_previous_zero_coverage_genomes2(
                                None,
                                current_genome,
                                tid,
                                coverage_estimators,
                                &target_names,
                                split_char,
                                coverage_taker,
                                &header,
                            );
                        }
                    } else if current_genome == last_genome.unwrap() {
                        debug!(
                            "Found {} reads mapped to tid {}",
                            num_mapped_reads_in_current_contig, last_tid
                        );
                        for ref mut coverage_estimator in coverage_estimators.iter_mut() {
                            coverage_estimator.add_contig(
                                &ups_and_downs,
                                num_mapped_reads_in_current_contig,
                                total_edit_distance_in_current_contig
                                    - total_indels_in_current_contig,
                            );
                        }
                        // Collect the length of reference sequences from this
                        // genome that had no hits that were just skipped over.
                        debug!("Filling unobserved from {} to {}", last_tid, tid);
                        unobserved_contig_length_and_first_tid
                            .unobserved_contig_lengths
                            .append(&mut fill_genome_length_backwards_to_last(
                                tid,
                                last_tid,
                                current_genome,
                            ));
                    } else {
                        debug!(
                            "Found {} reads mapped to tid {}",
                            num_mapped_reads_in_current_contig, last_tid
                        );
                        // Collect the length of refs from the end of the last genome that had no hits
                        debug!(
                            "Filling unobserved from {} to {} for {}",
                            last_tid,
                            tid,
                            &str::from_utf8(last_genome.unwrap()).unwrap()
                        );
                        unobserved_contig_length_and_first_tid
                            .unobserved_contig_lengths
                            .append(&mut fill_genome_length_backwards_to_last(
                                tid,
                                last_tid,
                                last_genome.unwrap(),
                            ));

                        let positive_coverage = print_last_genomes(
                            num_mapped_reads_in_current_contig,
                            last_genome,
                            &mut unobserved_contig_length_and_first_tid,
                            &ups_and_downs,
                            total_edit_distance_in_current_contig,
                            total_indels_in_current_contig,
                            current_genome,
                            coverage_estimators,
                            coverage_taker,
                            print_zero_coverage_genomes,
                            single_genome,
                            &target_names,
                            split_char,
                            &header,
                            tid,
                        );
                        if positive_coverage {
                            num_mapped_reads_total += num_mapped_reads_in_current_genome;
                        }
                        num_mapped_reads_in_current_genome = 0;
                        last_genome = Some(current_genome);

                        unobserved_contig_length_and_first_tid = fill_genome_length_backwards(
                            tid,
                            current_genome,
                            single_genome,
                            &target_names,
                            split_char,
                            &header,
                        );
                        debug!(
                            "Setting unobserved contig length to be {:?} for genome {}",
                            unobserved_contig_length_and_first_tid.unobserved_contig_lengths,
                            str::from_utf8(current_genome).unwrap()
                        );
                    }

                    ups_and_downs =
                        vec![0; header.target_len(tid).expect("Corrupt BAM file?") as usize];
                    num_mapped_reads_in_current_contig = 0;
                    total_edit_distance_in_current_contig = 0;
                    total_indels_in_current_contig = 0;
                    last_tid = tid;
                }

                // Add coverage info for the current record
                // for each chunk of the cigar string
                trace!(
                    "read name {:?}",
                    std::str::from_utf8(record.qname()).unwrap()
                );
                if !record.is_supplementary() {
                    // Supplementary reads are marked primary, so exclude
                    // supplementary mappings to avoid double counting.
                    num_mapped_reads_in_current_contig += 1;
                    num_mapped_reads_in_current_genome += 1;
                }
                let mut cursor: usize = record.pos() as usize;
                for cig in record.cigar().iter() {
                    trace!("Found cigar {:} from {}", cig, cursor);
                    match cig {
                        Cigar::Match(_) | Cigar::Diff(_) | Cigar::Equal(_) => {
                            // if M, X, or =, increment start and decrement end index
                            trace!(
                                "Adding M, X, or =, at {} and {}",
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
                total_edit_distance_in_current_contig += nm(&record);
            }
        }

        if doing_first && bam_generated.num_detected_primary_alignments() == 0 {
            warn!(
                "No primary alignments were observed for sample {} \
                   - perhaps something went wrong in the mapping?",
                stoit_name
            );
        } else {
            // Print the last genome
            // Give the single genome a dummy name
            if single_genome {
                last_genome = Some("genome1".as_bytes())
            }

            debug!(
                "Found {} reads mapped to tid {}",
                num_mapped_reads_in_current_contig, last_tid
            );
            // Collect the length of refs from the end of the last genome that had no hits
            debug!(
                "Filling unobserved from {} to end for {:?}",
                last_tid,
                match last_genome {
                    None => "No previous genome",
                    Some(g) => str::from_utf8(g).unwrap(),
                }
            );
            unobserved_contig_length_and_first_tid
                .unobserved_contig_lengths
                .append(&mut fill_genome_length_forwards(last_tid, last_genome));

            let positive_coverage = print_last_genomes(
                num_mapped_reads_in_current_contig,
                last_genome,
                &mut unobserved_contig_length_and_first_tid,
                &ups_and_downs,
                total_edit_distance_in_current_contig,
                total_indels_in_current_contig,
                b"",
                coverage_estimators,
                coverage_taker,
                print_zero_coverage_genomes,
                single_genome,
                &target_names,
                split_char,
                &header,
                header.target_count() - 1,
            );
            if positive_coverage {
                num_mapped_reads_total += num_mapped_reads_in_current_genome;
            }
        }

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

        bam_generated.finish();
    }
    reads_mapped_vector
}

fn extract_genome<'a>(tid: u32, target_names: &'a [&[u8]], split_char: u8) -> &'a [u8] {
    let target_name = target_names[tid as usize];
    trace!("target name {:?}, separator {:?}", target_name, split_char);
    let offset = find_first(target_name, split_char).unwrap_or_else(|_| panic!("Contig name {} does not contain split symbol, so cannot determine which genome it belongs to",
                 str::from_utf8(target_name).unwrap()));
    &target_name[0..offset]
}

fn fill_genome_length_backwards(
    current_tid: u32,
    target_genome: &[u8],
    single_genome: bool,
    target_names: &[&[u8]],
    split_char: u8,
    header: &bam::HeaderView,
) -> UnobservedLengthAndFirstTid {
    debug!(
        "At start of fill_genome_length_backwards, found current_tid {}, target_genome {}",
        current_tid,
        str::from_utf8(target_genome).unwrap()
    );
    if current_tid == 0 {
        return UnobservedLengthAndFirstTid {
            unobserved_contig_lengths: vec![],
            first_tid: current_tid as usize,
        };
    }

    let mut extras: Vec<u64> = vec![];
    let mut my_tid = current_tid - 1;
    while single_genome || extract_genome(my_tid, target_names, split_char) == target_genome {
        extras.push(
            header
                .target_len(my_tid)
                .expect("Malformed bam header or programming error encountered"),
        );
        if my_tid == 0 {
            return UnobservedLengthAndFirstTid {
                unobserved_contig_lengths: extras,
                first_tid: 0_usize,
            };
        } else {
            my_tid -= 1;
        }
    }
    debug!(
        "Returning UnobservedLengthAndFirstTid length {:?}, first_tid {}",
        extras,
        my_tid + 1
    );
    UnobservedLengthAndFirstTid {
        unobserved_contig_lengths: extras,
        first_tid: (my_tid + 1) as usize,
    }
}

// Print zero coverage for genomes that have no reads mapped. Genomes are
// detected from the header, counting backwards from the current tid until the
// last seen genome is encountered, or we reach the beginning of the tid array.
#[allow(clippy::too_many_arguments)]
fn print_previous_zero_coverage_genomes2<'a, T: CoverageTaker>(
    last_genome: Option<&[u8]>,
    current_genome: &[u8],
    current_tid: u32,
    pileup_coverage_estimators: &'a Vec<CoverageEstimator>,
    target_names: &[&[u8]],
    split_char: u8,
    coverage_taker: &mut T,
    header: &bam::HeaderView,
) -> &'a Vec<CoverageEstimator> {
    let mut my_current_genome = current_genome;
    let mut tid = current_tid;
    let mut genomes_to_print: Vec<&[u8]> = vec![];
    let mut genome_first_tids: Vec<usize> = vec![];
    let mut genomes_unobserved_length: Vec<u64> = vec![];
    let mut unobserved_length = 0;
    // Need to record the first TID from each genome, but we are iterating down.
    // Gah.
    let mut last_first_id = None;
    loop {
        let genome = extract_genome(tid, target_names, split_char);
        debug!(
            "in print_previous_zero_coverage_genomes2: tid {}, genome {}",
            tid,
            str::from_utf8(genome).unwrap()
        );
        if last_genome.is_some() && genome == last_genome.unwrap() {
            break;
        } else if genome != my_current_genome {
            // In-between genome encountered for the first time.
            // Push the last
            if let Some(id) = last_first_id {
                if last_genome.is_none() || genome != last_genome.unwrap() {
                    genome_first_tids.push(id as usize);
                    genomes_to_print.push(my_current_genome);
                    genomes_unobserved_length.push(unobserved_length);
                }
            }
            my_current_genome = genome;
            last_first_id = Some(tid);
            unobserved_length = header.target_len(tid).unwrap();
        } else if genome != current_genome {
            last_first_id = Some(tid);
            unobserved_length += header.target_len(tid).unwrap();
        }
        if tid == 0 {
            break;
        }
        tid -= 1;
    }
    if let Some(id) = last_first_id {
        genome_first_tids.push(id as usize);
        genomes_to_print.push(my_current_genome);
        genomes_unobserved_length.push(unobserved_length);
    }
    debug!(
        "genomes_to_print {:?}, genome_first_tids {:?}: unobserved: {:?}",
        genomes_to_print, genome_first_tids, genomes_unobserved_length
    );

    for i in (0..genomes_to_print.len()).rev() {
        coverage_taker.start_entry(
            genome_first_tids[i],
            str::from_utf8(genomes_to_print[i]).unwrap(),
        );
        for coverage_estimator in pileup_coverage_estimators {
            coverage_estimator.print_zero_coverage(coverage_taker, genomes_unobserved_length[i]);
        }
        coverage_taker.finish_entry();
    }
    pileup_coverage_estimators
}

#[cfg(test)]
mod tests {
    use super::*;
    use genome_exclusion::*;
    use rust_htslib::bam::Read;
    use shard_bam_reader::*;
    use std::collections::HashSet;
    use std::io::Read as _;
    use OutputWriter;

    fn test_streaming_with_stream<R: NamedBamReader, G: NamedBamReaderGenerator<R>>(
        expected: &str,
        bam_readers: Vec<G>,
        separator: u8,
        print_zero_coverage_contigs: bool,
        coverage_estimators: &mut Vec<CoverageEstimator>,
        proper_pairs_only: bool,
        single_genome: bool,
    ) -> Vec<ReadsMapped> {
        let res;
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        {
            let mut coverage_taker =
                CoverageTakerType::new_single_float_coverage_streaming_coverage_printer(
                    OutputWriter::generate(Some(t)),
                );
            let flags = FlagFilter {
                include_improper_pairs: !proper_pairs_only,
                include_secondary: true,
                include_supplementary: true,
            };
            res = mosdepth_genome_coverage(
                bam_readers,
                separator,
                &mut coverage_taker,
                print_zero_coverage_contigs,
                coverage_estimators,
                &flags,
                single_genome,
                1,
            );
        }
        let mut buf = vec![];
        std::fs::File::open(tf.path())
            .unwrap()
            .read_to_end(&mut buf)
            .unwrap();
        assert_eq!(expected, str::from_utf8(&buf).unwrap());
        res
    }

    fn test_streaming_with_stream_pileup_counts<
        R: NamedBamReader,
        G: NamedBamReaderGenerator<R>,
    >(
        expected: &str,
        bam_readers: Vec<G>,
        separator: u8,
        print_zero_coverage_contigs: bool,
        coverage_estimators: &mut Vec<CoverageEstimator>,
        proper_pairs_only: bool,
        single_genome: bool,
    ) -> Vec<ReadsMapped> {
        let res;
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        {
            let mut coverage_taker = CoverageTakerType::new_pileup_coverage_coverage_printer(
                OutputWriter::generate(Some(t)),
            );
            let flags = FlagFilter {
                include_improper_pairs: !proper_pairs_only,
                include_secondary: false,
                include_supplementary: false,
            };
            res = mosdepth_genome_coverage(
                bam_readers,
                separator,
                &mut coverage_taker,
                print_zero_coverage_contigs,
                coverage_estimators,
                &flags,
                single_genome,
                1,
            );
        }
        let mut buf = vec![];
        std::fs::File::open(tf.path())
            .unwrap()
            .read_to_end(&mut buf)
            .unwrap();
        assert_eq!(expected, str::from_utf8(&buf).unwrap());
        res
    }

    fn test_contig_names_with_stream<R: NamedBamReader, G: NamedBamReaderGenerator<R>>(
        expected: &str,
        bam_readers: Vec<G>,
        geco: &GenomesAndContigs,
        print_zero_coverage_contigs: bool,
        proper_pairs_only: bool,
        coverage_estimators: &mut [CoverageEstimator],
    ) -> Vec<ReadsMapped> {
        let res;
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        {
            let mut coverage_taker =
                CoverageTakerType::new_single_float_coverage_streaming_coverage_printer(
                    OutputWriter::generate(Some(t)),
                );
            let flags = FlagFilter {
                include_improper_pairs: !proper_pairs_only,
                include_secondary: false,
                include_supplementary: false,
            };
            res = mosdepth_genome_coverage_with_contig_names(
                bam_readers,
                geco,
                &mut coverage_taker,
                print_zero_coverage_contigs,
                &flags,
                coverage_estimators,
                1,
            );
        }
        let mut buf = vec![];
        std::fs::File::open(tf.path())
            .unwrap()
            .read_to_end(&mut buf)
            .unwrap();
        assert_eq!(expected, str::from_utf8(&buf).unwrap());
        res
    }

    fn test_contig_names_with_stream_pileup_counts<
        R: NamedBamReader,
        G: NamedBamReaderGenerator<R>,
    >(
        expected: &str,
        bam_readers: Vec<G>,
        geco: &GenomesAndContigs,
        print_zero_coverage_contigs: bool,
        proper_pairs_only: bool,
        coverage_estimators: &mut [CoverageEstimator],
    ) -> Vec<ReadsMapped> {
        let res;
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        {
            let mut coverage_taker = CoverageTakerType::new_pileup_coverage_coverage_printer(
                OutputWriter::generate(Some(t)),
            );
            let flags = FlagFilter {
                include_improper_pairs: !proper_pairs_only,
                include_secondary: false,
                include_supplementary: false,
            };
            res = mosdepth_genome_coverage_with_contig_names(
                bam_readers,
                geco,
                &mut coverage_taker,
                print_zero_coverage_contigs,
                &flags,
                coverage_estimators,
                1,
            );
        }
        let mut buf = vec![];
        std::fs::File::open(tf.path())
            .unwrap()
            .read_to_end(&mut buf)
            .unwrap();
        assert_eq!(expected, str::from_utf8(&buf).unwrap());
        res
    }

    #[test]
    fn test_one_genome_two_contigs_first_covered() {
        test_streaming_with_stream(
            "2seqs.reads_for_seq1\tse\t0.6\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            b'q',
            true,
            &mut vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
            false,
            false,
        );
    }

    #[test]
    fn test_one_genome_two_contigs_first_covered_contig_names() {
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("se".to_string());
        geco.insert("seq1".to_string(), genome1);
        geco.insert("seq2".to_string(), genome1);
        test_contig_names_with_stream(
            "2seqs.reads_for_seq1\tse\t0.6\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            &geco,
            true,
            false,
            &mut [CoverageEstimator::new_estimator_mean(0.0, 0, false)],
        );
    }

    #[test]
    fn test_one_genome_two_contigs_second_covered() {
        test_streaming_with_stream(
            "2seqs.reads_for_seq2\tse\t0.6\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq2.bam"]),
            b'q',
            true,
            &mut vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
            false,
            false,
        );
    }

    #[test]
    fn test_one_genome_two_contigs_second_covered_contig_names() {
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("se".to_string());
        geco.insert("seq1".to_string(), genome1);
        geco.insert("seq2".to_string(), genome1);
        test_contig_names_with_stream(
            "2seqs.reads_for_seq2\tse\t0.6\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq2.bam"]),
            &geco,
            true,
            false,
            &mut [CoverageEstimator::new_estimator_mean(0.0, 0, false)],
        );
    }

    #[test]
    fn test_one_genome_two_contigs_both_covered() {
        test_streaming_with_stream(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/2seqs.reads_for_seq1_and_seq2.bam",
            ]),
            b'e',
            true,
            &mut vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
            false,
            false,
        );
    }

    #[test]
    fn test_one_genome_two_contigs_both_covered_contig_names() {
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(), genome1);
        geco.insert("seq2".to_string(), genome1);
        test_contig_names_with_stream(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/2seqs.reads_for_seq1_and_seq2.bam",
            ]),
            &geco,
            true,
            false,
            &mut [CoverageEstimator::new_estimator_mean(0.0, 0, false)],
        );
    }

    #[test]
    fn test_one_genome_min_fraction_covered_under_min() {
        test_streaming_with_stream(
            "2seqs.reads_for_seq1_and_seq2\ts\t0\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/2seqs.reads_for_seq1_and_seq2.bam",
            ]),
            b'e',
            true,
            &mut vec![CoverageEstimator::new_estimator_mean(0.76, 0, false)],
            false,
            false,
        );
    }

    #[test]
    fn test_one_genome_min_fraction_covered_under_min_contig_names() {
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(), genome1);
        geco.insert("seq2".to_string(), genome1);
        test_contig_names_with_stream(
            "",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/2seqs.reads_for_seq1_and_seq2.bam",
            ]),
            &geco,
            false,
            false,
            &mut [CoverageEstimator::new_estimator_mean(0.76, 0, false)],
        );
    }

    #[test]
    fn test_one_genome_min_fraction_covered_just_ok() {
        test_streaming_with_stream(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/2seqs.reads_for_seq1_and_seq2.bam",
            ]),
            b'e',
            true,
            &mut vec![CoverageEstimator::new_estimator_mean(0.759, 0, false)],
            false,
            false,
        );
    }

    #[test]
    fn test_one_genome_min_fraction_covered_just_ok_contig_names() {
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(), genome1);
        geco.insert("seq2".to_string(), genome1);
        test_contig_names_with_stream(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/2seqs.reads_for_seq1_and_seq2.bam",
            ]),
            &geco,
            true,
            false,
            &mut [CoverageEstimator::new_estimator_mean(0.759, 0, false)],
        );
    }

    #[test]
    fn test_two_contigs_trimmed_mean() {
        test_streaming_with_stream(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.08875\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/2seqs.reads_for_seq1_and_seq2.bam",
            ]),
            b'e',
            true,
            &mut vec![CoverageEstimator::new_estimator_trimmed_mean(
                0.1, 0.9, 0.759, 0,
            )],
            false,
            false,
        );
    }

    #[test]
    fn test_two_contigs_trimmed_mean_contig_names() {
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(), genome1);
        geco.insert("seq2".to_string(), genome1);
        test_contig_names_with_stream(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.08875\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/2seqs.reads_for_seq1_and_seq2.bam",
            ]),
            &geco,
            true,
            false,
            &mut [CoverageEstimator::new_estimator_trimmed_mean(
                0.1, 0.9, 0.0, 0,
            )],
        );
    }

    #[test]
    fn test_two_contigs_pileup_counts_estimator() {
        test_streaming_with_stream_pileup_counts(
            "2seqs.reads_for_seq1_and_seq2\ts\t0\t482\n2seqs.reads_for_seq1_and_seq2\ts\t1\t922\n2seqs.reads_for_seq1_and_seq2\ts\t2\t371\n2seqs.reads_for_seq1_and_seq2\ts\t3\t164\n2seqs.reads_for_seq1_and_seq2\ts\t4\t61\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            b'e',
            true,
            &mut vec!(CoverageEstimator::new_estimator_pileup_counts(0.0,0)),
            false,
            false);
    }

    #[test]
    fn test_two_contigs_pileup_counts_estimator_contig_names() {
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("s".to_string());
        geco.insert("seq1".to_string(), genome1);
        geco.insert("seq2".to_string(), genome1);
        test_contig_names_with_stream_pileup_counts(
            "2seqs.reads_for_seq1_and_seq2\ts\t0\t482\n2seqs.reads_for_seq1_and_seq2\ts\t1\t922\n2seqs.reads_for_seq1_and_seq2\ts\t2\t371\n2seqs.reads_for_seq1_and_seq2\ts\t3\t164\n2seqs.reads_for_seq1_and_seq2\ts\t4\t61\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            true,
            false,
            &mut [CoverageEstimator::new_estimator_pileup_counts(0.0,0)]);
    }

    #[test]
    fn test_zero_coverage_genomes() {
        test_streaming_with_stream(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0\n7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome3\t0\n7seqs.reads_for_seq1_and_seq2\tgenome4\t0\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6\t0\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            b'~',
            true,
            &mut vec!(CoverageEstimator::new_estimator_mean(0.1,0,false)),
            false,
            false);

        test_streaming_with_stream(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            b'~',
            false,
            &mut vec!(CoverageEstimator::new_estimator_mean(0.1,0,false)),
            false,
            false);
    }

    #[test]
    fn test_sharded_bams_with_zero_coverage() {
        test_streaming_with_stream(
            "shard1|shard2\tgenome3\t0.10908099\nshard1|shard2\tgenome4\t0.109071076\nshard1|shard2\tgenome5\t0\nshard1|shard2\tgenome6\t0.10906117\nshard1|shard2\tgenome1\t0.10904135\nshard1|shard2\tgenome2\t0\n",
            generate_sharded_bam_reader_from_bam_files(
                vec!["tests/data/shard1.bam", "tests/data/shard2.bam"],
                4,
                &NoExclusionGenomeFilter{}),
            b'~',
            true,
            &mut vec!(CoverageEstimator::new_estimator_mean(0.1,0,false)),
            false,
            false);
    }

    #[test]
    fn test_sharded_bams_with_genome_exclusion() {
        let mut hashset: HashSet<&[u8]> = HashSet::new();
        hashset.insert(b"genome3");
        let ex = SeparatorGenomeExclusionFilter {
            split_char: b'~',
            excluded_genomes: hashset,
        };
        test_streaming_with_stream(
            "shard1|shard2\tgenome3\t0\n\
             shard1|shard2\tgenome4\t0.109071076\n\
             shard1|shard2\tgenome5\t0\n\
             shard1|shard2\tgenome6\t0.10906117\n\
             shard1|shard2\tgenome1\t0.10904135\n\
             shard1|shard2\tgenome2\t0\n",
            generate_sharded_bam_reader_from_bam_files(
                vec!["tests/data/shard1.bam", "tests/data/shard2.bam"],
                4,
                &ex,
            ),
            b'~',
            true,
            &mut vec![CoverageEstimator::new_estimator_mean(0.1, 0, false)],
            false,
            false,
        );
    }

    #[test]
    fn test_zero_coverage_genomes_after_min_fraction() {
        test_streaming_with_stream(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0\n7seqs.reads_for_seq1_and_seq2\tgenome2\t0\n7seqs.reads_for_seq1_and_seq2\tgenome3\t0\n7seqs.reads_for_seq1_and_seq2\tgenome4\t0\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6\t0\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            b'~',
            true,
            &mut vec!(CoverageEstimator::new_estimator_mean(0.759,0,false)),
            false,
            false);
    }

    #[test]
    fn test_single_genome() {
        test_streaming_with_stream(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0.04209345\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
            ]),
            b'~',
            true,
            &mut vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
            false,
            true,
        );
    }

    #[test]
    fn test_covered_bases_estimator() {
        // $ samtools depth 7seqs.reads_for_seq1_and_seq2.bam  |cut -f1 |us
        // 849 genome5~seq2
        // 669 genome2~seq1
        test_streaming_with_stream(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t669\n\
             7seqs.reads_for_seq1_and_seq2\tgenome5\t849\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
            ]),
            b'~',
            false,
            &mut vec![CoverageEstimator::new_estimator_covered_bases(0.0)],
            false,
            false,
        );
    }

    #[test]
    fn test_covered_bases_estimator_contig_end_exclusion() {
        // $ samtools depth 7seqs.reads_for_seq1_and_seq2.bam |awk '$2 > 75 && $2 < 926' |cut -f1 |us
        // 714 genome5~seq2
        // 669 genome2~seq1
        test_streaming_with_stream(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t669\n\
             7seqs.reads_for_seq1_and_seq2\tgenome5\t849\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
            ]),
            b'~',
            false,
            &mut vec![CoverageEstimator::new_estimator_covered_bases(0.0)],
            false,
            false,
        );
    }

    #[test]
    fn test_zero_coverage_genomes_contig_names() {
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
        geco.insert("genome1~random_sequence_length_11000".to_string(), genome1);
        geco.insert("genome1~random_sequence_length_11010".to_string(), genome1);
        geco.insert("genome2~seq1".to_string(), genome2);
        geco.insert("genome3~random_sequence_length_11001".to_string(), genome3);
        geco.insert("genome4~random_sequence_length_11002".to_string(), genome4);
        geco.insert("genome5~seq2".to_string(), genome5);
        geco.insert("genome6~random_sequence_length_11003".to_string(), genome6);
        test_contig_names_with_stream(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0\n\
            7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n\
            7seqs.reads_for_seq1_and_seq2\tgenome3\t0\n\
            7seqs.reads_for_seq1_and_seq2\tgenome4\t0\n\
            7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n\
            7seqs.reads_for_seq1_and_seq2\tgenome6\t0\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
            ]),
            &geco,
            true,
            false,
            &mut [CoverageEstimator::new_estimator_mean(0.1, 0, false)],
        );

        test_contig_names_with_stream(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            false,
            false,
            &mut [CoverageEstimator::new_estimator_mean(0.1,0,false)]);
    }

    #[test]
    fn test_zero_coverage_genomes_contig_names_with_multiple_methods() {
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
        geco.insert("genome1~random_sequence_length_11000".to_string(), genome1);
        geco.insert("genome1~random_sequence_length_11010".to_string(), genome1);
        geco.insert("genome2~seq1".to_string(), genome2);
        geco.insert("genome3~random_sequence_length_11001".to_string(), genome3);
        geco.insert("genome4~random_sequence_length_11002".to_string(), genome4);
        geco.insert("genome5~seq2".to_string(), genome5);
        geco.insert("genome6~random_sequence_length_11003".to_string(), genome6);
        test_contig_names_with_stream(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0\t0\n\
            7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\t1.3633634\n\
            7seqs.reads_for_seq1_and_seq2\tgenome3\t0\t0\n\
            7seqs.reads_for_seq1_and_seq2\tgenome4\t0\t0\n\
            7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\t0.6166166\n\
            7seqs.reads_for_seq1_and_seq2\tgenome6\t0\t0\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
            ]),
            &geco,
            true,
            false,
            &mut [
                CoverageEstimator::new_estimator_mean(0.1, 0, false),
                CoverageEstimator::new_estimator_variance(0.1, 0),
            ],
        );

        let reads_mapped = test_contig_names_with_stream(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\t1.3633634\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\t0.6166166\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/7seqs.reads_for_seq1_and_seq2.bam"]),
            &geco,
            false,
            false,
            &mut [CoverageEstimator::new_estimator_mean(0.1,0,false),
                CoverageEstimator::new_estimator_variance(0.1,0)]);
        assert_eq!(
            vec!(ReadsMapped {
                num_mapped_reads: 24,
                num_reads: 24
            }),
            reads_mapped
        );
    }

    #[test]
    fn test_genomes_and_contigs_reads_mapped() {
        let mut geco = GenomesAndContigs::new();
        // >genome1~random_sequence_length_11000
        //     >genome1~random_sequence_length_11010
        //     >genome2~seq1
        //     >genome3~random_sequence_length_11001
        //     >genome4~random_sequence_length_11002
        //     >genome5~seq2
        //     >genome6~random_sequence_length_11003
        let genome2 = geco.establish_genome("genome2".to_string());
        let genome3 = geco.establish_genome("genome3".to_string());
        geco.insert("genome2~seq1".to_string(), genome2);
        geco.insert("genome3~random_sequence_length_11001".to_string(), genome3);
        let reads_mapped = test_contig_names_with_stream(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\t1.3633634\n\
            7seqs.reads_for_seq1_and_seq2\tgenome3\t0\t0\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
            ]),
            &geco,
            true,
            false,
            &mut [
                CoverageEstimator::new_estimator_mean(0.1, 0, false),
                CoverageEstimator::new_estimator_variance(0.1, 0),
            ],
        );
        assert_eq!(
            vec!(ReadsMapped {
                num_mapped_reads: 12,
                num_reads: 24
            }),
            reads_mapped
        );
    }

    #[test]
    fn test_julian_error() {
        let reads_mapped = test_streaming_with_stream(
            "2seqs.reads_for_seq1.with_unmapped\tgenome1\t1.4985\n",
            // has unmapped reads, which caused problems with --proper-pairs-only
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/2seqs.reads_for_seq1.with_unmapped.bam",
            ]),
            b'\0',
            true,
            &mut vec![CoverageEstimator::new_estimator_mean(0.1, 0, true)],
            false,
            true,
        );
        assert_eq!(
            vec!(ReadsMapped {
                num_mapped_reads: 20,
                num_reads: 24
            }),
            reads_mapped
        );
    }

    #[test]
    fn test_multiple_outputs_one_zero_no_print_zeroes_single_genome() {
        test_streaming_with_stream(
            "2seqs.reads_for_seq1\tgenome1\t0.6\t0\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            b'q',
            true,
            &mut vec![
                CoverageEstimator::new_estimator_mean(0.0, 0, false),
                // covered fraction is 0.727, so go lower so trimmed mean is 0,
                // mean > 0.
                CoverageEstimator::new_estimator_trimmed_mean(0.0, 0.05, 0.0, 0),
            ],
            false,
            true,
        );
    }

    #[test]
    fn test_multiple_outputs_one_zero_no_print_zeroes_single_genome_reverse() {
        test_streaming_with_stream(
            "2seqs.reads_for_seq1\tgenome1\t0\t0.6\n",
            generate_named_bam_readers_from_bam_files(vec!["tests/data/2seqs.reads_for_seq1.bam"]),
            b'q',
            true,
            &mut vec![
                CoverageEstimator::new_estimator_trimmed_mean(0.0, 0.05, 0.0, 0),
                // covered fraction is 0.727, so go lower so trimmed mean is 0,
                // mean > 0.
                CoverageEstimator::new_estimator_mean(0.0, 0, false),
            ],
            false,
            true,
        );
    }

    #[test]
    fn test_multiple_outputs_one_zero_no_print_zeroes_separator() {
        test_streaming_with_stream(
            "7seqs.reads_for_seq1\tgenome1\t0\t0\n7seqs.reads_for_seq1\tgenome2\t1.2\t0\n7seqs.reads_for_seq1\tgenome3\t0\t0\n7seqs.reads_for_seq1\tgenome4\t0\t0\n7seqs.reads_for_seq1\tgenome5\t0\t0\n7seqs.reads_for_seq1\tgenome6\t0\t0\n",
            generate_named_bam_readers_from_bam_files(
                vec!["tests/data/7seqs.reads_for_seq1.bam"]),
            b'~',
            true,
            &mut vec!(
                CoverageEstimator::new_estimator_mean(0.0,0,false),
                // covered fraction is 0.727, so go lower so trimmed mean is 0,
                // mean > 0.
                CoverageEstimator::new_estimator_trimmed_mean(0.0,0.05,0.0,0)
            ),
            false,
            false);
    }

    #[test]
    fn test_multiple_outputs_one_zero_no_print_zeroes_separator_reverse() {
        test_streaming_with_stream(
            "7seqs.reads_for_seq1\tgenome1\t0\t0\n7seqs.reads_for_seq1\tgenome2\t0\t1.2\n7seqs.reads_for_seq1\tgenome3\t0\t0\n7seqs.reads_for_seq1\tgenome4\t0\t0\n7seqs.reads_for_seq1\tgenome5\t0\t0\n7seqs.reads_for_seq1\tgenome6\t0\t0\n",
            generate_named_bam_readers_from_bam_files(
                vec!["tests/data/7seqs.reads_for_seq1.bam"]),
            b'~',
            true,
            &mut vec!(
                // covered fraction is 0.727, so go lower so trimmed mean is 0,
                // mean > 0.
                CoverageEstimator::new_estimator_trimmed_mean(0.0,0.05,0.0,0),
                CoverageEstimator::new_estimator_mean(0.0,0,false),
            ),
            false,
            false);
    }

    #[test]
    fn test_fill_genome_length_backwards1() {
        let bam = bam::Reader::from_path("tests/data/7seqs.reads_for_seq1.bam").unwrap();
        let uf = fill_genome_length_backwards(
            2,
            "genome2".as_bytes(),
            false,
            &bam.header().target_names(),
            b'~',
            bam.header(),
        );
        assert_eq!(vec![0; 0], uf.unobserved_contig_lengths);
        assert_eq!(2, uf.first_tid);
    }

    #[test]
    fn test_fill_genome_length_backwards2() {
        let bam = bam::Reader::from_path("tests/data/7seqs.reads_for_seq1.bam").unwrap();
        let uf = fill_genome_length_backwards(
            1,
            "genome1".as_bytes(),
            false,
            &bam.header().target_names(),
            b'~',
            bam.header(),
        );
        assert_eq!(vec![11000], uf.unobserved_contig_lengths);
        assert_eq!(0, uf.first_tid);
    }

    #[test]
    fn test_fill_genome_length_backwards3() {
        let bam = bam::Reader::from_path("tests/data/7seqs.reads_for_seq1.bam").unwrap();
        let uf = fill_genome_length_backwards(
            5,
            "genome5".as_bytes(),
            false,
            &bam.header().target_names(),
            b'~',
            bam.header(),
        );
        assert_eq!(
            UnobservedLengthAndFirstTid {
                unobserved_contig_lengths: vec![],
                first_tid: 5
            },
            uf
        );
        assert_eq!(5, uf.first_tid);
    }

    #[test]
    fn test_print_previous_zero_coverage_genomes2() {
        let mut coverage_taker = CoverageTakerType::new_cached_single_float_coverage_taker(1);
        let estimator = CoverageEstimator::new_estimator_mean(0.0, 0, false);
        coverage_taker.start_stoit("dummy");
        let bam = bam::Reader::from_path("tests/data/7seqs.reads_for_seq1.bam").unwrap();
        print_previous_zero_coverage_genomes2(
            None,
            "genome2".as_bytes(),
            2,
            &vec![estimator],
            &bam.header().target_names(),
            b'~',
            &mut coverage_taker,
            bam.header(),
        );
        match coverage_taker {
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                stoit_names,
                entry_names,
                ..
            } => {
                assert_eq!(vec!("dummy"), stoit_names);
                assert_eq!(vec!(Some("genome1".to_string())), entry_names);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_print_previous_zero_coverage_genomes2_length_estimator() {
        let mut coverage_taker = CoverageTakerType::new_cached_single_float_coverage_taker(1);
        let estimator = CoverageEstimator::new_estimator_length();
        coverage_taker.start_stoit("dummy");
        let bam = bam::Reader::from_path("tests/data/7seqs.reads_for_seq1.bam").unwrap();
        print_previous_zero_coverage_genomes2(
            None,
            "genome2".as_bytes(),
            2,
            &vec![estimator],
            &bam.header().target_names(),
            b'~',
            &mut coverage_taker,
            bam.header(),
        );
        match coverage_taker {
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                stoit_names,
                entry_names,
                coverages,
                ..
            } => {
                assert_eq!(vec!("dummy"), stoit_names);
                assert_eq!(vec!(Some("genome1".to_string())), entry_names);
                assert_eq!(
                    vec!(vec!(CoverageEntry {
                        entry_index: 0,
                        coverage: 22010.0
                    })),
                    coverages
                )
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_read_count_calculator() {
        let res = test_streaming_with_stream(
            "7seqs.reads_for_seq1\tgenome1\t0\n\
             7seqs.reads_for_seq1\tgenome2\t12\n\
             7seqs.reads_for_seq1\tgenome3\t0\n\
             7seqs.reads_for_seq1\tgenome4\t0\n\
             7seqs.reads_for_seq1\tgenome5\t0\n\
             7seqs.reads_for_seq1\tgenome6\t0\n\
             7seqs.reads_for_seq1_and_seq2\tgenome1\t0\n\
             7seqs.reads_for_seq1_and_seq2\tgenome2\t12\n\
             7seqs.reads_for_seq1_and_seq2\tgenome3\t0\n\
             7seqs.reads_for_seq1_and_seq2\tgenome4\t0\n\
             7seqs.reads_for_seq1_and_seq2\tgenome5\t12\n\
             7seqs.reads_for_seq1_and_seq2\tgenome6\t0\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/7seqs.reads_for_seq1.bam",
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
            ]),
            b'~',
            true,
            &mut vec![
                // covered fraction is 0.727, so go lower so trimmed mean is 0,
                // mean > 0.
                CoverageEstimator::new_estimator_read_count(),
            ],
            false,
            false,
        );
        assert_eq!(
            vec!(
                ReadsMapped {
                    num_mapped_reads: 12,
                    num_reads: 12
                },
                ReadsMapped {
                    num_mapped_reads: 24,
                    num_reads: 24
                }
            ),
            res
        );
    }

    #[test]
    fn test_read_count_below_min_covered() {
        // First test with all genomes passing the covered fraction threshold
        let res = test_streaming_with_stream(
            "7seqs.reads_for_seq1\tgenome1\t0\t0\n\
             7seqs.reads_for_seq1\tgenome2\t12\t0.727\n\
             7seqs.reads_for_seq1\tgenome3\t0\t0\n\
             7seqs.reads_for_seq1\tgenome4\t0\t0\n\
             7seqs.reads_for_seq1\tgenome5\t0\t0\n\
             7seqs.reads_for_seq1\tgenome6\t0\t0\n\
             7seqs.reads_for_seq1_and_seq2\tgenome1\t0\t0\n\
             7seqs.reads_for_seq1_and_seq2\tgenome2\t12\t0.669\n\
             7seqs.reads_for_seq1_and_seq2\tgenome3\t0\t0\n\
             7seqs.reads_for_seq1_and_seq2\tgenome4\t0\t0\n\
             7seqs.reads_for_seq1_and_seq2\tgenome5\t12\t0.849\n\
             7seqs.reads_for_seq1_and_seq2\tgenome6\t0\t0\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/7seqs.reads_for_seq1.bam",
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
            ]),
            b'~',
            true,
            &mut vec![
                // covered fraction is 0.727, so go lower so trimmed mean is 0,
                // mean > 0.
                CoverageEstimator::new_estimator_read_count(),
                CoverageEstimator::new_estimator_covered_fraction(0.1),
            ],
            false,
            false,
        );
        assert_eq!(
            vec!(
                ReadsMapped {
                    num_mapped_reads: 12,
                    num_reads: 12
                },
                ReadsMapped {
                    num_mapped_reads: 24,
                    num_reads: 24
                }
            ),
            res
        );

        // Next test when the covered fraction fails
        let res = test_streaming_with_stream(
            "7seqs.reads_for_seq1\tgenome1\t0\n\
             7seqs.reads_for_seq1\tgenome2\t0\n\
             7seqs.reads_for_seq1\tgenome3\t0\n\
             7seqs.reads_for_seq1\tgenome4\t0\n\
             7seqs.reads_for_seq1\tgenome5\t0\n\
             7seqs.reads_for_seq1\tgenome6\t0\n\
             7seqs.reads_for_seq1_and_seq2\tgenome1\t0\n\
             7seqs.reads_for_seq1_and_seq2\tgenome2\t0\n\
             7seqs.reads_for_seq1_and_seq2\tgenome3\t0\n\
             7seqs.reads_for_seq1_and_seq2\tgenome4\t0\n\
             7seqs.reads_for_seq1_and_seq2\tgenome5\t0\n\
             7seqs.reads_for_seq1_and_seq2\tgenome6\t0\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/7seqs.reads_for_seq1.bam",
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
            ]),
            b'~',
            true,
            &mut vec![
                // All genomes do not pass the min covered fraction threshold.
                CoverageEstimator::new_estimator_covered_fraction(0.99),
            ],
            false,
            false,
        );
        assert_eq!(
            vec!(
                ReadsMapped {
                    num_mapped_reads: 0,
                    num_reads: 12
                },
                ReadsMapped {
                    num_mapped_reads: 0,
                    num_reads: 24
                }
            ),
            res
        );
    }

    #[test]
    fn test_read_count_below_min_covered_genomes_and_contigs() {
        let mut geco = GenomesAndContigs::new();
        let genome1 = geco.establish_genome("genome1".to_string());
        let genome2 = geco.establish_genome("genome2".to_string());
        let genome3 = geco.establish_genome("genome3".to_string());
        let genome4 = geco.establish_genome("genome4".to_string());
        let genome5 = geco.establish_genome("genome5".to_string());
        let genome6 = geco.establish_genome("genome6".to_string());
        geco.insert("genome1~random_sequence_length_11000".to_string(), genome1);
        geco.insert("genome1~random_sequence_length_11010".to_string(), genome1);
        geco.insert("genome2~seq1".to_string(), genome2);
        geco.insert("genome3~random_sequence_length_11001".to_string(), genome3);
        geco.insert("genome4~random_sequence_length_11002".to_string(), genome4);
        geco.insert("genome5~seq2".to_string(), genome5);
        geco.insert("genome6~random_sequence_length_11003".to_string(), genome6);

        // First test with all genomes passing the covered fraction threshold
        let mut reads_mapped = test_contig_names_with_stream(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\t1.3633634\n\
             7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\t0.6166166\n",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
            ]),
            &geco,
            false,
            false,
            &mut [
                CoverageEstimator::new_estimator_mean(0.1, 0, false),
                CoverageEstimator::new_estimator_variance(0.1, 0),
            ],
        );
        assert_eq!(
            vec!(ReadsMapped {
                num_mapped_reads: 24,
                num_reads: 24
            }),
            reads_mapped
        );

        // Then test when the reads do not make the threshold
        reads_mapped = test_contig_names_with_stream(
            "",
            generate_named_bam_readers_from_bam_files(vec![
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
            ]),
            &geco,
            false,
            false,
            &mut [
                CoverageEstimator::new_estimator_mean(0.99, 0, false),
                CoverageEstimator::new_estimator_variance(0.99, 0),
            ],
        );
        assert_eq!(
            vec!(ReadsMapped {
                num_mapped_reads: 0,
                num_reads: 24
            }),
            reads_mapped
        );
    }
}
