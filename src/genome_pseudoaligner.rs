use std::collections::{BTreeMap, BTreeSet};

use kmer_coverage::*;
use genomes_and_contigs::GenomesAndContigs;

use pseudoaligner::*;
use debruijn::Kmer;
use bio::io::fastq;
use log::Level;

pub fn calculate_genome_kmer_coverage<K>(
    forward_fastq: &str,
    num_threads: usize,
    print_zero_coverage_contigs: bool,
    index: &DebruijnIndex<K>,
    genomes_and_contigs: &GenomesAndContigs)
where K: Kmer + Sync + Send {

    // Do the mappings
    let reads = fastq::Reader::from_file(forward_fastq).expect("Failure to read reference sequences");
    let (eq_class_indices, eq_class_coverages, _eq_class_counts) = pseudoaligner::process_reads(
        reads, &index.index, num_threads)
        .expect("Failure during mapping process");
    info!("Finished mapping reads!");

    if log_enabled!(Level::Debug) {
        for (eq_class, index) in &eq_class_indices {
            debug!("eq_class_index: {:?}\t{}", eq_class, index);
        };
        for (i, eq_class_coverage) in eq_class_coverages.iter().enumerate() {
            debug!("eq_class_coverage {} {}", i, eq_class_coverage);
        }
    }

    let kmer_coverage_total: f64 = eq_class_coverages.iter().sum::<usize>() as f64;

    // Make a new vec containing the eq_class memberships
    let mut eq_classes = vec![None; eq_class_coverages.len()];
    for (eq_class, index) in &eq_class_indices {
        eq_classes[*index] = Some(eq_class)
    }

    // Generate some useful data structures that are used repeatedly in the EM
    // process.
    let genome_contigs = generate_genome_to_contig_indices_vec(
        &genomes_and_contigs,
        &index.tx_names
    );
    let contig_idx_to_genome_idx: Vec<usize> = index.tx_names.iter().map(
        |contig_name|
        genomes_and_contigs.contig_to_genome[contig_name]
    ).collect();
    let mut genome_lengths: Vec<usize> = vec![0; genomes_and_contigs.genomes.len()];
    for (genome_idx, contig_indices) in genome_contigs.iter().enumerate() {
        genome_lengths[genome_idx] = contig_indices.iter()
            .map(|i| index.seq_lengths[*i])
            .sum();
    }
    debug!("Genome contigs: {:?}", genome_contigs);
    debug!("contig_idx_to_genome_idx: {:?}", contig_idx_to_genome_idx);
    debug!("genome_lengths: {:?}", genome_lengths);

    // EM algorithm
    // Randomly assign random relative abundances to each of the genomes
    // TODO: remove: for dev, assign each the same abundance
    info!("Starting EM process ..");
    let mut genome_to_relative_abundance = vec![2.0; genomes_and_contigs.genomes.len()];
    let mut genome_to_read_count;
    let mut num_covergence_steps: u32 = 0;

    loop { // loop until converged
        // E-step: Determine the number of reads we expect come from each contig
        // given their relative abundance.

        // for each equivalence class / count pair, we expect for genome A the
        // abundance of A divided by the sum of abundances of genomes in the
        // equivalence class.
        genome_to_read_count = vec![0.0; genomes_and_contigs.genomes.len()];
        for (i, coverage) in eq_class_coverages.iter().enumerate() {
            match eq_classes[i] {
                None => unreachable!(),
                Some(ref eqs) => {
                    // Each genome is only counted once, instead of having it
                    // possible that 2 contigs from the same genome both match.
                    // This makes the analysis consistent between 2 contigs with
                    // a repeat vs 1 contig with 2 copies of the repeat.
                    let mut seen_genomes: BTreeSet<usize> = BTreeSet::new();

                    let relative_abundances: Vec<f64> = eqs.iter().map(|eq| {

                        let genome_index = contig_idx_to_genome_idx[*eq as usize];
                        match seen_genomes.insert(genome_index) {
                            true => {
                                // Genome not previously seen here
                                genome_to_relative_abundance[genome_index]
                            },
                            false => {
                                // Genome already seen, ignore it
                                0.0
                            }
                        }
                    }).collect();

                    // Add coverages divided by the sum of relative abundances
                    let total_abundance_of_matching_contigs: f64 = relative_abundances.iter().sum();
                    for (eq, relabund) in eqs.iter().zip(relative_abundances.iter()) {
                        genome_to_read_count[contig_idx_to_genome_idx[*eq as usize]] += (*coverage as f64) *
                            relabund / total_abundance_of_matching_contigs;
                    }
                }
            }
        }
        debug!("After E-step have genome coverages: {:?}", genome_to_read_count);

        // M-step: Work out the relative abundance given the number of reads
        // predicted in the E-step. Relative abundance is just the ratio of the
        // coverages, weighted by the inverse of each contig's length.
        //
        // Or, for genome mode, weighted by the inverse of the sum of the
        // genome's length.
        //TODO: zip here instead?
        // First determine the total scaled abundance
        //
        // TODO: Maybe doesn't make sense to loop over genomes that have
        // no coverage, can we not do that?
        let genome_coverages: Vec<f64> = genome_to_read_count.iter().zip(&genome_lengths).map(
            |(cov, l)|
            cov / (*l as f64)).collect();
        let total_coverage: f64 = genome_coverages.iter().sum();
        debug!("Found genome coverages {:?} for total {}", genome_coverages, total_coverage);

        // Next set the abundances so the total == 1.0
        //
        // Converge when all contigs with total abundance > 0.01 * kmer_coverage_total
        // change abundance by < 1%.
        let mut converge = true;
        for (i, cov) in genome_coverages.iter().enumerate() {
            let new_relabund = cov / total_coverage;
            if new_relabund * kmer_coverage_total > 0.01*100.0 {// Add 100 in there
                // as a rough mapped kmers per aligning fragment.
                // TODO: Use number of reads aligned instead?
                let delta = new_relabund / genome_to_relative_abundance[i];
                debug!("For testing convergence of index {}, found delta {}", i, delta);
                if delta < 0.99 || delta > 1.01 {
                    debug!("No converge for you");
                    converge = false;
                }
            }
            genome_to_relative_abundance[i] = new_relabund;
        }
        debug!("At end of M-step, have relative abundances: {:?}", genome_to_relative_abundance);

        num_covergence_steps += 1;
        if converge {
            info!("EM process converged after {} steps", num_covergence_steps);
            break;
        }
    }

    // Print results
    println!("contig\tkmer"); //TODO: Use the same methods as elsewhere for printing.
    for (i, total_coverage) in genome_to_read_count.iter().enumerate() {
        if print_zero_coverage_contigs || *total_coverage > 0.0 {
            // Print average coverage as total coverage divided by contig length.
            println!("{}\t{}", genomes_and_contigs.genomes[i], total_coverage / (index.seq_lengths[i] as f64));
        }
    }
    info!("Finished printing genome coverages");
}

fn generate_genome_to_contig_indices_vec(
    genomes_and_contigs: &GenomesAndContigs,
    tx_names: &Vec<String>
) -> Vec<Vec<usize>> {
    let mut tx_ids_of_genomes: Vec<Vec<usize>> = vec!(
        vec!(); genomes_and_contigs.genomes.len());

    let mut contig_to_index: BTreeMap<&str, usize> = BTreeMap::new();
    for (i, tx_name) in tx_names.iter().enumerate() {
        contig_to_index.insert(tx_name, i);
    }
    for (contig, genome_id) in genomes_and_contigs.contig_to_genome.iter() {
        let contig_str: &str = &contig;
        tx_ids_of_genomes[*genome_id].push(contig_to_index[contig_str]);
    }

    return tx_ids_of_genomes
}
