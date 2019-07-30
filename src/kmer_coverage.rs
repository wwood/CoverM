use pseudoaligner::build_index::build_index;
use pseudoaligner::*;

use bio::io::{fasta, fastq};

pub fn calculate_genome_kmer_coverage(
    reference_path: &str,
    forward_fastq: &str,
    print_zero_coverage_contigs: bool) {

    // Build index
    let reference_reader = fasta::Reader::from_file(reference_path).expect("reference reading failed.");
    // TODO: Remove duplication of gene and tax_id here - each contig is its own
    info!("Reading reference sequences in ..");
    let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(reference_reader)
        .expect("Failure to read transcripts in");
    info!("Building debruijn index ..");
    let index = build_index::<config::KmerType>(
        &seqs, &tx_names, &tx_gene_map
    ).expect("Failure to build index");
    info!("Finished building index!");

    // Do the mappings
    let reads = fastq::Reader::from_file(forward_fastq).expect("Failure to read reference sequences");
    let (eq_class_indices, eq_class_coverages) = pseudoaligner::process_reads::<config::KmerType>(reads, &index)
        .expect("Failure during mapping process");
    info!("Finished mapping reads!");
    debug!("eq_class_indices: {:?}, eq_class_coverages: {:?}", eq_class_indices, eq_class_coverages);

    let kmer_coverage_total: f64 = eq_class_coverages.iter().sum::<usize>() as f64;

    // EM algorithm
    // Randomly assign random relative abundances to each of the genomes
    let mut contig_to_relative_abundance = vec![1.0; seqs.len()]; // for dev, assign each the same abundance

    loop { // loop until converged
        // E-step: Determine the number of reads we expect come from each genome
        // given their relative abundance.

        // for each equivalence class / count pair, we expect for genome A the
        // abundance of A divided by the sum of abundances of genomes in the
        // equivalence class.
        let mut contig_to_read_count = vec![0.0; seqs.len()];
        // Make a new vec containing the eq_class memberships
        let mut eq_classes = vec![None; eq_class_coverages.len()];
        for (eq_class, index) in &eq_class_indices {
            eq_classes[*index] = Some(eq_class)
        }
        debug!("eq_classes: {:?}", eq_classes);
        for (i, coverage) in eq_class_coverages.iter().enumerate() {
            match eq_classes[i] {
                None => unreachable!(),
                Some(ref eqs) => {
                    let num_eqs = eqs.len() as f64;
                    for eq in *eqs {
                        // TODO: Remove the 'as usize' by making eq a usize throughout
                        contig_to_read_count[*eq as usize] += (*coverage as f64) / num_eqs;
                    }
                }
            }
        }
        debug!("After E-step have contig read counts: {:?}", contig_to_read_count);

        // M-step: Work out the relative abundance given the number of reads
        // predicted in the E-step. Relative abundance is just the ratio of the
        // coverages, weighted by the inverse of each contig's length.
        //TODO: zip here instead?
        // First determine the total scaled abundance
        let mut total_scaling_abundance: f64 = 0.0;
        for (i, read_count) in contig_to_read_count.iter().enumerate() {
            let to_add = read_count / (seqs[i].len() as f64);
            //debug!("Adding to_add {} for contig number {}", to_add, i);
            total_scaling_abundance += to_add;
        }
        // Next set the abundances so the total == 1.0
        //
        // Converge when all contigs with total abundance > 0.01 * kmer_coverage_total
        // change abundance by < 1%.
        let mut converge = true;
        for (i, read_count) in contig_to_read_count.iter().enumerate() {
            let to_add = read_count / (seqs[i].len() as f64);
            let new_relabund = to_add / total_scaling_abundance;
            if new_relabund * kmer_coverage_total > 0.01*100.0 { // Add 100 in there
                // as a rough mapped kmers per aligning fragment.

                // Enough abundance that this contig might stop convergence if it
                // changed by enough.
                let delta = new_relabund / contig_to_relative_abundance[i];
                debug!("For testing convergence, found delta {}", delta);
                if delta < 0.99 || delta > 1.01 {
                    debug!("No converge for you");
                    converge = false;
                }
            }
            contig_to_relative_abundance[i] = new_relabund;
        }
        debug!("At end of E-step, have coverages: {:?}", contig_to_relative_abundance);
        if converge { break; }
    }

    // Print results
    println!("contig\t{}", reference_path); //TODO: Use the same methods as elsewhere for printing.
    for (i, relabund) in contig_to_relative_abundance.iter().enumerate() {
        if print_zero_coverage_contigs || *relabund > 0.0 {
            println!("{}\t{}", tx_names[i], relabund*kmer_coverage_total / (seqs[i].len() as f64));
        }
    }
}
