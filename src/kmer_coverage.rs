use std::fs::File;
use std::io::BufWriter;

use pseudoaligner::build_index::build_index;
use pseudoaligner::*;
use debruijn::Kmer;

use bio::io::{fasta, fastq};
use bincode;
use serde::{Serialize, Deserialize, de::DeserializeOwned};

use log::Level;

#[derive(Serialize, Deserialize)]
pub struct DebruijnIndex<K>
where K: debruijn::Kmer {
    pub index: pseudoaligner::Pseudoaligner<K>,
    pub seq_lengths: Vec<usize>,
    pub tx_names: Vec<String>,
}

pub fn generate_debruijn_index<'a, K: Kmer + Sync + Send>(
    reference_path: &str,
    num_threads: usize)
    -> DebruijnIndex<K> {

    // Build index
    let reference_reader = fasta::Reader::from_file(reference_path).expect("reference reading failed.");
    // TODO: Remove duplication of gene and tax_id here - each contig is its own
    info!("Reading reference sequences in ..");
    let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(reference_reader)
        .expect("Failure to read contigs file");
    info!("Building debruijn index ..");
    let index = build_index::<K>(
        &seqs, &tx_names, &tx_gene_map, num_threads
    );

    match index {
        Err(e) => {
            error!("Error generating debruijn index: {:?}", e);
            std::process::exit(1);
        },
        Ok(true_index) => {
            info!("Successfully built DeBruijn index");
            return DebruijnIndex {
                index: true_index,
                seq_lengths: seqs.iter().map(|s| s.len()).collect(),
                tx_names: tx_names,
            }
        }
    }
}

pub fn save_index<K>(
    index: DebruijnIndex<K>,
    output_file: &str)
where K: Kmer + Serialize {

    let f = File::create(output_file).unwrap();
    let writer = BufWriter::new(f);
    let mut snapper = snap::Writer::new(writer);
    info!("Writing index to {} ..", output_file);
    bincode::serialize_into(&mut snapper, &index)
        .expect("Failure to serialize or write index");
    info!("Finished writing index");
}

pub fn restore_index<'a, K: Kmer + DeserializeOwned>(
    saved_index_path: &'a str)
    -> DebruijnIndex<K> {

    let f = File::open(saved_index_path).expect("file not found");
    let mut unsnapper = snap::Reader::new(f);
    debug!("Deserialising DB ..");
    return bincode::deserialize_from(&mut unsnapper)
        .expect("Error reading previously saved index - perhaps it was \
                 generated with a different version of CoverM?");
}


pub fn calculate_contig_kmer_coverage<K>(
    forward_fastq: &str,
    reverse_fastq: Option<&str>,
    num_threads: usize,
    print_zero_coverage_contigs: bool,
    index: DebruijnIndex<K>)
where K: Kmer + Sync + Send {

    // Do the mappings
    let reads = fastq::Reader::from_file(forward_fastq)
        .expect(&format!("Failure to read file {}", forward_fastq));
    let reverse_reads = match reverse_fastq {
        Some(s) => Some(
            fastq::Reader::from_file(s)
                .expect(&format!("Failure to read reverse read file {}", s))),
        None => None
    };
    let (eq_class_indices, eq_class_coverages, _eq_class_counts) = pseudoaligner::process_reads(
        reads, reverse_reads, &index.index, num_threads)
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

    // EM algorithm
    // Randomly assign random relative abundances to each of the genomes
    // TODO: remove: for dev, assign each the same abundance
    info!("Starting EM process ..");
    let mut contig_to_relative_abundance = vec![1.0; index.seq_lengths.len()];
    let mut contig_to_read_count;
    let mut num_covergence_steps: u32 = 0;

    loop { // loop until converged
        // E-step: Determine the number of reads we expect come from each contig
        // given their relative abundance.

        // for each equivalence class / count pair, we expect for genome A the
        // abundance of A divided by the sum of abundances of genomes in the
        // equivalence class.
        contig_to_read_count = vec![0.0; index.seq_lengths.len()];
        for (i, coverage) in eq_class_coverages.iter().enumerate() {
            match eq_classes[i] {
                None => unreachable!(),
                Some(ref eqs) => {
                    // Find coverages of each contig to add
                    let relative_abundances: Vec<f64> = eqs.iter().map(
                        |eq|
                        // TODO: Remove the 'as usize' by making eq a usize throughout
                        contig_to_relative_abundance[*eq as usize]).collect();

                    // Add coverages divided by the sum of relative abundances
                    let total_abundance_of_matching_contigs: f64 = relative_abundances.iter().sum();
                    for (eq, relabund) in eqs.iter().zip(relative_abundances.iter()) {
                        contig_to_read_count[*eq as usize] += (*coverage as f64) *
                            relabund / total_abundance_of_matching_contigs;
                    }
                }
            }
        }
        debug!("After E-step have contig read counts: {:?}", contig_to_read_count);

        // M-step: Work out the relative abundance given the number of reads
        // predicted in the E-step. Relative abundance is just the ratio of the
        // coverages, weighted by the inverse of each contig's length.
        //
        // Or, for genome mode, weighted by the inverse of the sum of the
        // genome's length.
        //TODO: zip here instead?
        // First determine the total scaled abundance
        let mut total_scaling_abundance: f64 = 0.0;
        let mut converge = true;

        for (i, read_count) in contig_to_read_count.iter().enumerate() {
            total_scaling_abundance += read_count / (index.seq_lengths[i] as f64)
        }

        // Next set the abundances so the total == 1.0
        //
        // Converge when all contigs with total abundance > 0.01 * kmer_coverage_total
        // change abundance by < 1%.
        for (i, read_count) in contig_to_read_count.iter().enumerate() {
            let to_add = read_count / (index.seq_lengths[i] as f64);
            let new_relabund = to_add / total_scaling_abundance;
            if new_relabund * kmer_coverage_total > 0.01*100.0 { // Add 100 in there
                // as a rough mapped kmers per aligning fragment.

                // Enough abundance that this contig might stop convergence if it
                // changed by enough.
                let delta = new_relabund / contig_to_relative_abundance[i];
                debug!("For testing convergence of index {}, found delta {}", i, delta);
                if delta < 0.99 || delta > 1.01 {
                    debug!("No converge for you");
                    converge = false;
                }
            }
            contig_to_relative_abundance[i] = new_relabund;
        }
        debug!("At end of M-step, have relative abundances: {:?}", contig_to_relative_abundance);

        num_covergence_steps += 1;
        if converge {
            info!("EM process converged after {} steps", num_covergence_steps);
            break;
        }
    }

    // Print results
    println!("contig\tkmer"); //TODO: Use the same methods as elsewhere for printing.
    for (i, total_coverage) in contig_to_read_count.iter().enumerate() {
        if print_zero_coverage_contigs || *total_coverage > 0.0 {
            // Print average coverage as total coverage divided by contig length.
            println!("{}\t{}", index.tx_names[i], total_coverage / (index.seq_lengths[i] as f64));
        }
    }
    info!("Finished printing contig coverages");
}
