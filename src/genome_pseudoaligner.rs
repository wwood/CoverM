use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use kmer_coverage::*;
use genomes_and_contigs::GenomesAndContigs;
use core_genome;
use core_genome::{CoreGenomePseudoaligner,CoreGenomicRegion};
use nucmer_core_genome_generator::nucmer_core_genomes_from_genome_fasta_files;

use pseudoaligner::*;
use debruijn::Kmer;
use debruijn::dna_string::DnaString;
use bio::io::{fasta,fastq};
use log::Level;
use csv;

/// Given a path to a file containing two columns (representative<tab>member),
/// where representative and member are paths to genome files, return a list of
/// groups, where the representative is the first in the 2nd dimension list. The
/// String data is not processed.
pub fn read_clade_definition_file(
    file_path: &str)
    -> Vec<Vec<String>> {

    // Open file as CSV

    let rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(Path::new(file_path));

    let mut current_rep: Option<String> = None;
    let mut current_members = vec![];
    let mut to_return = vec![];
    for result in rdr
        .expect(&format!("Failed to read clade definition file {}", file_path))
        .records() {
            let res = result.expect(&format!("Failed to read CSV record from clade definition file: {}", file_path));
            if res.len() != 2 {
                error!("Unexpectedly found a line in clade definition file with number of elements != 2: {:?}", res);
                std::process::exit(1);
            }
            let rep = &res[0];
            let member = &res[1];


            match current_rep {
                None => {
                    if rep != member {
                        error!("In clade definition file, the representative must be the first member, found {} and {}",
                               rep, member);
                        std::process::exit(1);
                    }
                    current_rep = Some(rep.to_string());
                    current_members.push(rep.to_string());
                },
                Some(ref previous_rep) => {
                    // Ensure that the first member is the rep
                    if rep == previous_rep {
                        current_members.push(member.to_string());
                    } else {
                        if rep != member {
                            error!("In clade definition file, the representative must be the first member, found {} and {}",
                                   rep, member);
                            std::process::exit(1);
                        }
                        to_return.push(current_members);
                        current_rep = Some(rep.to_string());
                        current_members = vec![rep.to_string()];
                    }
                }
            }
        }
    to_return.push(current_members);

    return to_return;
}

pub fn calculate_genome_kmer_coverage<K: Kmer + Sync + Send>(
    forward_fastq: &str,
    reverse_fastq: Option<&str>,
    num_threads: usize,
    print_zero_coverage_contigs: bool,
    core_genome_aligner: &CoreGenomePseudoaligner<K>,
    genomes_and_contigs: &GenomesAndContigs)
    -> Vec<(usize, f64)> {

    // Do the mappings
    let reads = fastq::Reader::from_file(forward_fastq)
        .expect(&format!("Failure to read file {}", forward_fastq));
    let reverse_reads = match reverse_fastq {
        Some(s) => Some(
            fastq::Reader::from_file(s)
                .expect(&format!("Failure to read reverse read file {}", s))),
        None => None
    };
    let (eq_class_indices, eq_class_coverages, _eq_class_counts) =
        pseudoaligner::process_reads::<K, CoreGenomePseudoaligner<K>>(
            reads, reverse_reads, &core_genome_aligner, num_threads)
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
    // process. TODO: Abstract this out so it isn't calculated for each sample.
    let genome_contigs = generate_genome_to_contig_indices_vec(
        &genomes_and_contigs,
        &core_genome_aligner.contig_names
    );
    let contig_idx_to_genome_idx: Vec<usize> = core_genome_aligner.contig_names.iter().map(
        |contig_name|
        genomes_and_contigs.contig_to_genome[contig_name]
    ).collect();
    let genome_lengths: &Vec<usize> = &core_genome_aligner.core_genome_sizes;
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
        let genome_coverages: Vec<f64> = genome_to_read_count.iter().zip(genome_lengths.iter()).map(
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
    let mut to_return = vec!();
    debug!("genome_to_read_count: {:?}", genome_to_read_count);
    for (i, total_coverage) in genome_to_read_count.iter().enumerate() {
        if print_zero_coverage_contigs || *total_coverage > 0.0 {
            // Print average coverage as total coverage divided by contig length.
            to_return.push((i, total_coverage / (genome_lengths[i] as f64)));
        }
    }
    info!("Finished printing genome coverages");
    return to_return;
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

fn report_core_genome_sizes(
    nucmer_core_genomes: &Vec<Vec<Vec<CoreGenomicRegion>>>,
    clades: &Vec<Vec<String>>) {

    let mut total_core_genome_size = 0u64;
    let mut minimum_core_genome_size = None;
    for (clade_i, core_genomes) in nucmer_core_genomes.iter().enumerate() {
        // minimum core genome size in the clade
        let minimum_of_clade = core_genomes
            .iter()
            .map(
                |core_genome_regions|
                core_genome_regions
                    .iter()
                    .fold(0, |acc,c| acc+c.stop-c.start))
            .min()
            .unwrap();
        total_core_genome_size += minimum_of_clade as u64;
        minimum_core_genome_size = match minimum_core_genome_size {
            None => Some((clade_i, minimum_of_clade)),
            Some((prev_i, prev_min)) => {
                if minimum_of_clade < prev_min {
                    Some((clade_i, minimum_of_clade))
                } else {
                    Some((prev_i,prev_min))
                }
            }
        }
    }
    info!("Found minimum core genome size {} from clade {}",
          minimum_core_genome_size.unwrap().1,
          clades[minimum_core_genome_size.unwrap().0][0]
    );
    info!("Found mean core genome size {}",
          total_core_genome_size as f64 / nucmer_core_genomes.len() as f64);
}

/// Read in each contig from each genome of each clade, given a path to their
/// fasta files.
fn read_clade_genome_strings(
    clades: &Vec<Vec<String>>)
    -> Vec<Vec<Vec<DnaString>>> {

    clades
        .iter()
        .map(
            |genome_fastas|
            genome_fastas
                .iter()
                .map(
                |genome_fasta| {
                    let fasta_reader = fasta::Reader::from_file(genome_fasta)
                        .expect(&format!("Error reading genome fasta file {}", genome_fasta));
                    fasta_reader
                        .records()
                        .map(
                            |r| {
                                let r2 = r.expect(&format!(
                                    "Error parsing fasta file entry in {}", genome_fasta));
                                // TODO: This conversion introduces randomness.
                                // How to treat N characters? Use
                                // DnaString#from_acgt_bytes_hashn ?
                                DnaString::from_acgt_bytes(
                                    r2.seq()
                                )
                            }
                        )
                        .collect()
                }
            ).collect()
        )
        .collect()
}

pub fn core_genome_coverage_pipeline<K: Kmer + Send + Sync>(
    read_inputs: &Vec<PseudoalignmentReadInput>,
    num_threads: usize,
    print_zero_coverage_contigs: bool,
    index: DebruijnIndex<K>,
    genomes_and_contigs: &GenomesAndContigs,
    clades: &Vec<Vec<String>>) {

    assert!(clades.len() > 0);

    // Write GFA TODO: debug
    // if log_enabled!(Level::Debug) { 
    //     info!("Writing GFA file ..");
    //     let mut gfa_writer = std::fs::File::create("/tmp/my.gfa").unwrap();
    //     index.index.dbg.write_gfa(&mut gfa_writer).unwrap();
    // }

    // For each clade, nucmer against the first genome.
    info!("Calculating core genomes ..");
    // TODO: ProgressBar?
    let nucmer_core_genomes: Vec<Vec<Vec<CoreGenomicRegion>>> = clades
        .iter()
        .enumerate()
        .map(
            |(i, clade_fastas)|
            nucmer_core_genomes_from_genome_fasta_files(
                &clade_fastas.iter().map(|s| &**s).collect::<Vec<&str>>()[..],
                i as u32)
        ).collect();
    info!("Finished calculating core genomes");

    // Check core genome sizes / report
    report_core_genome_sizes(&nucmer_core_genomes, clades);

    debug!("Found core genomes: {:#?}", nucmer_core_genomes);

    // Thread genomes recording the core genome nodes
    // TODO: These data are at least sometimes read in repeatedly, when they
    // maybe should just be cached or something.
    info!("Reading in genome FASTA files to thread graph");
    let dna_strings = read_clade_genome_strings(clades);
    info!("Threading DeBruijn graph");
    let core_genome_pseudoaligner = core_genome::generate_core_genome_pseudoaligner(
        &nucmer_core_genomes,
        &dna_strings,
        index
    );

    debug!("Found node_to_core_genomes: {:#?}",
           &core_genome_pseudoaligner.node_id_to_clade_cores);
    // if log_enabled!(Level::Debug) {
    //     // Write CSV data to be loaded into bandage
    //     use std::io::Write;
    //     let mut csv_writer = std::fs::File::create("/tmp/my.core_nodes.csv").unwrap();
    //     writeln!(csv_writer, "Node,Clades").unwrap();
    //     for (node_id, clades) in &core_genome_pseudoaligner.node_id_to_clade_cores {
    //         writeln!(csv_writer, "{},\"{:?}\"", node_id, clades).unwrap();
    //     }
    // }

    // Map / EM / Print
    println!("Sample\tGenome\tCoverage");
    for read_input in read_inputs {
        let covs = calculate_genome_kmer_coverage(
            &read_input.forward_fastq,
            match read_input.reverse_fastq {
                Some(ref s) => Some(&s),
                None => None
            },
            num_threads,
            print_zero_coverage_contigs,
            &core_genome_pseudoaligner,
            genomes_and_contigs);

        for res in covs {
            println!(
                "{}\t{}\t{}",
                read_input.sample_name,
                genomes_and_contigs.genomes[res.0],
                res.1);
        }
        info!("Finished printing genome coverages for sample {}",
              read_input.sample_name);
    }
    info!("Finished printing contig coverages");
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_read_clade_file() {
        let mut tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();

        writeln!(tf, "/path/g1.fna\t/path/g1.fna").unwrap();
        writeln!(tf, "/path/g1.fna\t/path/g3.fna").unwrap();
        writeln!(tf, "/path/g1.fna\t/path/g2.fna").unwrap();
        writeln!(tf, "/path/g10.fna\t/path/g10.fna").unwrap();
        writeln!(tf, "/path/g10.fna\t/path/g20.fna").unwrap();

        tf.flush().unwrap();
        assert_eq!(
            vec![vec!["/path/g1.fna","/path/g3.fna","/path/g2.fna"], vec!["/path/g10.fna","/path/g20.fna"]],
            read_clade_definition_file(tf.path().to_str().unwrap()));
    }
}
