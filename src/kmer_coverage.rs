use pseudoaligner::build_index::build_index;
use pseudoaligner::*;

use bio::io::{fasta, fastq};

pub fn calculate_genome_kmer_coverage(
    reference_path: &str,
    forward_fastq: &str) {
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

    // Print out the coverages divided by their length
    println!("Contig\t{}", reference_path); //TODO: Use the same methods as elsewhere for printing.
    for (i, coverage) in eq_class_coverages.iter().enumerate() {
        println!("{}\t{}", tx_names[i], *coverage as f64 / seqs[i].len() as f64)
    }
}
