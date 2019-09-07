use std::fs::File;
use std::io::BufWriter;

use pseudoaligner::*;
use pseudoaligner::build_index::build_index;
use genomes_and_contigs::GenomesAndContigs;
use bio::io::fasta;
use bincode;
use debruijn::dna_string::DnaString;
use serde::{Serialize, Deserialize, de::DeserializeOwned};
use debruijn::Kmer;

#[derive(Serialize, Deserialize)]
pub struct DebruijnIndex<K>
where K: debruijn::Kmer {
    pub index: pseudoaligner::Pseudoaligner<K>,
    pub seq_lengths: Vec<usize>,
    pub tx_names: Vec<String>,
}

/// Read a path into a DebruijnIndex, where each contig is considered
/// independently, and not e.g. associated with a genome.
pub fn generate_debruijn_index_without_groupings<K: Kmer + Sync + Send>(
    reference_path: &str,
    num_threads: usize)
    -> DebruijnIndex<K> {
    
    let reference_reader = fasta::Reader::from_file(reference_path).expect("reference reading failed.");
    // TODO: Remove duplication of gene and tax_id here - each contig is its own
    info!("Reading reference sequences in ..");
    let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(
        reference_reader,
        |contig_name| { (contig_name.to_string(), contig_name.to_string()) })
        .expect("Failure to read contigs file");

    return generate_debruijn_index(
        num_threads,
        &seqs,
        tx_names,
        &tx_gene_map);
}

/// Read a path into a DebruijnIndex, where each contig is considered
/// independently, and not e.g. associated with a genome.
pub fn generate_debruijn_index_grouping_via_separator<K: Kmer + Sync + Send>(
    reference_path: &str,
    separator: char,
    num_threads: usize)
    -> DebruijnIndex<K> {
    
    let reference_reader = fasta::Reader::from_file(reference_path).expect("reference reading failed.");
    // TODO: Remove duplication of gene and tax_id here - each contig is its own
    info!("Reading reference sequences in ..");
    let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(
        reference_reader,
        |contig_name| {
            debug!("target name {:?}, separator {:?}", contig_name, separator);
            let offset = contig_name.find(separator).expect(
                &format!("Contig name {} does not contain split symbol, so cannot determine which genome it belongs to",
                        contig_name));
            (contig_name[(0..offset)].to_string(), contig_name.to_string()) })
        .expect("Failure to read contigs file");

    return generate_debruijn_index(
        num_threads,
        &seqs,
        tx_names,
        &tx_gene_map);
}

pub fn generate_debruijn_index_grouping_via_genomes_and_contigs<K: Kmer + Sync + Send>(
    genomes_and_contigs: &GenomesAndContigs,
    reference_path: &str,
    num_threads: usize
) -> DebruijnIndex<K> {

    let reference_reader = fasta::Reader::from_file(reference_path).expect("reference reading failed.");
    // TODO: Remove duplication of gene and tax_id here - each contig is its own
    info!("Reading reference sequences in ..");
    let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(
        reference_reader,
        |contig_name|
            (contig_name.to_string(), 
            genomes_and_contigs.genome_of_contig(&contig_name.to_string())
            .expect(&format!("Contig name {} was not associated with any genome", contig_name))
            .to_string()))
        .expect("Failure to read contigs file");

    return generate_debruijn_index(
        num_threads,
        &seqs,
        tx_names,
        &tx_gene_map);
}

pub fn generate_debruijn_index<'a, K: Kmer + Sync + Send>(
    num_threads: usize,
    seqs: &[DnaString],
    sequence_names: Vec<String>,
    sequence_name_to_grouping: &std::collections::HashMap<String,String>)
    -> DebruijnIndex<K> {

    info!("Building debruijn index ..");
    let index = build_index::<K>(
        seqs, &sequence_names, sequence_name_to_grouping, num_threads
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
                tx_names: sequence_names,
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
