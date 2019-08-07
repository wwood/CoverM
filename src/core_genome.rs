use std::collections::{BTreeMap};

use pseudoaligner::pseudoaligner::Pseudoaligner;
use debruijn::{Kmer, Vmer, Dir};
use debruijn::dna_string::DnaString;

// A region marked as being core for a clade
pub struct CoreGenomicRegion {
    pub clade_id: usize,
    pub contig_id: usize,
    pub start: usize,
    pub stop: usize
}

/// Represent a Pseudoaligner that has extra annotations, specifically, some
/// regions are marked as being 'core' for a given clade, and so it is more
/// reliable to calculate abundance just based only on these regions.
pub struct CoreGenomePseudoaligner<K: Kmer> {
    pub index: Pseudoaligner<K>,

    /// Core genome size of each clade
    pub core_genome_sizes: Vec<usize>,

    /// Map of node_id in graph to list of clades where those nodes are in that
    /// clade's core. Nodes not in any core are not stored.
    pub node_id_to_clade_cores: BTreeMap<usize, Vec<usize>>
}


/// core_genome_regions are Vecs of all the core genome regions in each genome
/// grouped together. Contig_id in each refers to the index of that contig in
/// the contig_sequences list.
pub fn generate_core_genome_pseudoaligner<K: Kmer + Send + Sync>(
    core_genome_regions: &Vec<Vec<CoreGenomicRegion>>,
    contig_sequences: &Vec<DnaString>,
    aligner: Pseudoaligner<K>,
) -> CoreGenomePseudoaligner<K> {

    let node_id_to_clade_cores: BTreeMap<usize, Vec<usize>> =
        BTreeMap::new();

    // numbers above original_eq_class_len mark core genome regions, possibly
    for regions in core_genome_regions {
        let contig_id = regions[0].contig_id;
        // Find the start of the genome in the graph
        let contig = &contig_sequences[contig_id];
        //let read_length = contig.len();
        let mut kmer_pos: usize = 0;
        //let kmer_length = K::k();
        //let last_kmer_pos = read_length - kmer_length;

        // TODO: Be careful here. Because of the MPHF, hashing false positives
        // can occur. So need to check that the first matching kmer really is
        // that, or whether it should be hashed in rc instead.

        while kmer_pos <= 6 {//}last_kmer_pos {
            println!("\nDetermining kmer at position {}", kmer_pos);
            let kmer = contig.get_kmer::<K>(kmer_pos);
            println!("Got {:?}", kmer);

            match aligner.dbg_index.get(&kmer) {
                Some((nid, offset)) => {
                    println!("found kmer node {:?} offset {:?}", nid, offset);
                    let node = aligner.dbg.get_node(*nid as usize);
                    println!("seq of node: {:?}", node.sequence());
                    let exts = node.exts();
                    println!("exts of node: {:?}", exts);
                    // r_edges returns
                    //
                    // (target_node id, incoming side of target node, whether
                    // target node has is flipped)
                    let edges = node.r_edges();
                    println!("edges of node: {:?}", edges);

                    let node_r = aligner.dbg.get_node(edges[0].0);
                    println!("seq of node: {:?}", node_r.sequence());
                    let exts = node_r.exts();
                    println!("exts of node: {:?}", exts);
                    let edges = node_r.r_edges();
                    println!("edges of node: {:?}", edges);
                },
                None => {
                    println!("Faild");
                    println!("rc: {:?}", aligner.dbg_index.get(&kmer.rc()));
                    error!("Unable to find start of genome in Debruijn graph, when it should be there");
                    //std::process::exit(1);
                }
            };
            kmer_pos += 1;
        }
    }

    return CoreGenomePseudoaligner {
        index: aligner,
        core_genome_sizes: vec![], //TODO
        node_id_to_clade_cores: node_id_to_clade_cores,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pseudoaligner::*;
    use bio::io::fasta;

    #[test]
    fn test_core_genome_hello_world() {
        let cores = vec![vec![
            CoreGenomicRegion {
                clade_id: 7,
                contig_id: 0,
                start: 10,
                stop: 100
            }
        ]];

        // Build index
        let reference_reader = fasta::Reader::from_file(
            "tests/data/2_single_species_dummy_dataset/2genomes/genomes.fna").expect("reference reading failed.");
        info!("Reading reference sequences in ..");
        let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(reference_reader)
            .expect("Failure to read contigs file");
        info!("Building debruijn index ..");
        let index = build_index::build_index::<config::KmerType>(
            &seqs, &tx_names, &tx_gene_map, 1
        );
        let real_index = index.unwrap();
        println!("Graph was {:?}", &real_index.dbg);

        let _cores = generate_core_genome_pseudoaligner(
            &cores,
            &seqs,
            real_index
        );
        println!("done");
    }
}


