use std::collections::BTreeMap;

use pseudoaligner::pseudoaligner::Pseudoaligner;
use debruijn::{Kmer, Vmer, Dir};
use debruijn::dna_string::DnaString;

// A region marked as being core for a clade
pub struct CoreGenomicRegion {
    pub clade_id: u32,
    pub contig_id: usize,
    pub start: u32,
    pub stop: u32,
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
    pub node_id_to_clade_cores: BTreeMap<usize, Vec<u32>>
}


/// core_genome_regions are Vecs of all the core genome regions in each genome
/// grouped together. Contig_id in each refers to the index of that contig in
/// the contig_sequences list.
pub fn generate_core_genome_pseudoaligner<K: Kmer + Send + Sync>(
    core_genome_regions: &Vec<Vec<CoreGenomicRegion>>,
    contig_sequences: &Vec<DnaString>,
    aligner: Pseudoaligner<K>,
) -> CoreGenomePseudoaligner<K> {

    let mut node_id_to_clade_cores: BTreeMap<usize, Vec<u32>> =
        BTreeMap::new();

    // Function to extract the next tranch of core genome regions for the next
    // contig
    let indices_of_current_contig =
        |regions: &Vec<CoreGenomicRegion>, starting_index: usize|
                       -> (usize, usize) {

            let target_contig_id = regions[starting_index].contig_id;
            let mut i = starting_index + 1;
            while i < regions.len() {
                if regions[i].contig_id == target_contig_id {
                    i += 1;
                } else {
                    break;
                }
            }
            (starting_index, i)
        };

    for genome_regions in core_genome_regions {
        let clade_id = genome_regions[0].clade_id;
        let (mut region_index_start, mut region_index_stop) =
            indices_of_current_contig(genome_regions, 0);

        // While there are more contig tranches, process them
        loop {
            let contig_id = genome_regions[region_index_start].contig_id;
            println!("Marking contig {} .. ", contig_id);
            let core_node_ids = thread_and_find_core_nodes(
                &aligner,
                &contig_sequences[contig_id],
                contig_id,
                &genome_regions[region_index_start..region_index_stop]);
            for nid in core_node_ids {
                match node_id_to_clade_cores.get_mut(&nid) {
                    Some(clade_cores) => {
                        if clade_cores.iter().find(|&r| *r==clade_id).is_none() {
                            clade_cores.push(clade_id)
                        };
                    },
                    None => {
                        node_id_to_clade_cores.insert(nid, vec![clade_id]);
                    }
                }
            }

            // Update for next iteration
            if region_index_stop < genome_regions.len() {
                let nexts = indices_of_current_contig(genome_regions, region_index_stop+1);
                region_index_start = nexts.0;
                region_index_stop = nexts.1;
            } else {
                // No more tranches
                break;
            }
        }
    }

    return CoreGenomePseudoaligner {
        index: aligner,
        core_genome_sizes: vec![], //TODO
        node_id_to_clade_cores: node_id_to_clade_cores,
    }
}

#[derive(Debug)]
struct GraphPosition {
    pub node_id: usize,
    pub offset: u32,
    pub is_forward: bool,
    pub contig_position: u32,
}

fn get_starting_position<K: Kmer + Send + Sync>(
    aligner: &Pseudoaligner<K>,
    contig: &DnaString,
    contig_id: usize,
) -> GraphPosition {

    let kmer_length = K::k();

    // Find the start of the contig in genome space
    //
    // Be careful here. Because of the MPHF, hashing false positives can
    // occur. So need to check that the first matching kmer really is that,
    // or whether it should be hashed in rc instead.
    let fwd_first_node = match aligner.dbg_index.get(
        &contig.get_kmer::<K>(0)) {

        Some((nid, offset)) => {
            let found_slice = aligner
                .dbg
                .get_node(*nid as usize)
                .sequence()
                .slice(*offset as usize,kmer_length)
                .to_string();
            if found_slice == contig.get_kmer::<K>(0).to_string() {
                println!("Found forward node for kmer {}", found_slice);
                Some(GraphPosition {
                    node_id: *nid as usize,
                    offset: *offset,
                    is_forward: true,
                    contig_position: 0,
                })
            } else {
                // Kmer hash mismatch
                println!("kmer hash doesn't end up pointing to the correct kmer");
                None
            }
        },
        None => None
    };
    return match fwd_first_node {
        // TODO Possible corner case if the contig starts and ends with the
        // same kmer?
        Some(x) => x,
        None => {
            let first_contig_kmer_rc = &contig.get_kmer::<K>(0).rc();
            match aligner.dbg_index.get(first_contig_kmer_rc) {
                None => {
                    error!("Unable to find starting node for contig number {} in the graph",
                           contig_id);
                    std::process::exit(1);
                },
                Some((nid, offset)) => {
                    let found_slice = aligner
                        .dbg
                        .get_node(*nid as usize)
                        .sequence()
                        .slice(*offset as usize,kmer_length)
                        .to_string();
                    if found_slice == first_contig_kmer_rc.to_string() {
                        println!("Found rc node {}", found_slice);
                        GraphPosition {
                            node_id: *nid as usize,
                            offset: *offset,
                            is_forward: false,
                            contig_position: 0,
                        }
                    } else {
                        // Kmer hash mismatch
                        error!("Unable to find starting node for contig number {} in the graph",
                               contig_id);
                        std::process::exit(1);
                    }
                }
            }
        }
    }
}

// Mark nodes as being core genome for a single contig. All core_regions should
// be from that contig. Return a list of nodes to be marked as core.
fn thread_and_find_core_nodes<K: Kmer + Send + Sync>(
    aligner: &Pseudoaligner<K>,
    contig_sequence: &DnaString,
    contig_id: usize,
    core_regions: &[CoreGenomicRegion])
-> Vec<usize> {

    let mut current_position = get_starting_position(
        &aligner, contig_sequence, contig_id);
    println!("Starting with position: {:?}", current_position);

    let kmer_length = K::k();
    let last_kmer_pos = contig_sequence.len() - kmer_length;
    println!("Found last kmer index {}", last_kmer_pos);

    let mut marked_nodes = vec![];
    let mut last_node_id = None;
    let mut current_core_region_idx = 0;

    println!("Found entire contig sequence {:?}", contig_sequence);

    for _ in 1..last_kmer_pos {
        println!("\n===== Finding kmer at contig position {}", current_position.contig_position);
        // if position is in core genome, mark it.
        let mut current_core_region = &core_regions[current_core_region_idx];
        if current_core_region.start <= current_position.contig_position {
            if current_core_region.stop <= current_position.contig_position {
                // Finished the current core genome region. Move to the next
                // one.
                current_core_region_idx += 1;
                if current_core_region_idx >= core_regions.len() {
                    // No more core regions for this genome
                    break;
                }
            }
            println!("In core region or just coming out of one ..");

            // If we are in the range of the next core region, add this node
            current_core_region = &core_regions[current_core_region_idx];
            if current_core_region.start <= current_position.contig_position {
                println!("Marking the current node {}", current_position.node_id);
                match last_node_id {
                    Some(nid) => {
                        if current_position.node_id != nid {
                            last_node_id = Some(current_position.node_id);
                            marked_nodes.push(current_position.node_id);
                        }
                    },
                    None => {
                        last_node_id = Some(current_position.node_id);
                        marked_nodes.push(current_position.node_id);
                    }
                }
            }
        }

        let current_contig_position = current_position.contig_position as usize;
        let target_kmer = contig_sequence.get_kmer::<K>(current_contig_position+1);
        next_position(
            &mut current_position,
            &aligner,
            &target_kmer);
        println!("next_position: {:?}", current_position);

        // Double check that the sequence now has the right kmer in that
        // position.
        let found_sequence = aligner
            .dbg
            .get_node(current_position.node_id as usize)
            .sequence();
        // Found_kmer is a DnaStringSlice
        let mut found_kmer = found_sequence
            .get_kmer::<K>(current_position.offset as usize);
        println!("Before potential rc(), forward found was {:?}", found_kmer);
        if !current_position.is_forward {
            println!("not is_forward");
            found_kmer = found_kmer.rc();
        }
        if found_kmer != target_kmer {
            println!("Kmer returned from search was incorrect!, expected {:?}, found {:?}",
                     target_kmer, found_kmer);
            std::process::exit(1);
        }
        println!("Found kmer was correct");
    }

    return marked_nodes;
}

/// Modify the position in place to point to the next position defined by the
/// given kmer.
fn next_position<K: Kmer + Send + Sync>(
    position: &mut GraphPosition,
    aligner: &Pseudoaligner<K>,
    kmer: &K) {

    let k = K::k();
    println!("Finding kmer {:?}", kmer);

    // If we are in the middle of the node, then just update the offset
    let current_node = aligner.dbg.get_node(position.node_id as usize);
    if position.is_forward && position.offset as usize+k+1 < current_node.len() {
        println!("Just going forward on the same node");
        position.offset += 1;
    } else if !position.is_forward && position.offset > 0 {
        println!("Just going reverse on the same node");
        position.offset -= 1;
    } else {
        let edges = match position.is_forward {
            true => current_node.r_edges(),
            false => current_node.l_edges()
        };
        println!("Found potential edges: {:?}", edges);
        let correct_edge = edges.iter().find(|edge| {
            let (target_node_id, incoming_side, _is_flipped) = (edge.0, edge.1, edge.2);
            let target_node = aligner.dbg.get_node(target_node_id);

            let new_kmer = match incoming_side {
                Dir::Left => target_node.sequence().get_kmer::<K>(0),
                Dir::Right => target_node.sequence().get_kmer::<K>(target_node.len()-k).rc()
            };
            println!("Testing new kmer {:?} from entire sequence {:?}", new_kmer, target_node.sequence());
            new_kmer == *kmer
        });
        match correct_edge {
            Some(edge) => {
                println!("Found the right edge: {:?}", edge);
                position.node_id = edge.0;
                match edge.1 {
                    Dir::Left => {
                        position.offset = 0;
                        position.is_forward = true;
                    },
                    Dir::Right => {
                        let target_node = aligner.dbg.get_node(position.node_id);
                        position.offset = (target_node.len()-k) as u32;
                        position.is_forward = false;
                    }
                }
            },
            None => {
                panic!("Did not find the right edge")
            }
        }
        println!("Got as");
    }
    position.contig_position += 1;
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
            "tests/data/2_single_species_dummy_dataset/2genomes/genomes.fna")
            .expect("reference reading failed.");
        info!("Reading reference sequences in ..");
        let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(reference_reader)
            .expect("Failure to read contigs file");
        info!("Building debruijn index ..");
        let index = build_index::build_index::<config::KmerType>(
            &seqs, &tx_names, &tx_gene_map, 1
        );
        let real_index = index.unwrap();
        println!("Graph was {:?}", &real_index.dbg);

        let core_aligner = generate_core_genome_pseudoaligner(
            &cores,
            &seqs,
            real_index
        );
        println!("done");

        println!("core_aligner.node_id_to_clade_cores: {:?}",
                 core_aligner.node_id_to_clade_cores);
        let mut expected = BTreeMap::new();
        expected.insert(2, vec![7]);
        expected.insert(4, vec![7]);
        assert_eq!(expected, core_aligner.node_id_to_clade_cores);
    }

    #[test]
    fn test_core_genome_2_core_regions() {
        let cores = vec![vec![
            CoreGenomicRegion {
                clade_id: 11,
                contig_id: 0,
                start: 0,
                stop: 1
            },
            CoreGenomicRegion {
                clade_id: 11,
                contig_id: 0,
                start: 80,
                stop: 82
            },
        ]];

        // Build index
        let reference_reader = fasta::Reader::from_file(
            "tests/data/2_single_species_dummy_dataset/2genomes/genomes.fna")
            .expect("reference reading failed.");
        info!("Reading reference sequences in ..");
        let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(reference_reader)
            .expect("Failure to read contigs file");
        info!("Building debruijn index ..");
        let index = build_index::build_index::<config::KmerType>(
            &seqs, &tx_names, &tx_gene_map, 1
        );
        let real_index = index.unwrap();
        println!("Graph was {:?}", &real_index.dbg);

        let core_aligner = generate_core_genome_pseudoaligner(
            &cores,
            &seqs,
            real_index
        );
        println!("done");

        println!("core_aligner.node_id_to_clade_cores: {:?}",
                 core_aligner.node_id_to_clade_cores);
        let mut expected = BTreeMap::new();
        expected.insert(5, vec![11]);
        expected.insert(4, vec![11]);
        assert_eq!(expected, core_aligner.node_id_to_clade_cores);
    }
}


