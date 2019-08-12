// use std;
// use std::io::Read;
// use std::path::Path;
use std::collections::BTreeMap;

use run_command_safely;
use nucmer_runner;

use core_genome::CoreGenomicRegion;
use tempdir::TempDir;
// use bio::io::fasta;

use rstar::RTree;

#[derive(Debug, PartialOrd, PartialEq)]
struct AlignedSection {
    pub start: i64, // i64 since u32 isn't supported by rstar
    pub stop: i64,
    pub alignment_id: usize
}

impl rstar::RTreeObject for AlignedSection {
    // Hack here - always set the y dimension to 0, then the impl compiles.
    type Envelope = rstar::AABB<[i64; 2]>;

    fn envelope(&self) -> Self::Envelope {
        rstar::AABB::from_corners([self.start, 0], [self.stop, 0])
    }
}

#[derive(Debug, PartialOrd, PartialEq)]
struct RefAlignedSection {
    pub start: i64, // i64 since u32 isn't supported by rstar
    pub stop: i64,
}

impl rstar::RTreeObject for RefAlignedSection {
    // Hack here - always set the y dimension to 0, then the impl compiles.
    type Envelope = rstar::AABB<[i64; 2]>;

    fn envelope(&self) -> Self::Envelope {
        rstar::AABB::from_corners([self.start, 0], [self.stop, 0])
    }
}


pub fn nucmer_core_genomes_from_genome_fasta_files(
    genome_fasta_paths: &[&str],
    _clade_id: u32
) -> Vec<Vec<CoreGenomicRegion>> {

    // Use the first genome passed as the reference
    let reference_fasta = genome_fasta_paths[0];
    // Generate a suffix file for the reference
    // Make a tempdir
    let tmp_dir = TempDir::new("coverm-nucmer")
        .expect("Failed to make temporary directory for nucmer");
    let reference_suffix_prefix = tmp_dir.path().join("reference_suffixes");
    generate_nucmer_suffixes(reference_fasta, &reference_suffix_prefix);
    let reference_suffix_prefix_str = reference_suffix_prefix.to_str()
        .expect("Failed to convert prefix path to str");

    let reference_contig_name_to_id = read_reference_contig_names();

    //ref_contig name to all aligned regions
    let mut ref_aligned_regions: Option<BTreeMap<String, Vec<RTree<RefAlignedSection>>>>;
    //ref_contig name to poisoned regions
    let mut ref_poisoned_regions: BTreeMap<String, Vec<RTree<RefAlignedSection>>>;

    // Nucmer each non-first genome against the reference in a new function
    for (i, query_fasta) in genome_fasta_paths.iter().enumerate() {
        if i==0 { continue } // Don't compare reference to itself

        // Read the .delta file hits. The format is described at
        // http://mummer.sourceforge.net/manual/#nucmer
        let alignments = nucmer_runner::run_nucmer_and_show_coords(
            reference_suffix_prefix_str,
            reference_fasta,
            query_fasta);
        debug!("Without parsing, nucmer generated {} alignments for {}",
               alignments.len(), query_fasta);
        // Accept hits > 250bp
        let mut passable_alignments: Vec<&nucmer_runner::NucmerAlignment> = alignments
            .iter()
            .filter(|f| f.ref_length() >= 250 && f.identity > 80.)
            .collect();
        debug!("After removing low quality alignments, found {} alignments for {}",
               alignments.len(), query_fasta);

        // Split each of the alignments into per-reference-contig chunks
        let alignment_chunks = split_alignments_by_ref_contig(&mut passable_alignments);

        // insert the reference coordinates into the tree
        let mut regions_poisoned_by_query = BTreeMap::new(); //ref_contig name to position pairs
        for contig_alignments in alignment_chunks {
            if contig_alignments.ref_contig != "73.20120800_S1D.21_contig_8557" {
                continue; //debug
            }
            debug!("Determining alignments against contig {:?}", contig_alignments);
            let aligned_sections: Vec<AlignedSection> = contig_alignments
                .alignments
                .iter()
                .enumerate()
                .map(|(i, aln)|
                     AlignedSection {
                         start: aln.ref_start as i64,
                         stop: aln.ref_stop as i64,
                         alignment_id: i
                     })
                .collect();

            let rtree = RTree::bulk_load(aligned_sections);

            // Once all alignments are present, determine stretches that are covered
            let mut poisoned_regions_in_this_contig = vec![];
            for ase in rtree.iter() {
                debug!("Looking at ase: {:?}", ase);
                for intersect in rtree.locate_in_envelope_intersecting(
                    &rstar::AABB::from_corners([ase.start,0],[ase.stop,0])) {

                    // Use > not != here to avoid duplication
                    if ase.alignment_id > intersect.alignment_id {
                        poisoned_regions_in_this_contig.push(
                            (std::cmp::max(ase.start, intersect.start),
                             std::cmp::min(ase.stop, intersect.stop)));
                        debug!("Found poisonous intersect {:?}", intersect);
                    }
                }
            }

            if poisoned_regions_in_this_contig.len() > 0 {
                regions_poisoned_by_query.insert(
                    contig_alignments.ref_contig,
                    poisoned_regions_in_this_contig);
            }

            break;
        }
        debug!("Found poisoned regions: {:?}", regions_poisoned_by_query);

        // Finished processing all contigs in this query genome

        // Set aligned regions to be the intersection of the regions seen
        // aligned here
        ref_aligned_regions = Some(
            intersect(
                ref_aligned_regions,
                alignment_chunks,
                reference_contig_name_to_id));

        // Set poisoned regions to be the union of the regions poisoned
        // previously and those poisoned by this genome.
        ref_poisoned_regions = union(
            ref_poisoned_regions, regions_poisoned_by_query);

        break;
    }

    return generate_final_core_genomes(ref_aligned_regions, ref_poisoned_regions);
}

fn intersect(
    aligned_regions1: Option<BTreeMap<usize, Vec<RTree<RefAlignedSection>>>>,
    aligned_chunks: Vec<ContigAlignments>,
    reference_contig_name_to_id: BTreeMap<String, usize>)
    -> Option<BTreeMap<String, Vec<RTree<RefAlignedSection>>>> {

    match aligned_regions1 {
        None => {
            // No previous regions. Just add the regions in the aligned_chunks.
            for chunk in aligned_chunks {
                let ref_contig_id = reference_contig_name_to_it[aligned_chunks.ref_contig];

            }
        }
    }
}

// Generate the core genome regions, which is the union of the
// aligned sections minus the poisoned regions.
fn generate_final_core_genomes(
    ref_aligned_regions: Option<BTreeMap<String, Vec<RTree<RefAlignedSection>>>>,
    ref_poisoned_regions: BTreeMap<String, Vec<RTree<RefAlignedSection>>>)
    -> Vec<Vec<CoreGenomicRegion>> {
    panic!()
}

fn generate_nucmer_suffixes(fasta: &str, prefix: &std::path::Path) {
    let mut cmd = std::process::Command::new("nucmer");
    cmd
        .arg(&format!("--save={}", prefix.to_str().expect("Failed to convert path to str")))
        .arg(fasta);
    run_command_safely(cmd, "nucmer --save");
}

#[derive(Debug)]
struct ContigAlignments<'a> {
    ref_contig: &'a str,
    alignments: Vec<&'a nucmer_runner::NucmerAlignment>
}

struct RefAlignedRegion {
    start: u32,
    stop: u32,
}

// Split alignments by contig, then return alignments ordered by ref_contig,
// then ref_start, then undefined.
fn split_alignments_by_ref_contig<'a>(
    alignments: &'a mut Vec<&nucmer_runner::NucmerAlignment>)
    -> Vec<ContigAlignments<'a>> {

    alignments.sort_unstable_by(|a,b| a.partial_cmp(b).unwrap());

    let mut prev_contig = &alignments[0].ref_contig;
    let mut chunks: Vec<ContigAlignments> = vec![];
    let mut prev_chunk_id: usize = 0;
    for (i, aln) in alignments.iter().enumerate() {
        if i == 0 {
            chunks.push(ContigAlignments {
                ref_contig: prev_contig,
                alignments: vec![aln]
            });
        } else {
            if aln.ref_contig != *prev_contig {
                prev_chunk_id += 1;
                prev_contig = &aln.ref_contig;
                chunks.push(ContigAlignments {
                    ref_contig: prev_contig,
                    alignments: vec![aln]
                })
            } else {
                chunks[prev_chunk_id].alignments.push(aln)
            }
        }
    }
    return chunks;
}



#[cfg(test)]
mod tests {
    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn test_nucmer_generate_core_hello_world() {
        init();
        let parsed = nucmer_core_genomes_from_genome_fasta_files(
            &["tests/data/parsnp/1_first_group/73.20120800_S1D.21.fna",
              "tests/data/parsnp/1_first_group/73.20120600_E3D.30.fna"],
            8
        );
    }
}
