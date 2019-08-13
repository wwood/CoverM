// use std;
// use std::io::Read;
// use std::path::Path;
use std::collections::BTreeMap;

use run_command_safely;
use nucmer_runner;

use core_genome::CoreGenomicRegion;
use tempdir::TempDir;
use bio::io::fasta;

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

#[derive(Debug, PartialOrd, PartialEq, Clone, Eq, Ord)]
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

#[derive(Debug)]
struct ContigAlignments<'a> {
    ref_contig: &'a str,
    alignments: Vec<&'a nucmer_runner::NucmerDeltaAlignment>
}


pub fn nucmer_core_genomes_from_genome_fasta_files(
    genome_fasta_paths: &[&str],
    clade_id: u32
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

    let reference_contig_name_to_id = read_reference_contig_names(&reference_fasta);
    debug!("Read {} contigs from reference", reference_contig_name_to_id.len());

    debug!("Generating nucmer alignments ..");
    let mut all_alignments = generate_nucmer_alignments(
        &reference_suffix_prefix_str,
        &reference_fasta,
        &genome_fasta_paths[1..genome_fasta_paths.len()]);

    //ref_contig index to all aligned regions (whole thing is None for first)
    let mut ref_aligned_regions: Option<Vec<Option<RTree<RefAlignedSection>>>> = None;
    //ref_contig index to poisoned regions
    let mut ref_poisoned_regions: Vec<Option<RTree<RefAlignedSection>>> =
        vec![None; reference_contig_name_to_id.len()];

    // Nucmer each non-first genome against the reference in a new function
    for mut alignments in all_alignments.iter_mut() {
        // Split each of the alignments into per-reference-contig chunks
        let alignment_chunks = split_alignments_by_ref_contig(&mut alignments);

        // insert the reference coordinates into the tree
        let mut regions_poisoned_by_query = vec![vec![]; reference_contig_name_to_id.len()]; //ref_contig index to position pairs
        for contig_alignments in alignment_chunks.iter() {
            let reference_contig_index = reference_contig_name_to_id[contig_alignments.ref_contig];
            debug!("Determining alignments against contig {:?}", contig_alignments);
            let aligned_sections: Vec<AlignedSection> = contig_alignments
                .alignments
                .iter()
                .enumerate()
                .map(|(i, aln)|
                     AlignedSection {
                         start: aln.seq1_start as i64,
                         stop: aln.seq1_stop as i64,
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
                regions_poisoned_by_query[reference_contig_index] =
                    poisoned_regions_in_this_contig;
            }
        }
        debug!("Found poisoned regions: {:?}", regions_poisoned_by_query);

        // Finished processing all contigs in this query genome

        // Set aligned regions to be the intersection of the regions seen
        // aligned here
        ref_aligned_regions = Some(
            intersect(
                &ref_aligned_regions,
                &alignment_chunks,
                &reference_contig_name_to_id));
        debug!("ref_aligned_regions remaining: {:?}", ref_aligned_regions);

        // Set poisoned regions to be the union of the regions poisoned
        // previously and those poisoned by this genome.
        ref_poisoned_regions = union(
            &ref_poisoned_regions,
            &regions_poisoned_by_query);
    }

    let ref_alignments = generate_final_ref_core_genome_regions(
        clade_id,
        ref_aligned_regions,
        ref_poisoned_regions);

    let mut query_alignments = match_query_alignments_to_ref_alignments(
        clade_id,
        &reference_contig_name_to_id,
        &genome_fasta_paths[1..],
        &ref_alignments,
        &all_alignments);

    query_alignments.insert(0, ref_alignments);
    return query_alignments;
}

fn read_reference_contig_names(
    fasta_path: &str)
    -> BTreeMap<String, usize> {

    let reader = fasta::Reader::from_file(fasta_path)
        .expect(&format!("Failed to open reference fasta file {}", fasta_path));
    let mut tree = BTreeMap::new();
    for (i, record) in reader.records().enumerate() {
        tree.insert(record.expect("Failed to parse ref fasta").id().to_string(), i);
    }
    return tree;
}

fn generate_nucmer_alignments(
    reference_suffix_prefix_str: &str,
    reference_fasta: &str,
    query_fasta_paths: &[&str]
) -> Vec<Vec<nucmer_runner::NucmerDeltaAlignment>> {

    return query_fasta_paths
        .iter()
        .map(|q| nucmer_runner::nucmer_to_deltas(
            reference_suffix_prefix_str,
            reference_fasta,
            q))
        .collect()
}

fn intersect(
    aligned_regions1: &Option<Vec<Option<RTree<RefAlignedSection>>>>,
    aligned_chunks: &Vec<ContigAlignments>,
    reference_contig_name_to_id: &BTreeMap<String, usize>)
    -> Vec<Option<RTree<RefAlignedSection>>> {

    return match aligned_regions1 {
        None => {
            // No previous regions. Just add the regions in the aligned_chunks.
            let mut regions1 = vec![None; reference_contig_name_to_id.len()];
            for contig_alignments in aligned_chunks {
                let ref_contig_id = reference_contig_name_to_id[contig_alignments.ref_contig];
                regions1[ref_contig_id] = Some(
                    RTree::bulk_load(
                        contig_alignments
                            .alignments
                            .iter()
                            .map(|a| RefAlignedSection {
                                start: a.seq1_start as i64,
                                stop: a.seq1_stop as i64,
                            })
                            .collect()))
            }
            regions1
        },
        Some(to_return) => {
            let mut regions_vec = vec![None; reference_contig_name_to_id.len()];
            for contig_alignments in aligned_chunks {
                let ref_contig_id = reference_contig_name_to_id[contig_alignments.ref_contig];

                // Intersect the alignments for the current region against the
                // previous alignments.
                match &to_return[ref_contig_id] {
                    None => {
                        // No previous alignments in this contig, so ignore.
                    },
                    Some(prev_rtree) => {
                        let mut new_alignments = vec![];
                        for alignment in &contig_alignments.alignments {
                            for intersect in prev_rtree.locate_in_envelope_intersecting(
                                &rstar::AABB::from_corners(
                                    [alignment.seq1_start as i64,0],
                                    [alignment.seq1_stop as i64,0])) {
                                new_alignments.push(
                                    RefAlignedSection {
                                        start: std::cmp::max(alignment.seq1_start as i64, intersect.start),
                                        stop: std::cmp::min(alignment.seq1_stop as i64, intersect.stop),
                                    }
                                )
                            }
                        }
                        regions_vec[ref_contig_id] = Some(RTree::bulk_load(new_alignments))
                    }
                }
            }
            regions_vec
        }
    }
}

fn union(
    ref_poisoned_regions: &Vec<Option<RTree<RefAlignedSection>>>,
    regions_poisoned_by_query: &Vec<Vec<(i64, i64)>>)
    -> Vec<Option<RTree<RefAlignedSection>>> {

    // The two vectors should have the same length, so zip iter.
    return ref_poisoned_regions.iter().zip(regions_poisoned_by_query.iter())
        .map(|(prev_poisons_opt, cur_poisons)| {
            // Given two sorted lists of coordinates, return their union, for each contig
            match prev_poisons_opt {
                None => {
                    if cur_poisons.is_empty() {
                        None
                    } else {
                        let mut current_start = cur_poisons[0].0;
                        let mut current_stop = cur_poisons[0].1;
                        let mut poisoned_regions = vec![];
                        for coord in cur_poisons {
                            if current_start <= coord.0 {
                                // If fully contained, do nothing
                                // If start is within, but stop is beyond, extend
                                current_stop = std::cmp::max(current_stop, coord.1);
                            } else {
                                // If start is beyond, save last and start anew
                                poisoned_regions.push(RefAlignedSection {
                                    start: current_start,
                                    stop: current_stop
                                });
                                current_start = coord.0;
                                current_stop = coord.1;
                            }
                        }
                        // Push the last
                        poisoned_regions.push(RefAlignedSection {
                            start: current_start,
                            stop: current_stop
                        });
                        Some(RTree::bulk_load(poisoned_regions))
                    }
                },
                Some(_prev_rtree) => panic!()
            }
        })
        .collect();
}

// Generate the core genome regions, which is the union of the
// aligned sections minus the poisoned regions.
fn generate_final_ref_core_genome_regions(
    clade_id: u32,
    ref_aligned_regions: Option<Vec<Option<RTree<RefAlignedSection>>>>,
    ref_poisoned_regions: Vec<Option<RTree<RefAlignedSection>>>)
    -> Vec<CoreGenomicRegion> {

    let mut core_genome_regions = vec![];
    ref_aligned_regions
        .expect("Found no aligned regions for this species, giving up.")
        .iter()
        .zip(ref_poisoned_regions.iter())
        .enumerate()
        .for_each( |(contig_idx, (aligned_opt, poison_opt))| {
            match aligned_opt {
                None => {}, // No alignments in this contig.
                Some(alignment_rtree) => {
                    // TODO: It seems like iterating over the rtree in sorted
                    // order should be possible, but I've not thought too deeply
                    // about it.
                    let mut aligns_sorted: Vec<&RefAlignedSection> = alignment_rtree.iter().collect();
                    aligns_sorted.sort_unstable();

                    // We do not need to flatten the alignments because marking
                    // something as core is just a yes/no thing later on in this
                    // whole process.
                    for align in aligns_sorted.iter() {
                        match poison_opt {
                            None => {
                                core_genome_regions.push(CoreGenomicRegion {
                                    clade_id: clade_id,
                                    contig_id: contig_idx,
                                    start: align.start as u32,
                                    stop: align.stop as u32,
                                })
                            },
                            Some(poison_rtree) => {
                                // Find intersecting poison alignments, and sort them
                                let mut poisons: Vec<&RefAlignedSection> = poison_rtree
                                    .locate_in_envelope_intersecting(
                                        &rstar::AABB::from_corners([align.start,0],[align.stop,0]))
                                    .collect();
                                poisons.sort_unstable();

                                // Add the remaining parts
                                let mut current_align_start = align.start;
                                for poison in poisons {
                                    // If align start < poison_start, log an
                                    // aligned section, and reset align_start
                                    if current_align_start < poison.start {
                                        core_genome_regions.push (CoreGenomicRegion {
                                            clade_id: clade_id,
                                            contig_id: contig_idx,
                                            start: current_align_start as u32,
                                            stop: poison.start as u32 - 1,
                                        });
                                        current_align_start = poison.stop + 1;
                                    } else {
                                        // Else set next aligned start to be max
                                        // of current_poison.stop+1 and
                                        // current_align_start. Have to include
                                        // the latter because this poison might
                                        // be fully contained within the last
                                        // one.
                                        current_align_start = std::cmp::max(
                                            poison.stop + 1,
                                            current_align_start)
                                    }
                                }
                                // Add the final aligned section if one exists
                                if current_align_start < align.stop {
                                    core_genome_regions.push(CoreGenomicRegion {
                                        clade_id: clade_id,
                                        contig_id: contig_idx,
                                        start: current_align_start as u32,
                                        stop: align.stop as u32 - 1,
                                    })
                                }
                            }
                        }
                    }
                }
            }
        });
    return core_genome_regions;
}

/// Given core regions defined in the reference sequence, return the analogous
/// regions from the query alignments.
fn match_query_alignments_to_ref_alignments(
    clade_id: u32,
    reference_contig_name_to_id: &BTreeMap<String, usize>,
    query_fasta_paths: &[&str],
    ref_alignments: &Vec<CoreGenomicRegion>,
    all_alignments: &Vec<Vec<nucmer_runner::NucmerDeltaAlignment>>
) -> Vec<Vec<CoreGenomicRegion>> {
    let num_reference_contigs = reference_contig_name_to_id.len();
    let mut ref_contig_tmp: Vec<Vec<&CoreGenomicRegion>> = vec![vec![]; num_reference_contigs];
    for ref_align in ref_alignments {
        ref_contig_tmp[ref_align.contig_id].push(ref_align);
    }
    let ref_contig_rtrees: Vec<RTree<RefAlignedSection>> = ref_contig_tmp
        .iter()
        .map(|aligns| RTree::bulk_load(
            aligns
                .iter()
                .map(|a| RefAlignedSection {
                    start: a.start as i64,
                    stop: a.stop as i64,
                })
                .collect()
        ))
        .collect();

    // Iterate each alignment fro the query
    return all_alignments
        .iter()
        .enumerate()
        .map( |(query_index, query_genome_alignments)| {
            let mut aligning_regions = vec![];
            let query_contig_name_to_idx = read_reference_contig_names(query_fasta_paths[query_index]);
            for query_align in query_genome_alignments {
                let ref_contig_id = reference_contig_name_to_id[&query_align.seq1_name];
                for aligning_ref_region in ref_contig_rtrees[ref_contig_id]
                    .locate_in_envelope_intersecting(
                        &rstar::AABB::from_corners(
                            [query_align.seq1_start as i64,0],
                            [query_align.seq1_stop as i64,0])) {
                        // debug!("Adding core genome derived from ref {:#?} and query {:#?}",
                        //        aligning_ref_region, query_align);

                        aligning_regions.push(CoreGenomicRegion {
                            clade_id: clade_id,
                            contig_id: query_contig_name_to_idx[&query_align.seq2_name],
                            start: query_align.query_coord_at(
                                std::cmp::max(
                                    aligning_ref_region.start as u32,
                                    query_align.seq1_start,
                                )
                            ),
                            stop: query_align.query_coord_at(
                                std::cmp::min(
                                    aligning_ref_region.stop as u32,
                                    query_align.seq1_stop,
                                )
                            ),
                        })
                    }
            }
            aligning_regions
        })
        .collect()
}


fn generate_nucmer_suffixes(fasta: &str, prefix: &std::path::Path) {
    let mut cmd = std::process::Command::new("nucmer");
    cmd
        .arg(&format!("--save={}", prefix.to_str().expect("Failed to convert path to str")))
        .arg(fasta);
    run_command_safely(cmd, "nucmer --save");
}

// Split alignments by contig, then return alignments ordered by ref_contig,
// then ref_start, then undefined.
fn split_alignments_by_ref_contig<'a>(
    alignments: &'a mut Vec<nucmer_runner::NucmerDeltaAlignment>)
    -> Vec<ContigAlignments<'a>> {

    alignments.sort_unstable_by(|a,b| a.partial_cmp(b).unwrap());

    let mut prev_contig = &alignments[0].seq1_name;
    let mut chunks: Vec<ContigAlignments> = vec![];
    let mut prev_chunk_id: usize = 0;
    for (i, aln) in alignments.iter().enumerate() {
        if i == 0 {
            chunks.push(ContigAlignments {
                ref_contig: prev_contig,
                alignments: vec![aln]
            });
        } else {
            if aln.seq1_name != *prev_contig {
                prev_chunk_id += 1;
                prev_contig = &aln.seq1_name;
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
    fn test_nucmer_generate_core_hello_world_dummy() {
        init();
        let parsed = nucmer_core_genomes_from_genome_fasta_files(
            &["tests/data/2_single_species_dummy_dataset/2genomes/parsnp/g2.fasta",
              "tests/data/2_single_species_dummy_dataset/2genomes/parsnp/g1_duplicate.fasta",
              "tests/data/2_single_species_dummy_dataset/2genomes/parsnp/g1_split.fasta",
              ],
            8
        );
        println!("{:#?}", parsed);
        assert_eq!(vec![
            vec![
                CoreGenomicRegion {
                    clade_id: 8,
                    contig_id: 0,
                    start: 0,
                    stop: 86,
                },
                CoreGenomicRegion {
                    clade_id: 8,
                    contig_id: 0,
                    start: 87,
                    stop: 199,
                },
            ],
            vec![
                CoreGenomicRegion {
                    clade_id: 8,
                    contig_id: 0,
                    start: 88,
                    stop: 200,
                },
                CoreGenomicRegion {
                    clade_id: 8,
                    contig_id: 0,
                    start: 1,
                    stop: 87,
                },
            ],
            vec![
                CoreGenomicRegion {
                    clade_id: 8,
                    contig_id: 1,
                    start: 0,
                    stop: 86,
                },
                CoreGenomicRegion {
                    clade_id: 8,
                    contig_id: 0,
                    start: 0,
                    stop: 112,
                },
            ]
        ], parsed);
    }

    #[test]
    fn test_nucmer_generate_core_hello_world_real() {
        init();
        let parsed = nucmer_core_genomes_from_genome_fasta_files(
            &["tests/data/parsnp/1_first_group/73.20120800_S1D.21.fna",
              "tests/data/parsnp/1_first_group/73.20120600_E3D.30.fna"],
            8
        );
    }
}
