
use std::io::{BufReader, Read, Write};
use std::collections::{BTreeMap,BTreeSet};

use finish_command_safely;

use tempfile;
use csv;

/// Given a list of genomes, return them clustered.
// TODO: Test whether this is a good enough procedure or if there are nasties
// e.g. failure to cluster bad quality genomes.
pub fn minhash_clusters(
    genomes: &[&str],
    min_ani: f32,
) -> Vec<Vec<usize>> {

    // Generate sketches for all input files
    let mut filter = finch::filtering::FilterParams { // dummy, no filtering is applied.
        filter_on: false,
        abun_filter: (None, None),
        err_filter: 0f32,
        strand_filter: 0f32,
    }
    info!("Sketching genomes for clustering ..");
    finch::mash_files(
        genomes,
        1000,
        1000,
        21,
        &mut filter,
        true,
        0,
    );
    info!("Finished sketching genomes for clustering.");


    // Mash each against to 

    // TODO: start here....... aligns =





    // Write genomes to tempfile
    let mut genomes_file = tempfile::Builder::new()
        .prefix("coverm-fastani-tmp")
        .tempfile()
        .expect("Failed to make tempfile");
    let mut genome_order = BTreeMap::new();
    for (i, genome) in genomes.iter().enumerate() {
        writeln!(genomes_file, "{}", genome)
            .expect("Failed to write genome fasta path to tempfile");
        genome_order.insert(genome.to_string(), i);
    }
    genomes_file.flush().expect("Failed to flush file input to fastANI");
    let genomes_file_str = &genomes_file.path().to_str()
        .expect("tempfile path to str failed");

    // FastANI to stdout
    let mut cmd = std::process::Command::new("fastANI");
    cmd
        .arg("-o")
        .arg("/dev/stdout")
        .arg("--queryList")
        .arg(&genomes_file_str)
        .arg("--refList")
        .arg(&genomes_file_str)
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    debug!("Running fastANI command: {:?}", &cmd);
    let mut process = cmd.spawn().expect(&format!("Failed to spawn {}", "fastANI"));
    let stdout = process.stdout.as_mut().unwrap();
    let stdout_reader = BufReader::new(stdout);
    let mut aligns = parse_fastani(stdout_reader, &genome_order, min_ani);
    finish_command_safely(process, "fastANI");
    aligns.sort_unstable_by(|a,b| a.partial_cmp(b).unwrap());

    debug!("Found alignments from fastANI: {:#?}", aligns);
    std::fs::copy(genomes_file.path(), std::path::Path::new("/tmp/in")).unwrap();

    // 2 pass clustering - first choose reps greedily
    let first_reps = choose_reps(&aligns);
    debug!("Found rep genome ids: {:?}", first_reps);

    // Second stage of clustering - re-assign each genome to its nearest rep
    aligns.sort_unstable_by(|a,b| b.ani.partial_cmp(&a.ani).unwrap());
    return assign_to_reps(first_reps, &aligns, genomes.len());
}

fn assign_to_reps(
    first_reps: Vec<usize>,
    aligns: &Vec<GenomeAlignment>, // must be sorted by ani largest to smallest
    num_genomes: usize,
) -> Vec<Vec<usize>> {

    let mut first_reps_idx = BTreeSet::new();
    for rep in &first_reps {
        first_reps_idx.insert(rep);
    }

    let rep_alignments: Vec<&GenomeAlignment> = aligns
        .iter()
        .filter(
            |align|
            (first_reps_idx.contains(&align.ref_genome_index) || first_reps_idx.contains(&align.query_genome_index)))
        .collect();
    debug!("After removing non-rep hits, alignments found: {:#?}", rep_alignments);

    let mut clusters: Vec<Vec<usize>> = vec![vec![]; first_reps.len()];
    let mut singletons = vec![];

    // For each genome, linear search the aligns. The first encounter with a rep is the best.
    // TODO: This seems substandard. Is it actually slow though?
    for i in 0..num_genomes {
        match first_reps_idx.get(&i) {
            Some(_) => {
                // This genome is a rep. Add it to its own list
                clusters[i].push(i)
            },
            None => {
                match rep_alignments.iter().find(
                    |align|
                    (align.ref_genome_index == i || align.query_genome_index == i)
                ) {
                    Some(align) => {
                        if align.ref_genome_index == i {
                            clusters[align.query_genome_index].push(i)
                        } else {
                            clusters[align.ref_genome_index].push(i)
                        }
                    },
                    None => {
                        // This genome matches no rep.
                        singletons.push(i)
                    }
                }
            }
        }
    }

    let mut singleton_vec: Vec<Vec<usize>> = singletons.into_iter().map(|s| vec![s]).collect();
    clusters.append(&mut singleton_vec);

    return clusters;
}

fn choose_reps(
    aligns: &Vec<GenomeAlignment>)
    -> Vec<usize> {

    let mut member_to_rep: BTreeMap<usize, usize> = BTreeMap::new();

    for align in aligns {
        // If r has been seen, ignore
        if member_to_rep.get(&align.ref_genome_index).is_some() {
            continue;
        }

        // Else if q has been seen, ignore (may not be right but we are greedy)
        if member_to_rep.get(&align.query_genome_index).is_some() {
            continue;
        }

        member_to_rep.insert(align.query_genome_index, align.ref_genome_index);
    }

    let mut reps: Vec<usize> = member_to_rep.values().map(|a| *a).collect();
    reps.sort_unstable();
    reps.dedup();

    return reps;
}

#[derive(PartialEq, PartialOrd, Debug)]
struct GenomeAlignment {
    ref_genome_index: usize,
    query_genome_index: usize,
    ani: f32, // Must be last in the struct so derive of PartialEq works
}

fn parse_fastani<R: Read>(
    reader: R,
    genome_order: &BTreeMap<String, usize>,
    min_ani: f32,
) -> Vec<GenomeAlignment> {

    let mut to_return = vec![];

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(reader);

    for record_res in rdr.records() {
        match record_res {
            Ok(record) => {
                assert!(record.len() == 5);
                let ani = record[2].parse().expect("Failed to convert fastani ANI to float value");
                if ani >= min_ani {
                    let ref_genome = genome_order[&record[0].to_string()];
                    let query_genome = genome_order[&record[1].to_string()];
                    if ref_genome < query_genome {
                        to_return.push(GenomeAlignment {
                            ref_genome_index: ref_genome,
                            query_genome_index: query_genome,
                            ani: ani
                        })
                    }
                }
            },
            Err(e) => {
                error!("Error parsing fastani output: {}", e);
                std::process::exit(1);
            }
        }
    }
    return to_return;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn test_fastani_hello_world() {
        init();
        let clusters = fastani_clusters(
            &["tests/data/parsnp/1_first_group/73.20120800_S1X.13.fna",
              "tests/data/parsnp/1_first_group/73.20120600_S2D.19.fna",
              "tests/data/parsnp/1_first_group/73.20120700_S3X.12.fna",
              "tests/data/parsnp/1_first_group/73.20110800_S2D.13.fna",
            ],
            95.0
        );
        assert_eq!(
            vec![vec![0,1,2,3]],
            clusters
        )
    }

    #[test]
    fn test_fastani_two_clusters() {
        init();
        let clusters = fastani_clusters(
            &["tests/data/parsnp/1_first_group/73.20120800_S1X.13.fna",
              "tests/data/parsnp/1_first_group/73.20120600_S2D.19.fna",
              "tests/data/parsnp/1_first_group/73.20120700_S3X.12.fna",
              "tests/data/parsnp/1_first_group/73.20110800_S2D.13.fna",
            ],
            98.0
        );
        assert_eq!(
            vec![vec![0,1,3],vec![2]],
            clusters
        )
    }
}
