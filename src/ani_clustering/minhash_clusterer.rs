
use std::collections::{BTreeSet,BTreeMap};

use finch::distance::distance;
use finch::serialization::Sketch;

/// Given a list of genomes, return them clustered.
// TODO: Test whether this is a good enough procedure or if there are nasties
// e.g. failure to cluster bad quality genomes.
pub fn minhash_clusters(
    genomes: &[&str],
    min_ani: f32,
) -> Vec<Vec<usize>> {

    // Generate sketches for all input files
    let mut filter = finch::filtering::FilterParams { // dummy, no filtering is applied.
        filter_on: None,
        abun_filter: (None, None),
        err_filter: 0f32,
        strand_filter: 0f32,
    };
    info!("Sketching genomes for clustering ..");
    let sketches = finch::mash_files(
        genomes,
        1000,
        1000,
        21,
        &mut filter,
        true,
        0)
        .expect("Failed to create finch sketches for input genomes");
    info!("Finished sketching genomes for clustering.");

    let distance_threshold: f64 = (100.0 - min_ani as f64)/100.0;
    assert!(distance_threshold >= 0.0);
    assert!(distance_threshold < 1.0);

    // Greedily find reps
    let clusters = find_minhash_representatives(&sketches.sketches, distance_threshold);

    // Reassign non-reps based so they are assigned to the nearest
    // representative.
    return find_minhash_memberships(&clusters, &sketches.sketches);
}

/// Choose representatives, greedily assigning based on the min_ani threshold.
fn find_minhash_representatives(
    sketches: &[Sketch],
    ani_threshold: f64)
    -> BTreeSet<usize> {

    let mut to_return: BTreeSet<usize> = BTreeSet::new();

    for (i, sketch1) in sketches.iter().enumerate() {
        let mut is_rep = true;
        for j in &to_return {
            let sketch2: &Sketch = &sketches[*j];
            if distance(&sketch1.hashes, &sketch2.hashes, "", "", true)
                .expect("Failed to calculate distance by sketch comparison")
                .mashDistance
                <= ani_threshold {

                is_rep = false;
                break;
            }
        }
        if is_rep {
            to_return.insert(i);
        }
    }
    return to_return;
}

/// For each genome (sketch) assign it to the closest representative genome:
fn find_minhash_memberships(
    representatives: &BTreeSet<usize>,
    sketches: &[Sketch],
) -> Vec<Vec<usize>> {

    let mut rep_to_index = BTreeMap::new();
    for (i, rep) in representatives.iter().enumerate() {
        rep_to_index.insert(rep, i);
    }

    let mut to_return: Vec<Vec<usize>> = vec![vec![]; representatives.len()];
    for (i, sketch1) in sketches.iter().enumerate() {
        if representatives.contains(&i) {
            to_return[rep_to_index[&i]].push(i);
        } else {
            let mut best_rep_min_ani = None;
            let mut best_rep = None;
            for rep in representatives.iter() {
                let dist = distance(&sketch1.hashes, &sketches[*rep].hashes, "", "", true)
                    .expect("Failed to calculate distance by sketch comparison")
                    .mashDistance;
                if best_rep_min_ani.is_none() || dist < best_rep_min_ani.unwrap() {
                    best_rep = Some(rep);
                    best_rep_min_ani = Some(dist);
                }
            }
            to_return[rep_to_index[best_rep.unwrap()]].push(i);
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
    fn test_minhash_hello_world() {
        init();
        let clusters = minhash_clusters(
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
    fn test_minhash_two_clusters() {
        init();
        let clusters = minhash_clusters(
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
