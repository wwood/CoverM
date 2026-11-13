use std::convert::TryFrom;
use std::io::Write;
use std::process;

use OutputWriter;

pub enum CoverageTakerType {
    SingleFloatCoverageStreamingCoveragePrinter {
        print_stream: OutputWriter,
        current_stoit: Option<String>,
    },
    PileupCoverageCoveragePrinter {
        print_stream: OutputWriter,
        current_stoit: Option<String>,
        current_entry: Option<String>,
    },
    CachedSingleFloatCoverageTaker {
        stoit_names: Vec<String>,
        entry_names: Vec<Option<String>>,
        // per-stoit recorded coverages
        coverages: Vec<StoitCoverages>,
        current_stoit_index: Option<usize>,
        current_entry_index: Option<usize>,
        num_coverages: usize, // number of different coverage calculations
    },
}

pub trait CoverageTaker {
    fn start_stoit(&mut self, stoit_name: &str);
    fn start_entry(&mut self, entry_order_id: usize, entry_name: &str);
    fn add_single_coverage(&mut self, coverage: f32);
    // This function is only used with
    // CoverageEstimator::PileupCountsGenomeCoverageEstimator, and is a bit of a
    // hack to have it in this generalised trait, really.
    fn add_coverage_entry(&mut self, num_reads: usize, num_bases: u64);
    fn finish_entry(&mut self);
}

impl CoverageTakerType {
    pub fn new_single_float_coverage_streaming_coverage_printer(
        print_stream: OutputWriter,
    ) -> CoverageTakerType {
        CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter {
            print_stream,
            current_stoit: None,
        }
    }
    pub fn new_pileup_coverage_coverage_printer(print_stream: OutputWriter) -> CoverageTakerType {
        CoverageTakerType::PileupCoverageCoveragePrinter {
            print_stream,
            current_stoit: None,
            current_entry: None,
        }
    }
    pub fn new_cached_single_float_coverage_taker(num_coverages: usize) -> CoverageTakerType {
        CoverageTakerType::CachedSingleFloatCoverageTaker {
            stoit_names: vec![],
            entry_names: vec![],
            coverages: vec![],
            current_stoit_index: None,
            current_entry_index: None,
            num_coverages,
        }
    }
}

impl CoverageTaker for CoverageTakerType {
    fn start_stoit(&mut self, stoit_name: &str) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter {
                print_stream: _,
                ref mut current_stoit,
            } => *current_stoit = Some(stoit_name.to_owned()),
            CoverageTakerType::PileupCoverageCoveragePrinter {
                print_stream: _,
                ref mut current_stoit,
                current_entry: _,
            } => *current_stoit = Some(stoit_name.to_owned()),
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                ref mut stoit_names,
                entry_names: _,
                ref mut coverages,
                ref mut current_stoit_index,
                current_entry_index: _,
                num_coverages,
            } => {
                // Release over-allocation from the previous stoit, and size this
                // one to match it, since stoits usually have the same entries.
                let expected_num_entries = match coverages.last_mut() {
                    Some(previous) => {
                        previous.shrink_to_fit();
                        previous.len()
                    }
                    None => 0,
                };
                stoit_names.push(stoit_name.to_owned());
                coverages.push(StoitCoverages::with_capacity(
                    expected_num_entries,
                    *num_coverages,
                ));
                *current_stoit_index = Some(stoit_names.len() - 1);
            }
        }
    }

    fn start_entry(&mut self, entry_order_id: usize, entry_name: &str) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter {
                print_stream,
                current_stoit,
            } => match current_stoit {
                Some(stoit) => write!(print_stream, "{stoit}\t{entry_name}").unwrap(),
                None => unreachable!(),
            },
            CoverageTakerType::PileupCoverageCoveragePrinter {
                print_stream: _,
                current_stoit: _,
                ref mut current_entry,
            } => {
                *current_entry = Some(entry_name.to_owned());
            }
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                stoit_names: _,
                ref mut entry_names,
                ref mut coverages,
                ref current_stoit_index,
                ref mut current_entry_index,
                ..
            } => {
                debug!("Starting an entry with ID {entry_order_id} and name {entry_name}");
                // if the first time this entry has been seen, record its name.
                if entry_order_id >= entry_names.len() {
                    debug!(
                        "Adding 1 to entry order id {}, where names size is {}",
                        entry_order_id,
                        entry_names.len()
                    );
                    entry_names.resize(entry_order_id + 1, None)
                }
                if entry_names[entry_order_id].is_none() {
                    entry_names[entry_order_id] = Some(entry_name.to_owned());
                }
                match &entry_names[entry_order_id] {
                    Some(prev) => {
                        if prev != entry_name {
                            error!(
                                "Found a difference amongst the reference sets used for \
                                    mapping. For this (non-streaming) usage of CoverM, all \
                                    BAM files must have the same set of reference sequences. \
                                    Previous entry was {prev}, new is {entry_name}"
                            );
                            process::exit(1);
                        }
                    }
                    None => unreachable!(),
                }
                coverages[current_stoit_index.unwrap()].push_entry(entry_order_id);
                *current_entry_index = Some(entry_order_id);
            }
        }
    }

    fn add_single_coverage(&mut self, coverage: f32) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter {
                print_stream, ..
            } => {
                if coverage == 0.0 {
                    write!(print_stream, "\t0").unwrap()
                } else {
                    write!(print_stream, "\t{coverage}").unwrap()
                }
            }
            CoverageTakerType::PileupCoverageCoveragePrinter { .. } => {
                unreachable!();
            }
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                stoit_names: _,
                entry_names: _,
                ref mut coverages,
                ref current_stoit_index,
                ..
            } => coverages[current_stoit_index.unwrap()]
                .values
                .push(coverage),
        }
    }

    // This function is only used with
    // CoverageEstimator::PileupCountsGenomeCoverageEstimator, and is a bit of a
    // hack, really.
    fn add_coverage_entry(&mut self, num_reads: usize, num_bases: u64) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter { .. }
            | CoverageTakerType::CachedSingleFloatCoverageTaker { .. } => unreachable!(),
            CoverageTakerType::PileupCoverageCoveragePrinter {
                ref mut print_stream,
                ref mut current_stoit,
                ref mut current_entry,
            } => {
                let stoit = match current_stoit {
                    Some(ref stoit) => stoit,
                    None => unreachable!(),
                };
                let entry = match current_entry {
                    Some(ref entry) => entry,
                    None => unreachable!(),
                };
                writeln!(print_stream, "{stoit}\t{entry}\t{num_reads}\t{num_bases}").unwrap();
            }
        }
    }

    fn finish_entry(&mut self) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter {
                print_stream, ..
            } => writeln!(print_stream).unwrap(),
            CoverageTakerType::PileupCoverageCoveragePrinter { .. }
            | CoverageTakerType::CachedSingleFloatCoverageTaker { .. } => {}
        }
    }
}

/// Coverages recorded for one stoit by a CachedSingleFloatCoverageTaker,
/// stored as `num_coverages` consecutive f32 values per recorded entry. While
/// entries arrive as 0, 1, 2, .. (the usual case when zero coverage entries
/// are printed) their indices are implicit, so only the values are stored. A
/// u32 index per entry is kept once an entry is skipped.
#[derive(PartialEq, Debug, Default)]
pub struct StoitCoverages {
    entry_indices: Option<Vec<u32>>,
    num_entries: usize,
    values: Vec<f32>,
}

impl StoitCoverages {
    fn with_capacity(num_entries: usize, num_coverages: usize) -> StoitCoverages {
        StoitCoverages {
            entry_indices: None,
            num_entries: 0,
            values: Vec::with_capacity(num_entries * num_coverages),
        }
    }

    fn push_entry(&mut self, entry_index: usize) {
        // Printing walks each stoit with a StoitCoverageCursor, which relies
        // on entries being recorded in strictly increasing order. Fail loudly
        // rather than print coverages against the wrong entry.
        if !self.is_empty() {
            let previous = self.entry_index(self.num_entries - 1);
            assert!(
                entry_index > previous,
                "Coverage entries must be recorded in increasing order, but entry {entry_index} \
                followed entry {previous}. Please report this as a bug in CoverM."
            );
        }
        let entry_index_u32 = u32::try_from(entry_index)
            .expect("Too many reference sequences / genomes to cache coverage values");
        match self.entry_indices {
            Some(ref mut indices) => indices.push(entry_index_u32),
            None => {
                if entry_index != self.num_entries {
                    let mut indices: Vec<u32> = (0..self.num_entries as u32).collect();
                    indices.push(entry_index_u32);
                    self.entry_indices = Some(indices);
                }
            }
        }
        self.num_entries += 1;
    }

    fn shrink_to_fit(&mut self) {
        self.values.shrink_to_fit();
        if let Some(ref mut indices) = self.entry_indices {
            indices.shrink_to_fit();
        }
    }

    /// Number of entries recorded for this stoit.
    pub fn len(&self) -> usize {
        self.num_entries
    }

    pub fn is_empty(&self) -> bool {
        self.num_entries == 0
    }

    /// The entry index of the i'th recorded entry.
    pub fn entry_index(&self, i: usize) -> usize {
        match &self.entry_indices {
            Some(indices) => indices[i] as usize,
            None => i,
        }
    }

    /// The coverage values of the i'th recorded entry.
    pub fn coverages(&self, i: usize, num_coverages: usize) -> &[f32] {
        &self.values[i * num_coverages..(i + 1) * num_coverages]
    }
}

/// Walks the entries recorded for one stoit in entry order, returning the
/// coverages for each requested entry, or None if the stoit has no record of
/// it (i.e. its coverages are all zero). Requested entry indices must be
/// increasing.
pub struct StoitCoverageCursor<'a> {
    stoit: &'a StoitCoverages,
    num_coverages: usize,
    next: usize,
}

impl<'a> StoitCoverageCursor<'a> {
    pub fn new(stoit: &'a StoitCoverages, num_coverages: usize) -> StoitCoverageCursor<'a> {
        StoitCoverageCursor {
            stoit,
            num_coverages,
            next: 0,
        }
    }

    pub fn coverages_for(&mut self, entry_index: usize) -> Option<&'a [f32]> {
        while self.next < self.stoit.len() && self.stoit.entry_index(self.next) < entry_index {
            self.next += 1;
        }
        if self.next < self.stoit.len() && self.stoit.entry_index(self.next) == entry_index {
            self.next += 1;
            Some(self.stoit.coverages(self.next - 1, self.num_coverages))
        } else {
            None
        }
    }
}

/// Indices of entries observed in any stoit, in increasing order. These are
/// the rows printed from a CachedSingleFloatCoverageTaker, with zeros filled in
/// for stoits which have no record of a given entry.
pub fn observed_entry_indices(entry_names: &[Option<String>]) -> impl Iterator<Item = usize> + '_ {
    entry_names
        .iter()
        .enumerate()
        .filter(|(_, name)| name.is_some())
        .map(|(i, _)| i)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn recorded(stoit: &StoitCoverages, num_coverages: usize) -> Vec<(usize, Vec<f32>)> {
        (0..stoit.len())
            .map(|i| {
                (
                    stoit.entry_index(i),
                    stoit.coverages(i, num_coverages).to_vec(),
                )
            })
            .collect()
    }

    #[test]
    fn test_cached_hello_world() {
        let mut c = CoverageTakerType::new_cached_single_float_coverage_taker(2);
        c.start_stoit("stoit1");
        c.start_entry(0, "contig1");
        c.add_single_coverage(1.1);
        c.add_single_coverage(1.2);
        match c {
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                stoit_names,
                entry_names,
                coverages,
                current_stoit_index,
                current_entry_index,
                num_coverages,
            } => {
                assert_eq!(vec!("stoit1".to_string()), stoit_names);
                assert_eq!(vec!(Some("contig1".to_string())), entry_names);
                assert_eq!(1, coverages.len());
                assert_eq!(vec![(0, vec![1.1, 1.2])], recorded(&coverages[0], 2));
                // Contiguous entries do not need their indices stored
                assert_eq!(None, coverages[0].entry_indices);
                assert_eq!(0, current_stoit_index.unwrap());
                assert_eq!(0, current_entry_index.unwrap());
                assert_eq!(2, num_coverages);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_cached_two_samples_matching() {
        let mut c = CoverageTakerType::new_cached_single_float_coverage_taker(2);
        c.start_stoit("stoit1");
        c.start_entry(0, "contig1");
        c.add_single_coverage(1.1);
        c.add_single_coverage(1.2);
        c.start_entry(3, "contig2");
        c.add_single_coverage(2.1);
        c.add_single_coverage(2.2);
        c.start_stoit("stoit2");
        c.start_entry(0, "contig1");
        c.add_single_coverage(10.1);
        c.add_single_coverage(10.2);
        c.start_entry(3, "contig2");
        c.add_single_coverage(20.1);
        c.add_single_coverage(20.2);
        match c {
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                stoit_names,
                entry_names,
                coverages,
                current_stoit_index,
                current_entry_index,
                ..
            } => {
                assert_eq!(
                    vec!("stoit1".to_string(), "stoit2".to_string()),
                    stoit_names
                );
                assert_eq!(
                    vec!(
                        Some("contig1".to_string()),
                        None,
                        None,
                        Some("contig2".to_string())
                    ),
                    entry_names
                );
                assert_eq!(
                    vec![(0, vec![1.1, 1.2]), (3, vec![2.1, 2.2])],
                    recorded(&coverages[0], 2)
                );
                assert_eq!(
                    vec![(0, vec![10.1, 10.2]), (3, vec![20.1, 20.2])],
                    recorded(&coverages[1], 2)
                );
                assert_eq!(1, current_stoit_index.unwrap());
                assert_eq!(3, current_entry_index.unwrap());
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_cached_two_samples_mismatching() {
        let mut c = CoverageTakerType::new_cached_single_float_coverage_taker(2);
        c.start_stoit("stoit1");
        c.start_entry(0, "contig1");
        c.add_single_coverage(1.1);
        c.add_single_coverage(1.2);
        c.start_entry(3, "contig2");
        c.add_single_coverage(2.1);
        c.add_single_coverage(2.2);
        c.start_stoit("stoit2");
        c.start_entry(1, "contig1.5");
        c.add_single_coverage(10.1);
        c.add_single_coverage(10.2);
        c.start_entry(3, "contig2");
        c.add_single_coverage(20.1);
        c.add_single_coverage(20.2);
        c.start_entry(5, "contig5");
        c.add_single_coverage(20.1);
        c.add_single_coverage(20.2);
        match c {
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                stoit_names,
                entry_names,
                coverages,
                current_stoit_index,
                current_entry_index,
                ..
            } => {
                assert_eq!(
                    vec!("stoit1".to_string(), "stoit2".to_string()),
                    stoit_names
                );
                assert_eq!(
                    vec!(
                        Some("contig1".to_string()),
                        Some("contig1.5".to_string()),
                        None,
                        Some("contig2".to_string()),
                        None,
                        Some("contig5".to_string())
                    ),
                    entry_names
                );
                assert_eq!(
                    vec![(0, vec![1.1, 1.2]), (3, vec![2.1, 2.2])],
                    recorded(&coverages[0], 2)
                );
                assert_eq!(
                    vec![
                        (1, vec![10.1, 10.2]),
                        (3, vec![20.1, 20.2]),
                        (5, vec![20.1, 20.2])
                    ],
                    recorded(&coverages[1], 2)
                );
                assert_eq!(1, current_stoit_index.unwrap());
                assert_eq!(5, current_entry_index.unwrap());
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_cached_cursor_fills_zeros() {
        let mut c = CoverageTakerType::new_cached_single_float_coverage_taker(2);
        c.start_stoit("stoit1");
        c.start_entry(0, "contig1");
        c.add_single_coverage(1.1);
        c.add_single_coverage(1.2);
        c.start_entry(3, "contig2");
        c.add_single_coverage(2.1);
        c.add_single_coverage(2.2);

        c.start_stoit("stoit2");
        c.start_entry(1, "contig1.5");
        c.add_single_coverage(10.1);
        c.add_single_coverage(10.2);
        c.start_entry(3, "contig2");
        c.add_single_coverage(20.1);
        c.add_single_coverage(20.2);
        c.start_entry(5, "contig5");
        c.add_single_coverage(20.1);
        c.add_single_coverage(20.2);

        match c {
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                entry_names,
                coverages,
                ..
            } => {
                let rows: Vec<usize> = observed_entry_indices(&entry_names).collect();
                assert_eq!(vec![0, 1, 3, 5], rows);

                let mut cursor = StoitCoverageCursor::new(&coverages[0], 2);
                let got: Vec<Option<&[f32]>> =
                    rows.iter().map(|e| cursor.coverages_for(*e)).collect();
                assert_eq!(
                    vec![Some(&[1.1, 1.2][..]), None, Some(&[2.1, 2.2][..]), None],
                    got
                );

                let mut cursor = StoitCoverageCursor::new(&coverages[1], 2);
                let got: Vec<Option<&[f32]>> =
                    rows.iter().map(|e| cursor.coverages_for(*e)).collect();
                assert_eq!(
                    vec![
                        None,
                        Some(&[10.1, 10.2][..]),
                        Some(&[20.1, 20.2][..]),
                        Some(&[20.1, 20.2][..])
                    ],
                    got
                );
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_cached_cursor_one_coverage() {
        let mut c = CoverageTakerType::new_cached_single_float_coverage_taker(1);
        c.start_stoit("stoit1");
        c.start_entry(0, "contig1");
        c.add_single_coverage(1.1);
        c.start_entry(3, "contig2");
        c.add_single_coverage(2.1);

        c.start_stoit("stoit2");
        c.start_entry(1, "contig1.5");
        c.add_single_coverage(10.1);
        c.start_entry(3, "contig2");
        c.add_single_coverage(20.1);
        c.start_entry(5, "contig5");
        c.add_single_coverage(20.1);

        match c {
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                entry_names,
                coverages,
                ..
            } => {
                let rows: Vec<usize> = observed_entry_indices(&entry_names).collect();
                let mut cursor = StoitCoverageCursor::new(&coverages[0], 1);
                let got: Vec<f32> = rows
                    .iter()
                    .map(|e| cursor.coverages_for(*e).map_or(0.0, |c| c[0]))
                    .collect();
                assert_eq!(vec![1.1, 0.0, 2.1, 0.0], got);
                let mut cursor = StoitCoverageCursor::new(&coverages[1], 1);
                let got: Vec<f32> = rows
                    .iter()
                    .map(|e| cursor.coverages_for(*e).map_or(0.0, |c| c[0]))
                    .collect();
                assert_eq!(vec![0.0, 10.1, 20.1, 20.1], got);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_stoit_coverages_becomes_sparse() {
        let mut s = StoitCoverages::default();
        s.push_entry(0);
        s.push_entry(1);
        assert_eq!(None, s.entry_indices);
        s.push_entry(4);
        assert_eq!(Some(vec![0, 1, 4]), s.entry_indices);
        s.push_entry(5);
        assert_eq!(
            vec![0, 1, 4, 5],
            (0..4).map(|i| s.entry_index(i)).collect::<Vec<_>>()
        );
    }

    #[test]
    #[should_panic(expected = "increasing order")]
    fn test_stoit_coverages_rejects_out_of_order_entry() {
        let mut s = StoitCoverages::default();
        s.push_entry(0);
        s.push_entry(4);
        s.push_entry(2);
    }

    #[test]
    #[should_panic(expected = "increasing order")]
    fn test_stoit_coverages_rejects_repeated_entry() {
        let mut s = StoitCoverages::default();
        s.push_entry(0);
        s.push_entry(0);
    }
}
