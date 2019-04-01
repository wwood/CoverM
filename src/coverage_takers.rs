use std;
use std::fmt;

pub enum CoverageTakerType<'a> {
    SingleFloatCoverageStreamingCoveragePrinter {
        print_stream: &'a mut std::io::Write,
        current_stoit: Option<String>
    },
    PileupCoverageCoveragePrinter {
        print_stream: &'a mut std::io::Write,
        current_stoit: Option<String>,
        current_entry: Option<String>,
    },
    CachedSingleFloatCoverageTaker {
        stoit_names: Vec<String>,
        entry_names: Vec<Option<String>>,
        // list of per-stoit lists of recorded coverages
        coverages: Vec<Vec<CoverageEntry>>,
        current_stoit_index: Option<usize>,
        current_entry_index: Option<usize>,
        num_coverages: usize, // number of different coverage calculations
    }
}



pub trait CoverageTaker {
    fn start_stoit(&mut self, stoit_name: &str);
    fn start_entry(&mut self, entry_order_id: usize, entry_name: &str);
    fn add_single_coverage(&mut self, coverage: f32);
    // This function is only used with
    // CoverageEstimator::PileupCountsGenomeCoverageEstimator, and is a bit of a
    // hack to have it in this generalised trait, really.
    fn add_coverage_entry(&mut self, num_reads: usize, num_bases: u32);
    fn finish_entry(&mut self);
}

#[derive(PartialEq, Debug)]
pub struct CoverageEntry {
    pub entry_index: usize,
    pub coverage: f32
}

impl<'a> CoverageTakerType<'a> {
    pub fn new_single_float_coverage_streaming_coverage_printer(
        print_stream: &mut std::io::Write) -> CoverageTakerType {
        CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter {
            print_stream: print_stream,
            current_stoit: None
        }
    }
    pub fn new_pileup_coverage_coverage_printer(
        print_stream: &mut std::io::Write) -> CoverageTakerType {
        CoverageTakerType::PileupCoverageCoveragePrinter {
            print_stream: print_stream,
            current_stoit: None,
            current_entry: None
        }
    }
    pub fn new_cached_single_float_coverage_taker(num_coverages: usize) -> CoverageTakerType<'a> {
        CoverageTakerType::CachedSingleFloatCoverageTaker {
            stoit_names: vec!(),
            entry_names: vec!(),
            coverages: vec!(),
            current_stoit_index: None,
            current_entry_index: None,
            num_coverages: num_coverages,
        }
    }
}



impl<'a> CoverageTaker for CoverageTakerType<'a> {
    fn start_stoit(&mut self, stoit_name: &str) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter {
                print_stream:_,
                ref mut current_stoit,
            } => {
                *current_stoit = Some(stoit_name.to_owned())
            },
            CoverageTakerType::PileupCoverageCoveragePrinter {
                print_stream: _,
                ref mut current_stoit,
                current_entry: _
            } => {
                *current_stoit = Some(stoit_name.to_owned())
            }
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                ref mut stoit_names,
                entry_names: _,
                ref mut coverages,
                ref mut current_stoit_index,
                current_entry_index: _, ..
            } => {
                stoit_names.push(stoit_name.to_owned());
                coverages.push(vec!());
                *current_stoit_index = Some(stoit_names.len() - 1);
            }
        }
    }

    fn start_entry(&mut self, entry_order_id: usize, entry_name: &str) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter{
                print_stream,
                current_stoit
            } => {
                match current_stoit {
                    Some(stoit) => {
                        write!(
                            print_stream, "{}\t{}", stoit, entry_name).unwrap()},
                    None => unreachable!()
                }
            },
            CoverageTakerType::PileupCoverageCoveragePrinter{
                print_stream: _,
                current_stoit: _,
                ref mut current_entry} => {
                *current_entry = Some(entry_name.to_owned());
            },
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                stoit_names: _,
                ref mut entry_names,
                coverages: _,
                current_stoit_index: _,
                ref mut current_entry_index, ..
            } => {
                debug!("Starting an entry with ID {} and name {}",
                      entry_order_id, entry_name);
                // if the first time this entry has been seen, record its name.
                if entry_order_id >= entry_names.len() {
                    debug!("Adding 1 to entry order id {}, where names size is {}",
                           entry_order_id, entry_names.len());
                    entry_names.resize(entry_order_id+1, None)
                }
                if entry_names[entry_order_id].is_none() {
                    entry_names[entry_order_id] = Some(entry_name.to_owned());
                }
                match &entry_names[entry_order_id] {
                    Some(prev) => {
                        if prev != entry_name {
                            panic!("Found a difference amongst the reference sets used for \
                                    mapping. For this (non-streaming) usage of CoverM, all \
                                    BAM files must have the same set of reference sequences. \
                                    Previous entry was {}, new is {}",
                            prev, entry_name)
                        }
                    },
                    None => {} // should never happen. There's probably a more
                    // idiomatic way to structure this.
                }
                *current_entry_index = Some(entry_order_id);
            }
        }
    }

    fn add_single_coverage(&mut self, coverage: f32) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter{
                print_stream, .. } => {
                if coverage == 0.0 {
                    write!(print_stream, "\t0").unwrap()
                } else {
                    write!(print_stream, "\t{}", coverage).unwrap()
                }
            },
            CoverageTakerType::PileupCoverageCoveragePrinter{..} => {
                unreachable!();
            },
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                stoit_names: _,
                entry_names: _,
                ref mut coverages,
                ref current_stoit_index,
                ref current_entry_index, ..
            } => {
                coverages[current_stoit_index.unwrap()].push(CoverageEntry {
                    entry_index: current_entry_index.unwrap(),
                    coverage: coverage
                })
            }
        }
    }

    // This function is only used with
    // CoverageEstimator::PileupCountsGenomeCoverageEstimator, and is a bit of a
    // hack, really.
    fn add_coverage_entry(&mut self, num_reads: usize, num_bases: u32) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter{..} |
            CoverageTakerType::CachedSingleFloatCoverageTaker{..} => {
                unreachable!()
            },
            CoverageTakerType::PileupCoverageCoveragePrinter{
                ref mut print_stream,
                ref mut current_stoit,
                ref mut current_entry
            } => {
                let stoit = match current_stoit {
                    Some(ref stoit) => stoit,
                    None => unreachable!()
                };
                let entry = match current_entry {
                    Some(ref entry) => entry,
                    None => unreachable!()
                };
                writeln!(
                    print_stream,
                    "{}\t{}\t{}\t{}",
                    stoit,
                    entry,
                    num_reads,
                    num_bases).unwrap();
            }
        }
    }

    fn finish_entry(&mut self) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter{
                print_stream, .. } => {
                writeln!(print_stream).unwrap()
            },
            CoverageTakerType::PileupCoverageCoveragePrinter{..} |
            CoverageTakerType::CachedSingleFloatCoverageTaker{..} => {}
        }
    }
}


#[derive(PartialEq, Debug)]
pub struct EntryAndCoverages {
    pub entry_index: usize,
    pub stoit_index: usize,
    pub coverages: Vec<f32>
}


pub struct CoverageTakerTypeIterator<'a> {
    coverage_taker_type: &'a CoverageTakerType<'a>,
    // indices for iterating
    iter_next_entry_indices: Vec<usize>, // indexes into coverages[stoit]
    iter_current_stoit_index: usize, // indexes into coverages
    iter_last_entry_order_index: Option<usize>, // index into entry_names
}


impl<'a> std::fmt::Debug for CoverageTakerTypeIterator<'a> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.debug_struct("CoverageTakerTypeIterator")
            .field("iter_next_entry_indices", &self.iter_next_entry_indices)
            .field("iter_current_stoit_index", &self.iter_current_stoit_index)
            .field("iter_last_entry_order_index", &self.iter_last_entry_order_index)
            .finish()
    }
}


impl<'a> CoverageTakerType<'a> {
    pub fn generate_iterator(&'a self) -> CoverageTakerTypeIterator<'a> {
        match self {
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                stoit_names, ..} => {
                return CoverageTakerTypeIterator {
                    coverage_taker_type: &self,
                    iter_next_entry_indices: vec![0; stoit_names.len()], // indexes into coverages[stoit]
                    iter_current_stoit_index: 0, // indexes into coverages
                    iter_last_entry_order_index: None, // index into entry_names
                }
            },
            _ => unreachable!()
        }
    }
}

impl<'a> Iterator for CoverageTakerTypeIterator<'a> {
    type Item = EntryAndCoverages;

    fn next(&mut self) -> Option<EntryAndCoverages> {
        match self.coverage_taker_type {
            CoverageTakerType::CachedSingleFloatCoverageTaker {
                ref stoit_names,
                entry_names: _,
                ref coverages,
                current_stoit_index: _,
                current_entry_index: _,
                ref num_coverages,
            } => {
                while self.iter_current_stoit_index <= stoit_names.len() {
                    let mut lowest_entry_order = None;
                    // Search for the lowest entry_index amongst the head of the
                    // queues of entries for each stoit
                    for (stoit_i, coverage_i) in (self.iter_next_entry_indices)
                        .iter().enumerate() {
                            // if there are not more entries in this stoit, do nothing
                            if *coverage_i < coverages[stoit_i].len() {
                                // if current lowest is None or higher than current, update
                                let current_coverage_entry = &coverages[stoit_i][*coverage_i];
                                if self.iter_last_entry_order_index.is_none() ||
                                    current_coverage_entry.entry_index >
                                    self.iter_last_entry_order_index.unwrap() {
                                    match lowest_entry_order {
                                        None => {
                                            lowest_entry_order = Some(current_coverage_entry.entry_index);
                                        },
                                        Some(old_best_entry_order) => {
                                            if current_coverage_entry.entry_index < old_best_entry_order {
                                                // new winner
                                                lowest_entry_order = Some(current_coverage_entry.entry_index);
                                            }
                                            // else must be the same or a loser, do nothing
                                        }
                                    }
                                }
                            }
                        }
                    // Winning entry picked, update state and return it
                    match lowest_entry_order {
                        Some(lowest_entry_i) => {
                            let mut chosen_stoit_entry_id =
                                self.iter_next_entry_indices[self.iter_current_stoit_index];
                            let potential_next_stoit_list = &coverages[self.iter_current_stoit_index];
                            let mut struct_to_return;
                            if chosen_stoit_entry_id >= potential_next_stoit_list.len() ||
                                potential_next_stoit_list[chosen_stoit_entry_id].entry_index != lowest_entry_i {
                                    // There are no more coverages from this stoit,
                                    // now we are just returning zeroes to fill out
                                    // the larger matrix.
                                    struct_to_return = Some(EntryAndCoverages{
                                        entry_index: lowest_entry_i,
                                        stoit_index: self.iter_current_stoit_index,
                                        coverages: vec![0.0; *num_coverages]
                                    })
                                } else {
                                    let mut to_return = vec!();
                                    // collect the coverages to return from the
                                    // stoit currently being iterated.
                                    for _ in 0..*num_coverages {
                                        to_return.push(potential_next_stoit_list[chosen_stoit_entry_id].coverage);
                                        chosen_stoit_entry_id += 1;
                                    }

                                    struct_to_return = Some(EntryAndCoverages{
                                        entry_index: lowest_entry_i,
                                        stoit_index: self.iter_current_stoit_index,
                                        coverages: to_return
                                    })
                                }

                            // update the pointers for each stoit
                            for stoit_i in 0..stoit_names.len() {
                                if coverages[stoit_i].len() > self.iter_next_entry_indices[stoit_i] &&
                                    coverages[stoit_i][self.iter_next_entry_indices[stoit_i]].entry_index
                                    == lowest_entry_i {
                                        self.iter_next_entry_indices[stoit_i] += *num_coverages
                                    }
                            }
                            self.iter_last_entry_order_index = Some(lowest_entry_i);
                            return struct_to_return;
                        },
                        None => {
                            // all coverages from the current stoit have been returned
                            self.iter_current_stoit_index += 1;
                            if self.iter_current_stoit_index >= stoit_names.len() {
                                return None // Finished all iteration now.
                            }
                            self.iter_next_entry_indices = vec![0; stoit_names.len()];
                            self.iter_last_entry_order_index = None;
                        }
                    }
                }
                return None
            },
            _ => unreachable!()
        }
    }
}




#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cached_hello_world() {
        let mut c = CoverageTakerType::new_cached_single_float_coverage_taker(2);
        c.start_stoit("stoit1");
        c.start_entry(0, "contig1");
        c.add_single_coverage(1.1);
        c.add_single_coverage(1.2);
        match c {
            CoverageTakerType::CachedSingleFloatCoverageTaker{
                stoit_names,
                entry_names,
                coverages,
                current_stoit_index,
                current_entry_index,
                num_coverages,
            } => {
                assert_eq!(vec!("stoit1".to_string()), stoit_names);
                assert_eq!(vec!(Some("contig1".to_string())), entry_names);
                assert_eq!(vec![vec![
                    CoverageEntry { entry_index: 0, coverage: 1.1},
                    CoverageEntry { entry_index: 0, coverage: 1.2}]],
                           coverages);
                assert_eq!(0, current_stoit_index.unwrap());
                assert_eq!(0, current_entry_index.unwrap());
                assert_eq!(2, num_coverages);
            },
            _ => panic!()
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
            CoverageTakerType::CachedSingleFloatCoverageTaker{
                stoit_names,
                entry_names,
                coverages,
                current_stoit_index,
                current_entry_index,
                ..
            } => {
                assert_eq!(vec!("stoit1".to_string(), "stoit2".to_string()), stoit_names);
                assert_eq!(vec!(
                    Some("contig1".to_string()), None, None, Some("contig2".to_string())),
                    entry_names);
                assert_eq!(vec![vec![
                    CoverageEntry { entry_index: 0, coverage: 1.1},
                    CoverageEntry { entry_index: 0, coverage: 1.2},
                    CoverageEntry { entry_index: 3, coverage: 2.1},
                    CoverageEntry { entry_index: 3, coverage: 2.2},
                ], vec![
                    CoverageEntry { entry_index: 0, coverage: 10.1},
                    CoverageEntry { entry_index: 0, coverage: 10.2},
                    CoverageEntry { entry_index: 3, coverage: 20.1},
                    CoverageEntry { entry_index: 3, coverage: 20.2},
                ]],
                           coverages);
                assert_eq!(1, current_stoit_index.unwrap());
                assert_eq!(3, current_entry_index.unwrap());
            },
            _ => panic!()
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
            CoverageTakerType::CachedSingleFloatCoverageTaker{
                stoit_names,
                entry_names,
                coverages,
                current_stoit_index,
                current_entry_index,
                ..
            } => {
                assert_eq!(vec!("stoit1".to_string(), "stoit2".to_string()), stoit_names);
                assert_eq!(vec!(
                    Some("contig1".to_string()),
                    Some("contig1.5".to_string()),
                    None,
                    Some("contig2".to_string()),
                    None,
                    Some("contig5".to_string())),
                    entry_names);
                assert_eq!(vec![vec![
                    CoverageEntry { entry_index: 0, coverage: 1.1},
                    CoverageEntry { entry_index: 0, coverage: 1.2},
                    CoverageEntry { entry_index: 3, coverage: 2.1},
                    CoverageEntry { entry_index: 3, coverage: 2.2},
                ], vec![
                    CoverageEntry { entry_index: 1, coverage: 10.1},
                    CoverageEntry { entry_index: 1, coverage: 10.2},
                    CoverageEntry { entry_index: 3, coverage: 20.1},
                    CoverageEntry { entry_index: 3, coverage: 20.2},
                    CoverageEntry { entry_index: 5, coverage: 20.1},
                    CoverageEntry { entry_index: 5, coverage: 20.2},
                ]],
                           coverages);
                assert_eq!(1, current_stoit_index.unwrap());
                assert_eq!(5, current_entry_index.unwrap());
            },
            _ => panic!()
        }
    }


    #[test]
    fn test_cached_next() {
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

        let mut it = c.generate_iterator();

        assert_eq!(Some(EntryAndCoverages{
            entry_index: 0,
            stoit_index: 0,
            coverages: vec![1.1,1.2]
        }), it.next());
        assert_eq!(Some(EntryAndCoverages{
            entry_index: 1,
            stoit_index: 0,
            coverages: vec![0.0,0.0]
        }), it.next());
        assert_eq!(Some(EntryAndCoverages{
            entry_index: 3,
            stoit_index: 0,
            coverages: vec![2.1,2.2]
        }), it.next());
        assert_eq!(Some(EntryAndCoverages{
            entry_index: 5,
            stoit_index: 0,
            coverages: vec![0.0,0.0]
        }), it.next());

        assert_eq!(Some(EntryAndCoverages{
            entry_index: 0,
            stoit_index: 1,
            coverages: vec![0.0,0.0]
        }), it.next());
        assert_eq!(Some(EntryAndCoverages{
            entry_index: 1,
            stoit_index: 1,
            coverages: vec![10.1,10.2]
        }), it.next());
        assert_eq!(Some(EntryAndCoverages{
            entry_index: 3,
            stoit_index: 1,
            coverages: vec![20.1,20.2]
        }), it.next());
        assert_eq!(Some(EntryAndCoverages{
            entry_index: 5,
            stoit_index: 1,
            coverages: vec![20.1,20.2]
        }), it.next());
        assert_eq!(None, it.next());
    }

    #[test]
    fn test_cached_next_one_coverage() {
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

        let mut it = c.generate_iterator();

        println!("it: {:?}", it);
        assert_eq!(Some(EntryAndCoverages{
            entry_index: 0,
            stoit_index: 0,
            coverages: vec![1.1]
        }), it.next());
        println!("it: {:?}", it);
        assert_eq!(Some(EntryAndCoverages{
            entry_index: 1,
            stoit_index: 0,
            coverages: vec![0.0]
        }), it.next());
        assert_eq!(Some(EntryAndCoverages{
            entry_index: 3,
            stoit_index: 0,
            coverages: vec![2.1]
        }), it.next());
        assert_eq!(Some(EntryAndCoverages{
            entry_index: 5,
            stoit_index: 0,
            coverages: vec![0.0]
        }), it.next());

        assert_eq!(Some(EntryAndCoverages{
            entry_index: 0,
            stoit_index: 1,
            coverages: vec![0.0]
        }), it.next());
        assert_eq!(Some(EntryAndCoverages{
            entry_index: 1,
            stoit_index: 1,
            coverages: vec![10.1]
        }), it.next());
        assert_eq!(Some(EntryAndCoverages{
            entry_index: 3,
            stoit_index: 1,
            coverages: vec![20.1]
        }), it.next());
        assert_eq!(Some(EntryAndCoverages{
            entry_index: 5,
            stoit_index: 1,
            coverages: vec![20.1]
        }), it.next());
        assert_eq!(None, it.next());
    }
}


