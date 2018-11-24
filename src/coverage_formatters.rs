use std;

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
}

pub trait CoverageTaker {
    fn start_stoit(&mut self, stoit_name: &str);
    fn start_entry(&mut self, entry_name: &str);
    fn add_single_coverage(&mut self, coverage: f32);
    // This function is only used with
    // CoverageEstimator::PileupCountsGenomeCoverageEstimator, and is a bit of a
    // hack to have it in this generalised trait, really.
    fn add_coverage_entry(&mut self, num_reads: usize, num_bases: u32);
    fn finish_entry(&mut self);
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
        }
    }

    fn start_entry(&mut self, entry_name: &str) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter{
                print_stream,
                current_stoit
            } => {
                match current_stoit {
                    Some(stoit) => {
                        write!(
                            print_stream, "{}\t{}", stoit, entry_name).unwrap()},
                    None => {panic!("programming error")}
                }
            },
            CoverageTakerType::PileupCoverageCoveragePrinter{
                print_stream: _,
                current_stoit: _,
                ref mut current_entry} => {
                *current_entry = Some(entry_name.to_owned());
            },
        }
    }

    fn add_single_coverage(&mut self, coverage: f32) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter{
                print_stream, .. } => {
                if coverage == 0.0 {
                    write!(print_stream, "\t0.0").unwrap()
                } else {
                    write!(print_stream, "\t{}", coverage).unwrap()
                }
            },
            CoverageTakerType::PileupCoverageCoveragePrinter{..} => {
                panic!("programming error");
            }
        }
    }

    // This function is only used with
    // CoverageEstimator::PileupCountsGenomeCoverageEstimator, and is a bit of a
    // hack, really.
    fn add_coverage_entry(&mut self, num_reads: usize, num_bases: u32) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter{..} => {
                panic!("Programming error")
            },
            CoverageTakerType::PileupCoverageCoveragePrinter{
                ref mut print_stream,
                ref mut current_stoit,
                ref mut current_entry
            } => {
                let stoit = match current_stoit {
                    Some(ref stoit) => stoit,
                    None => panic!("programming error")
                };
                let entry = match current_entry {
                    Some(ref entry) => entry,
                    None => panic!("programming error")
                };
                writeln!(
                    print_stream,
                    "{}\t{}\t{}\t{}",
                    stoit,
                    entry,
                    num_reads,
                    num_bases);
            }
        }
    }

    fn finish_entry(&mut self) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter{
                print_stream, .. } => {
                writeln!(print_stream).unwrap()
            },
            CoverageTakerType::PileupCoverageCoveragePrinter{..} => {}
        }
    }
}
