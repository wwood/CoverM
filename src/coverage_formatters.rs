use std;

pub enum CoverageTakerType<'a> {
    SingleFloatCoverageStreamingCoveragePrinter {
        print_stream: &'a mut std::io::Write
    },
    PileupCoverageCoveragePrinter {
        print_stream: &'a mut std::io::Write,
        current_stoit: Option<String>,
        current_entry: Option<String>,
    }
}

pub trait CoverageTaker {
    fn start_entry(&mut self, entry_name: &str, stoit_name: &str);
    fn add_single_coverage(&mut self, coverage: f32);
    // This function is only used with
    // CoverageEstimator::PileupCountsGenomeCoverageEstimator, and is a bit of a
    // hack to have it in this generalised trait, really.
    fn add_coverage_entry(&mut self, num_reads: usize, num_bases: u32);
    fn finish_entry(&mut self);
}

impl<'a> CoverageTaker for CoverageTakerType<'a> {
    fn start_entry(&mut self, entry_name: &str, stoit_name: &str) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter{print_stream} => {
                write!(print_stream, "{}\t{}", stoit_name, entry_name).unwrap()
            },
            CoverageTakerType::PileupCoverageCoveragePrinter{
                print_stream: _,
                ref mut current_stoit,
                ref mut current_entry} => {
                *current_entry = Some(entry_name.to_owned());
                *current_stoit = Some(stoit_name.to_owned());
            }
        }
    }

    fn add_single_coverage(&mut self, coverage: f32) {
        match self {
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter{print_stream} => {
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
            CoverageTakerType::SingleFloatCoverageStreamingCoveragePrinter{print_stream} => {
                writeln!(print_stream).unwrap()
            },
            CoverageTakerType::PileupCoverageCoveragePrinter{..} => {}
        }
    }
}
