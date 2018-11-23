use std;

pub trait CoverageTaker {
    fn start_entry(&mut self, entry_name: &str, stoit_name: &str);
    fn add_single_coverage(&mut self, coverage: f32);
    // This function is only used with
    // CoverageEstimator::PileupCountsGenomeCoverageEstimator, and is a bit of a
    // hack to have it in this generalised trait, really.
    fn add_coverage_entry(&mut self, num_reads: usize, num_bases: u32);
    fn finish_entry(&mut self);
}

pub struct SingleFloatCoverageStreamingCoveragePrinter<'a> {
    pub print_stream: &'a mut std::io::Write,
}

impl<'a> CoverageTaker for SingleFloatCoverageStreamingCoveragePrinter<'a> {
    fn start_entry(&mut self, entry_name: &str, stoit_name: &str) {
        write!(self.print_stream, "{}\t{}", stoit_name, entry_name).unwrap()
    }

    fn add_single_coverage(&mut self, coverage: f32) {
        if coverage == 0.0 {
            write!(self.print_stream, "\t0.0").unwrap()
        } else {
            write!(self.print_stream, "\t{}", coverage).unwrap()
        }
    }

    // This function is only used with
    // CoverageEstimator::PileupCountsGenomeCoverageEstimator, and is a bit of a
    // hack, really.
    fn add_coverage_entry(&mut self, _num_reads: usize, _num_bases: u32) {
        panic!("Programming error")
    }

    fn finish_entry(&mut self) {
        writeln!(self.print_stream).unwrap();
    }
}
