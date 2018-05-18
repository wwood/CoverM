use std;




pub trait MosdepthGenomeCoverageEstimator {
    fn setup(&mut self);

    fn add_contig(&mut self, ups_and_downs: &Vec<i32>);

    fn calculate_coverage(&mut self, unobserved_contig_length: u32) -> f32;

    fn print_genome<'a >(&self, stoit_name: &str, genome: &str, coverage: &f32,
                         print_stream: &'a mut std::io::Write) -> &'a mut std::io::Write {
        writeln!(print_stream, "{}\t{}\t{}",
                 stoit_name,
                 genome,
                 coverage).unwrap();
        return print_stream;
    }

    fn print_zero_coverage<'a>(&self, stoit_name: &str, genome: &str,
                               print_stream: &'a mut std::io::Write) -> &'a mut std::io::Write {
        writeln!(print_stream, "{}\t{}\t0.0",
                 stoit_name,
                 genome).unwrap();
        return print_stream;
    }
}

pub struct MeanGenomeCoverageEstimator {
    total_count: u32,
    total_bases: u32,
    num_covered_bases: u32,
    min_fraction_covered_bases: f32
}
impl MeanGenomeCoverageEstimator {
    pub fn new(min_fraction_covered_bases: f32) -> MeanGenomeCoverageEstimator {
        MeanGenomeCoverageEstimator {
            total_count: 0,
            total_bases: 0,
            num_covered_bases: 0,
            min_fraction_covered_bases: min_fraction_covered_bases
        }
    }
}
impl MosdepthGenomeCoverageEstimator for MeanGenomeCoverageEstimator {
    fn setup(&mut self) {
        self.total_count = 0;
        self.total_bases = 0;
    }

    fn add_contig(&mut self, ups_and_downs: &Vec<i32>) {
        let len = ups_and_downs.len();
        self.total_bases += len as u32;
        let mut cumulative_sum: i32 = 0;
        for i in 0..len {
            let current = ups_and_downs[i as usize];
            cumulative_sum += current;
            if cumulative_sum > 0 {
                self.num_covered_bases += 1
            }
            self.total_count += cumulative_sum as u32;
        }
    }

    fn calculate_coverage(&mut self, unobserved_contig_length: u32) -> f32 {
        let final_total_bases = self.total_bases + unobserved_contig_length;
        if final_total_bases == 0 ||
            (self.num_covered_bases as f32 / final_total_bases as f32) < self.min_fraction_covered_bases {
            return 0.0
        } else {
            return self.total_count as f32 / final_total_bases as f32
        }
    }
}

pub struct TrimmedMeanGenomeCoverageEstimator {
    counts: Vec<u32>,
    observed_contig_length: u32,
    num_covered_bases: u32,
    min: f32,
    max: f32,
    min_fraction_covered_bases: f32
}
impl TrimmedMeanGenomeCoverageEstimator {
    pub fn new(min: f32, max: f32, min_fraction_covered_bases: f32) -> TrimmedMeanGenomeCoverageEstimator {
        TrimmedMeanGenomeCoverageEstimator {
            counts: vec!(),
            observed_contig_length: 0,
            num_covered_bases: 0,
            min_fraction_covered_bases: min_fraction_covered_bases,
            min: min,
            max: max
        }
    }
}
impl MosdepthGenomeCoverageEstimator for TrimmedMeanGenomeCoverageEstimator {
    fn setup(&mut self) {
        self.observed_contig_length = 0;
        self.num_covered_bases = 0;
        self.counts = vec!();
    }

    fn add_contig(&mut self, ups_and_downs: &Vec<i32>) {
        let len1 = ups_and_downs.len();
        debug!("Adding len1 {}", len1);
        self.observed_contig_length += len1 as u32;
        let mut cumulative_sum: i32 = 0;
        for current in ups_and_downs {
            if *current != 0 {
                debug!("cumulative sum {} and current {}", cumulative_sum, current);
                debug!("At i some, ups and downs {:?}", ups_and_downs);
            }
            cumulative_sum += current;
            if cumulative_sum > 0 {
                self.num_covered_bases += 1
            }
            if self.counts.len() <= cumulative_sum as usize {
                self.counts.resize(cumulative_sum as usize +1, 0);
            }
            self.counts[cumulative_sum as usize] += 1
        }
    }

    fn calculate_coverage(&mut self, unobserved_contig_length: u32) -> f32 {
        let total_bases = self.observed_contig_length + unobserved_contig_length;
        debug!("Calculating coverage with observed length {}, unobserved_length {} and counts {:?}", self.num_covered_bases, unobserved_contig_length, self.counts);
        let answer = match total_bases {
            0 => 0.0,
            _ => {
                if (self.num_covered_bases as f32 / total_bases as f32) < self.min_fraction_covered_bases {
                    0.0
                } else {
                    let min_index: usize = (self.min * total_bases as f32).floor() as usize;
                    let max_index: usize = (self.max * total_bases as f32).ceil() as usize;
                    if self.num_covered_bases == 0 {return 0.0;}
                    self.counts[0] += unobserved_contig_length;

                    let mut num_accounted_for: usize = 0;
                    let mut total: usize = 0;
                    let mut started = false;
                    let mut i = 0;
                    for num_covered in self.counts.iter() {
                        num_accounted_for += *num_covered as usize;
                        debug!("start: i {}, num_accounted_for {}, total {}, min {}, max {}", i, num_accounted_for, total, min_index, max_index);
                        if num_accounted_for >= min_index {
                            debug!("inside");
                            if started {
                                if num_accounted_for > max_index {
                                    let num_wanted = max_index - (num_accounted_for - *num_covered as usize) + 1;
                                    debug!("num wanted1: {}", num_wanted);
                                    total += num_wanted * i;
                                    break;
                                } else {
                                    total += *num_covered as usize * i;
                                }
                            } else {
                                if num_accounted_for > max_index {
                                    // all coverages are the same in the trimmed set
                                    total = (max_index-min_index+1) * i;
                                    started = true
                                } else if num_accounted_for < min_index {
                                    debug!("too few on first")
                                } else {
                                    let num_wanted = num_accounted_for - min_index + 1;
                                    debug!("num wanted2: {}", num_wanted);
                                    total = num_wanted * i;
                                    started = true;
                                }
                            }
                        }
                        debug!("end i {}, num_accounted_for {}, total {}", i, num_accounted_for, total);

                        i += 1;
                    }
                    total as f32 / (max_index-min_index) as f32
                }
            }
        };
        return answer
    }
}






pub struct PileupCountsGenomeCoverageEstimator {
    counts: Vec<u32>,
    observed_contig_length: u32,
    num_covered_bases: u32,
    min_fraction_covered_bases: f32
}
impl PileupCountsGenomeCoverageEstimator {
    pub fn new(min_fraction_covered_bases: f32) -> PileupCountsGenomeCoverageEstimator {
        PileupCountsGenomeCoverageEstimator {
            counts: vec!(),
            observed_contig_length: 0,
            num_covered_bases: 0,
            min_fraction_covered_bases: min_fraction_covered_bases
        }
    }
}

impl MosdepthGenomeCoverageEstimator for PileupCountsGenomeCoverageEstimator {
    fn setup(&mut self) {
        self.observed_contig_length = 0;
        self.num_covered_bases = 0;
        self.counts = vec!();
    }

    // Method directly copied from Trimmed mean estimator.
    fn add_contig(&mut self, ups_and_downs: &Vec<i32>) {
        let len1 = ups_and_downs.len();
        debug!("Adding len1 {}", len1);
        self.observed_contig_length += len1 as u32;
        let mut cumulative_sum: i32 = 0;
        for current in ups_and_downs {
            if *current != 0 {
                debug!("cumulative sum {} and current {}", cumulative_sum, current);
                debug!("At i some, ups and downs {:?}", ups_and_downs);
            }
            cumulative_sum += current;
            if cumulative_sum > 0 {
                self.num_covered_bases += 1
            }
            if self.counts.len() <= cumulative_sum as usize {
                self.counts.resize(cumulative_sum as usize +1, 0);
            }
            self.counts[cumulative_sum as usize] += 1
        }
    }

    fn calculate_coverage(&mut self, unobserved_contig_length: u32) -> f32 {
        // No need to actually calculate any kind of coverage, just return
        // whether any coverage was detected
        match self.observed_contig_length {
            0 => 0.0,
            _ => {
                let total_bases = self.observed_contig_length + unobserved_contig_length;
                if (self.num_covered_bases as f32 / total_bases as f32) < self.min_fraction_covered_bases {
                    0.0
                } else {
                    // Hack: Return the number of zero coverage bases as the
                    // coverage, plus 1 so it is definitely non-zero, so that
                    // the print_genome function knows this info.
                    (total_bases - self.num_covered_bases + 1) as f32
                }
            }
        }
    }

    fn print_genome<'a >(&self, stoit_name: &str, genome: &str, coverage: &f32,
                         print_stream: &'a mut std::io::Write) -> &'a mut std::io::Write {
        let mut i = 0;
        debug!("starting to print {}", genome);
        debug!("{:?}",self.counts);
        for num_covered in self.counts.iter() {
            let cov: u32 = match i {
                0 => coverage.floor() as u32 - 1,
                _ => *num_covered
            };
            writeln!(print_stream, "{}\t{}\t{:}\t{:}", stoit_name, genome, i, cov).unwrap();
            i += 1
        }
        return print_stream;
    }

    fn print_zero_coverage<'a>(&self, _stoit_name: &str, _genome: &str,
                               print_stream: &'a mut std::io::Write) -> &'a mut std::io::Write {
        // zeros are not printed usually, so do not print the length.
        return print_stream;
    }
}
