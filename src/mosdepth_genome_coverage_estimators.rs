use std;

#[derive(Clone, Debug)]
pub enum CoverageEstimator {
    MeanGenomeCoverageEstimator{
        total_count: u32,
        total_bases: u32,
        num_covered_bases: u32,
        min_fraction_covered_bases: f32
    },
    TrimmedMeanGenomeCoverageEstimator {
        counts: Vec<u32>,
        observed_contig_length: u32,
        num_covered_bases: u32,
        min: f32,
        max: f32,
        min_fraction_covered_bases: f32
    },
    PileupCountsGenomeCoverageEstimator {
        counts: Vec<u32>,
        observed_contig_length: u32,
        num_covered_bases: u32,
        min_fraction_covered_bases: f32
    },
    CoverageFractionGenomeCoverageEstimator {
        total_bases: u32,
        num_covered_bases: u32,
        min_fraction_covered_bases: f32
    },
    VarianceGenomeCoverageEstimator {
        counts: Vec<u32>,
        observed_contig_length: u32,
        num_covered_bases: u32,
        min_fraction_covered_bases: f32
    }
}

impl CoverageEstimator {
    pub fn column_headers(&self) -> Vec<&str> {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator{..} => {vec!("Mean")},
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator{..} => {vec!("TrimmedMean")},
            CoverageEstimator::PileupCountsGenomeCoverageEstimator{..} => {vec!("Coverage","Bases")},
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator{..} => {vec!("CoveredFraction")},
            CoverageEstimator::VarianceGenomeCoverageEstimator{..} => {vec!("Variance")}
        }
    }
}

impl CoverageEstimator {
    pub fn new_estimator_mean(
        min_fraction_covered_bases: f32)
        -> CoverageEstimator {
        CoverageEstimator::MeanGenomeCoverageEstimator {
            total_count: 0,
            total_bases: 0,
            num_covered_bases: 0,
            min_fraction_covered_bases: min_fraction_covered_bases
        }
    }
    pub fn new_estimator_trimmed_mean(
        min: f32, max: f32, min_fraction_covered_bases: f32)
        -> CoverageEstimator {
        CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
            counts: vec!(),
            observed_contig_length: 0,
            num_covered_bases: 0,
            min_fraction_covered_bases: min_fraction_covered_bases,
            min: min,
            max: max
        }
    }
    pub fn new_estimator_pileup_counts(
        min_fraction_covered_bases: f32)
        -> CoverageEstimator {
        CoverageEstimator::PileupCountsGenomeCoverageEstimator {
            counts: vec!(),
            observed_contig_length: 0,
            num_covered_bases: 0,
            min_fraction_covered_bases: min_fraction_covered_bases
        }
    }
    pub fn new_estimator_covered_fraction(
        min_fraction_covered_bases: f32)
        -> CoverageEstimator {
        CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
            total_bases: 0,
            num_covered_bases: 0,
            min_fraction_covered_bases: min_fraction_covered_bases
        }
    }
    pub fn new_estimator_variance(
        min_fraction_covered_bases: f32)
        -> CoverageEstimator {
        CoverageEstimator::VarianceGenomeCoverageEstimator {
            counts: vec!(),
            observed_contig_length: 0,
            num_covered_bases: 0,
            min_fraction_covered_bases: min_fraction_covered_bases
        }
    }
}

pub trait MosdepthGenomeCoverageEstimator {

    fn setup(&mut self);

    fn add_contig(&mut self, ups_and_downs: &Vec<i32>);

    fn calculate_coverage(&mut self, unobserved_contig_length: u32) -> f32;

    fn print_coverage<'a >(&self,
                           stoit_name: &str, genome: &str, coverage: &f32,
                           print_stream: &'a mut std::io::Write) -> &'a mut std::io::Write;

    fn print_zero_coverage<'a>(&self,
                               print_stream: &'a mut std::io::Write) -> &'a mut std::io::Write;

    fn copy(&self) -> CoverageEstimator;
}

impl MosdepthGenomeCoverageEstimator for CoverageEstimator {
    fn setup(&mut self) {
        debug!("Running setup..");
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                ref mut total_count,
                ref mut total_bases,
                ref mut num_covered_bases, ..
            } => {
                *total_count = 0;
                *total_bases = 0;
                *num_covered_bases = 0;
            },
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases, ..
            } | CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases, ..
            } | CoverageEstimator::VarianceGenomeCoverageEstimator {
                ref mut observed_contig_length,
                ref mut counts,
                ref mut num_covered_bases, ..
            } => {
                *counts = vec!();
                *observed_contig_length = 0;
                *num_covered_bases = 0;
            },
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases, ..
            } => {
                *total_bases = 0;
                *num_covered_bases = 0;
            }
        }
    }

    fn add_contig(&mut self, ups_and_downs: &Vec<i32>) {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                ref mut total_count,
                ref mut total_bases,
                ref mut num_covered_bases, ..
            } => {
                let len = ups_and_downs.len();
                *total_bases += len as u32;
                let mut cumulative_sum: i32 = 0;
                for i in 0..len {
                    let current = ups_and_downs[i as usize];
                    cumulative_sum += current;
                    if cumulative_sum > 0 {
                        *num_covered_bases += 1
                    }
                    *total_count += cumulative_sum as u32;
                }
                debug!("After adding contig, have total_count {}, total_bases {}, num_covered_bases {}",
                       total_count, total_bases, num_covered_bases);
            },
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases, ..
            } | CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases, ..
            } | CoverageEstimator::VarianceGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases, ..
            } => {
                let len1 = ups_and_downs.len();
                debug!("Adding len1 {}", len1);
                *observed_contig_length += len1 as u32;
                let mut cumulative_sum: i32 = 0;
                for current in ups_and_downs {
                    if *current != 0 {
                        debug!("cumulative sum {} and current {}", cumulative_sum, current);
                        debug!("At i some, ups and downs {:?}", ups_and_downs);
                    }
                    cumulative_sum += current;
                    if cumulative_sum > 0 {
                        *num_covered_bases += 1
                    }
                    if counts.len() <= cumulative_sum as usize {
                        (*counts).resize(cumulative_sum as usize +1, 0);
                    }
                    (*counts)[cumulative_sum as usize] += 1
                }
            },
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases, ..
            } => {
                let len = ups_and_downs.len();
                *total_bases += len as u32;
                let mut cumulative_sum: i32 = 0;
                for i in 0..len {
                    let current = ups_and_downs[i as usize];
                    cumulative_sum += current;
                    if cumulative_sum > 0 {
                        *num_covered_bases += 1
                    }
                }
            }
        }
    }

    fn calculate_coverage(&mut self, unobserved_contig_length: u32) -> f32 {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                total_count,
                total_bases,
                num_covered_bases,
                min_fraction_covered_bases
            } => {
                debug!("Calculating coverage with unobserved {}, total bases {}, num_covered_bases {}, total_count {}",
                       unobserved_contig_length, total_bases, num_covered_bases, total_count);
                let final_total_bases = *total_bases + unobserved_contig_length;
                if final_total_bases == 0 ||
                    (*num_covered_bases as f32 / final_total_bases as f32) < *min_fraction_covered_bases {
                    return 0.0
                } else {
                    return *total_count as f32 / final_total_bases as f32
                }
            },
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
                counts,
                observed_contig_length,
                num_covered_bases,
                min_fraction_covered_bases,
                min,
                max
            } => {
                let total_bases = *observed_contig_length + unobserved_contig_length;
                debug!("Calculating coverage with observed length {}, unobserved_length {} and counts {:?}",
                       num_covered_bases, unobserved_contig_length, counts);
                let answer = match total_bases {
                    0 => 0.0,
                    _ => {
                        if (*num_covered_bases as f32 / total_bases as f32) < *min_fraction_covered_bases {
                            0.0
                        } else {
                            let min_index: usize = (*min * total_bases as f32).floor() as usize;
                            let max_index: usize = (*max * total_bases as f32).ceil() as usize;
                            if *num_covered_bases == 0 {return 0.0;}
                            counts[0] += unobserved_contig_length;

                            let mut num_accounted_for: usize = 0;
                            let mut total: usize = 0;
                            let mut started = false;
                            let mut i = 0;
                            for num_covered in counts.iter() {
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
            },
            CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                counts: _,
                observed_contig_length,
                num_covered_bases,
                min_fraction_covered_bases
            } => {
                // No need to actually calculate any kind of coverage, just return
                // whether any coverage was detected
                match observed_contig_length {
                    0 => 0.0,
                    _ => {
                        let total_bases = *observed_contig_length + unobserved_contig_length;
                        if (*num_covered_bases as f32 / total_bases as f32) < *min_fraction_covered_bases {
                            0.0
                        } else {
                            // Hack: Return the number of zero coverage bases as the
                            // coverage, plus 1 so it is definitely non-zero, so that
                            // the print_genome function knows this info.
                            (total_bases - *num_covered_bases + 1) as f32
                        }
                    }
                }
            },
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
                total_bases,
                num_covered_bases,
                min_fraction_covered_bases
            } => {
                let final_total_bases = *total_bases + unobserved_contig_length;
                if final_total_bases == 0 ||
                    (*num_covered_bases as f32 / final_total_bases as f32) < *min_fraction_covered_bases {
                    return 0.0
                } else {
                    return *num_covered_bases as f32 / final_total_bases as f32
                }
            },
            CoverageEstimator::VarianceGenomeCoverageEstimator {
                observed_contig_length,
                ref mut counts,
                num_covered_bases,
                min_fraction_covered_bases
            } => {
                let total_bases = *observed_contig_length + unobserved_contig_length;
                debug!("Calculating coverage with observed length {}, unobserved_length {} and counts {:?}", num_covered_bases, unobserved_contig_length, counts);
                match total_bases {
                    0 => 0.0,
                    _ => {
                        if (*num_covered_bases as f32 / total_bases as f32) < *min_fraction_covered_bases {
                            0.0
                        } else if total_bases < 3 {
                            0.0
                        } else {
                            counts[0] += unobserved_contig_length;
                            // Calculate variance using the shifted method
                            let mut k = 0;
                            // Ensure K is within the range of coverages - take the
                            // lowest coverage.
                            while counts[k] == 0 {
                                k += 1;
                            }
                            let mut ex = 0;
                            let mut ex2 = 0;
                            for (x, num_covered) in counts.iter().enumerate() {
                                let nc = *num_covered as usize;
                                ex += (x-k) * nc;
                                ex2 += (x-k)*(x-k) * nc;
                            }
                            // Return sample variance not population variance since
                            // almost all MAGs are incomplete.
                            (ex2 as f32 - (ex*ex) as f32/total_bases as f32) / (total_bases - 1) as f32
                        }
                    }
                }
            }
        }
    }

    fn copy(&self) -> CoverageEstimator {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                total_count: _,
                total_bases: _,
                num_covered_bases: _,
                min_fraction_covered_bases
            } =>  {
                CoverageEstimator::new_estimator_mean(*min_fraction_covered_bases)
            },
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
                counts: _,
                observed_contig_length: _,
                num_covered_bases: _,
                min_fraction_covered_bases,
                min,
                max
            } => {
                CoverageEstimator::new_estimator_trimmed_mean(*min, *max, *min_fraction_covered_bases)
            },
            CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                counts: _,
                observed_contig_length: _,
                num_covered_bases: _,
                min_fraction_covered_bases
            } => {
                CoverageEstimator::new_estimator_pileup_counts(*min_fraction_covered_bases)
            },
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                min_fraction_covered_bases
            } => {
                CoverageEstimator::new_estimator_covered_fraction(*min_fraction_covered_bases)
            },
            CoverageEstimator::VarianceGenomeCoverageEstimator {
                observed_contig_length: _,
                counts: _,
                num_covered_bases: _,
                min_fraction_covered_bases
            } => {
                CoverageEstimator::new_estimator_variance(*min_fraction_covered_bases)
            }
        }
    }

    fn print_coverage<'a >(&self,
                           stoit_name: &str, genome: &str, coverage: &f32,
                           print_stream: &'a mut std::io::Write) -> &'a mut std::io::Write{
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {..} |
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator{..} |
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator{..} |
            CoverageEstimator::VarianceGenomeCoverageEstimator{..} => {
                write!(print_stream, "\t{}", coverage).unwrap();
                return print_stream;
            },
            CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                counts, ..
            } => {
                let mut i = 0;
                debug!("starting to print {}", genome);
                debug!("{:?}", counts);
                for num_covered in counts.iter() {
                    let cov: u32 = match i {
                        0 => {
                            let c = coverage.floor() as u32;
                            match c {
                                0 => 0,
                                _ => c - 1
                            }
                        },
                        _ => *num_covered
                    };
                    if i == 0 {
                        writeln!(print_stream, "\t{:}\t{:}", i, cov).unwrap();
                    } else if i+1 == counts.clone().iter().len(){
                        write!(print_stream, "{}\t{}\t{:}\t{:}", stoit_name, genome, i, cov).unwrap();
                    } else {
                        writeln!(print_stream, "{}\t{}\t{:}\t{:}", stoit_name, genome, i, cov).unwrap();
                    }
                    i += 1
                }
                return print_stream;
            }
        }
    }

    fn print_zero_coverage<'a >(&self,
                                print_stream: &'a mut std::io::Write) -> &'a mut std::io::Write {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator{..} |
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator{..} |
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator{..} |
            CoverageEstimator::VarianceGenomeCoverageEstimator{..} => {
                write!(print_stream, "\t0.0").unwrap();
                return print_stream;
            },
            CoverageEstimator::PileupCountsGenomeCoverageEstimator{..} => {
                return print_stream
            }
        }
    }
}
