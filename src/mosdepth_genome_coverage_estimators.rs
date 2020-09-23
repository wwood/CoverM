use coverage_takers::CoverageTaker;

#[derive(Clone, Debug)]
pub enum CoverageEstimator {
    MeanGenomeCoverageEstimator {
        total_count: u64,
        total_bases: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        total_mismatches: u64,
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
        exclude_mismatches: bool,
    },
    TrimmedMeanGenomeCoverageEstimator {
        counts: Vec<u64>,
        observed_contig_length: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min: f32,
        max: f32,
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
    },
    PileupCountsGenomeCoverageEstimator {
        counts: Vec<u64>,
        observed_contig_length: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
    },
    CoverageFractionGenomeCoverageEstimator {
        total_bases: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min_fraction_covered_bases: f32,
    },
    NumCoveredBasesCoverageEstimator {
        total_bases: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min_fraction_covered_bases: f32,
    },
    RPKMCoverageEstimator {
        total_bases: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min_fraction_covered_bases: f32,
    },
    TPMCoverageEstimator {
        total_bases: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min_fraction_covered_bases: f32,
    },
    VarianceGenomeCoverageEstimator {
        counts: Vec<u64>,
        observed_contig_length: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
    },
    ReferenceLengthCalculator {
        observed_contig_length: u64,
        num_mapped_reads: u64,
    },
    ReadCountCalculator {
        num_mapped_reads: u64,
    },
    ReadsPerBaseCalculator {
        observed_contig_length: u64,
        num_mapped_reads: u64,
    },
}

impl CoverageEstimator {
    pub fn column_headers(&self) -> Vec<&str> {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator { .. } => vec!["Mean"],
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator { .. } => vec!["Trimmed Mean"],
            CoverageEstimator::PileupCountsGenomeCoverageEstimator { .. } => {
                vec!["Coverage", "Bases"]
            }
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator { .. } => {
                vec!["Covered Fraction"]
            }
            CoverageEstimator::NumCoveredBasesCoverageEstimator { .. } => vec!["Covered Bases"],
            CoverageEstimator::RPKMCoverageEstimator { .. } => vec!["RPKM"],
            CoverageEstimator::TPMCoverageEstimator { .. } => vec!["TPM"],
            CoverageEstimator::VarianceGenomeCoverageEstimator { .. } => vec!["Variance"],
            CoverageEstimator::ReferenceLengthCalculator { .. } => vec!["Length"],
            CoverageEstimator::ReadCountCalculator { .. } => vec!["Read Count"],
            CoverageEstimator::ReadsPerBaseCalculator { .. } => vec!["Reads per base"],
        }
    }
}

impl CoverageEstimator {
    pub fn new_estimator_mean(
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
        exclude_mismatches: bool,
    ) -> CoverageEstimator {
        CoverageEstimator::MeanGenomeCoverageEstimator {
            total_count: 0,
            total_bases: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            total_mismatches: 0,
            min_fraction_covered_bases: min_fraction_covered_bases,
            contig_end_exclusion: contig_end_exclusion,
            exclude_mismatches: exclude_mismatches,
        }
    }
    pub fn new_estimator_trimmed_mean(
        min: f32,
        max: f32,
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
    ) -> CoverageEstimator {
        CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
            counts: vec![],
            observed_contig_length: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases: min_fraction_covered_bases,
            min: min,
            max: max,
            contig_end_exclusion: contig_end_exclusion,
        }
    }
    pub fn new_estimator_pileup_counts(
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
    ) -> CoverageEstimator {
        CoverageEstimator::PileupCountsGenomeCoverageEstimator {
            counts: vec![],
            observed_contig_length: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases: min_fraction_covered_bases,
            contig_end_exclusion: contig_end_exclusion,
        }
    }
    pub fn new_estimator_covered_fraction(min_fraction_covered_bases: f32) -> CoverageEstimator {
        CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
            total_bases: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases: min_fraction_covered_bases,
        }
    }
    pub fn new_estimator_rpkm(min_fraction_covered_bases: f32) -> CoverageEstimator {
        CoverageEstimator::RPKMCoverageEstimator {
            total_bases: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases: min_fraction_covered_bases,
        }
    }
    pub fn new_estimator_tpm(min_fraction_covered_bases: f32) -> CoverageEstimator {
        CoverageEstimator::TPMCoverageEstimator {
            total_bases: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases: min_fraction_covered_bases,
        }
    }
    pub fn new_estimator_covered_bases(min_fraction_covered_bases: f32) -> CoverageEstimator {
        CoverageEstimator::NumCoveredBasesCoverageEstimator {
            total_bases: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases: min_fraction_covered_bases,
        }
    }
    pub fn new_estimator_variance(
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
    ) -> CoverageEstimator {
        CoverageEstimator::VarianceGenomeCoverageEstimator {
            counts: vec![],
            observed_contig_length: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases: min_fraction_covered_bases,
            contig_end_exclusion: contig_end_exclusion,
        }
    }
    pub fn new_estimator_length() -> CoverageEstimator {
        CoverageEstimator::ReferenceLengthCalculator {
            observed_contig_length: 0,
            num_mapped_reads: 0,
        }
    }
    pub fn new_estimator_read_count() -> CoverageEstimator {
        CoverageEstimator::ReadCountCalculator {
            num_mapped_reads: 0,
        }
    }
    pub fn new_estimator_reads_per_base() -> CoverageEstimator {
        CoverageEstimator::ReadsPerBaseCalculator {
            observed_contig_length: 0,
            num_mapped_reads: 0,
        }
    }

    fn calculate_unobserved_bases(
        unobserved_contig_lengths: &Vec<u64>,
        contig_end_exclusion: u64,
    ) -> u64 {
        let unobserved_not_excluded = unobserved_contig_lengths
            .iter()
            .map(|l| {
                let e = &(2 * contig_end_exclusion);
                if l < e {
                    *l
                } else {
                    l - e
                }
            })
            .sum();
        unobserved_not_excluded
    }
}

pub trait MosdepthGenomeCoverageEstimator {
    fn setup(&mut self);

    fn add_contig(
        &mut self,
        ups_and_downs: &Vec<i32>,
        num_mapped_reads: u64,
        total_mismatches: u64,
    );

    fn calculate_coverage(&mut self, unobserved_contig_lengths: &Vec<u64>) -> f32;

    fn print_coverage<T: CoverageTaker>(&self, coverage: &f32, coverage_taker: &mut T);

    fn print_zero_coverage<T: CoverageTaker>(&self, coverage_taker: &mut T, entry_length: u64);

    fn copy(&self) -> CoverageEstimator;

    fn num_mapped_reads(&self) -> u64;
}

impl MosdepthGenomeCoverageEstimator for CoverageEstimator {
    fn setup(&mut self) {
        debug!("Running setup..");
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                ref mut total_count,
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ref mut total_mismatches,
                ..
            } => {
                *total_count = 0;
                *total_bases = 0;
                *num_covered_bases = 0;
                *num_mapped_reads = 0;
                *total_mismatches = 0;
            }
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::VarianceGenomeCoverageEstimator {
                ref mut observed_contig_length,
                ref mut counts,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            } => {
                *counts = vec![];
                *observed_contig_length = 0;
                *num_covered_bases = 0;
                *num_mapped_reads = 0;
            }
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::NumCoveredBasesCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::RPKMCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::TPMCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            } => {
                *total_bases = 0;
                *num_covered_bases = 0;
                *num_mapped_reads = 0;
            }
            CoverageEstimator::ReferenceLengthCalculator {
                ref mut observed_contig_length,
                ref mut num_mapped_reads,
            }
            | CoverageEstimator::ReadsPerBaseCalculator {
                ref mut observed_contig_length,
                ref mut num_mapped_reads,
            } => {
                *observed_contig_length = 0;
                *num_mapped_reads = 0;
            }
            CoverageEstimator::ReadCountCalculator {
                ref mut num_mapped_reads,
            } => {
                *num_mapped_reads = 0;
            }
        }
    }

    fn add_contig(
        &mut self,
        ups_and_downs: &Vec<i32>,
        num_mapped_reads_in_contig: u64,
        total_mismatches_in_contig: u64,
    ) {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                ref mut total_count,
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ref mut total_mismatches,
                contig_end_exclusion,
                ..
            } => {
                *num_mapped_reads += num_mapped_reads_in_contig;
                *total_mismatches += total_mismatches_in_contig;
                let len = ups_and_downs.len();
                match *contig_end_exclusion * 2 < len as u64 {
                    true => *total_bases += len as u64 - 2 * *contig_end_exclusion,
                    false => {
                        debug!("Contig too short - less than twice the contig-end-exclusion");
                        return; //contig is all ends, too short
                    }
                }
                let mut cumulative_sum: i32 = 0;
                let start_from = *contig_end_exclusion as usize;
                let end_at = len - *contig_end_exclusion as usize - 1;
                for i in 0..len {
                    let current = ups_and_downs[i as usize];
                    cumulative_sum += current;
                    if i >= start_from && i <= end_at {
                        if cumulative_sum > 0 {
                            *num_covered_bases += 1
                        }
                        *total_count += cumulative_sum as u64;
                    }
                }
                debug!(
                    "After adding contig, have total_count {}, total_bases {}, \
                        num_covered_bases {}, mismatches {}",
                    total_count, total_bases, num_covered_bases, total_mismatches_in_contig
                );
            }
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                contig_end_exclusion,
                ..
            }
            | CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                contig_end_exclusion,
                ..
            }
            | CoverageEstimator::VarianceGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                contig_end_exclusion,
                ..
            } => {
                *num_mapped_reads = num_mapped_reads_in_contig;
                let len1 = ups_and_downs.len();
                match *contig_end_exclusion * 2 < len1 as u64 {
                    true => {
                        debug!("Adding len1 {}", len1);
                        *observed_contig_length += len1 as u64 - 2 * *contig_end_exclusion
                    }
                    false => {
                        debug!("Contig too short - less than twice the contig-end-exclusion");
                        return; //contig is all ends, too short
                    }
                }
                debug!("Total observed length now {}", *observed_contig_length);
                let mut cumulative_sum: i32 = 0;
                let start_from = *contig_end_exclusion as usize;
                let end_at = len1 - *contig_end_exclusion as usize - 1;
                debug!("ups and downs {:?}", ups_and_downs);
                for (i, current) in ups_and_downs.iter().enumerate() {
                    if *current != 0 {
                        debug!(
                            "At i {}, cumulative sum {} and current {}",
                            i, cumulative_sum, current
                        );
                    }
                    cumulative_sum += current;
                    if i >= start_from && i <= end_at {
                        if cumulative_sum > 0 {
                            *num_covered_bases += 1
                        }
                        if counts.len() <= cumulative_sum as usize {
                            (*counts).resize(cumulative_sum as usize + 1, 0);
                        }
                        (*counts)[cumulative_sum as usize] += 1
                    }
                }
            }
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::NumCoveredBasesCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::RPKMCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::TPMCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            } => {
                *num_mapped_reads += num_mapped_reads_in_contig;
                let len = ups_and_downs.len();
                *total_bases += len as u64;
                let mut cumulative_sum: i32 = 0;

                for i in 0..len {
                    let current = ups_and_downs[i as usize];
                    cumulative_sum += current;
                    if cumulative_sum > 0 {
                        *num_covered_bases += 1
                    }
                }
            }
            CoverageEstimator::ReferenceLengthCalculator {
                ref mut observed_contig_length,
                ref mut num_mapped_reads,
            }
            | CoverageEstimator::ReadsPerBaseCalculator {
                ref mut observed_contig_length,
                ref mut num_mapped_reads,
            } => {
                *observed_contig_length += ups_and_downs.len() as u64;
                *num_mapped_reads += num_mapped_reads_in_contig;
            }
            CoverageEstimator::ReadCountCalculator {
                ref mut num_mapped_reads,
            } => {
                *num_mapped_reads += num_mapped_reads_in_contig;
            }
        }
    }

    fn calculate_coverage(&mut self, unobserved_contig_lengths: &Vec<u64>) -> f32 {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                total_count,
                total_bases,
                num_covered_bases,
                num_mapped_reads: _,
                total_mismatches,
                contig_end_exclusion,
                min_fraction_covered_bases,
                exclude_mismatches,
            } => {
                let final_total_bases = *total_bases
                    + CoverageEstimator::calculate_unobserved_bases(
                        unobserved_contig_lengths,
                        *contig_end_exclusion,
                    );
                debug!(
                    "Calculating coverage with unobserved {:?}, \
                        total bases {}, num_covered_bases {}, total_count {}, \
                        total_mismatches {}, final_total_bases {}",
                    unobserved_contig_lengths,
                    total_bases,
                    num_covered_bases,
                    total_count,
                    total_mismatches,
                    final_total_bases,
                );
                if final_total_bases == 0
                    || (*num_covered_bases as f32 / final_total_bases as f32)
                        < *min_fraction_covered_bases
                {
                    return 0.0;
                } else {
                    let calculated_coverage = match exclude_mismatches {
                        true => (*total_count - (*total_mismatches as u64)) as f32,
                        false => *total_count as f32,
                    } / final_total_bases as f32;
                    debug!("Found mean coverage {}", calculated_coverage);
                    return calculated_coverage;
                }
            }
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
                counts,
                observed_contig_length,
                num_covered_bases,
                num_mapped_reads: _,
                contig_end_exclusion,
                min_fraction_covered_bases,
                min,
                max,
            } => {
                let unobserved_contig_length = CoverageEstimator::calculate_unobserved_bases(
                    unobserved_contig_lengths,
                    *contig_end_exclusion,
                );
                let total_bases = *observed_contig_length + unobserved_contig_length;
                debug!("Calculating coverage with num_covered_bases {}, observed_length {}, unobserved_length {:?} and counts {:?}",
                       num_covered_bases, observed_contig_length, unobserved_contig_lengths, counts);
                let answer = match total_bases {
                    0 => 0.0,
                    _ => {
                        if (*num_covered_bases as f32 / total_bases as f32)
                            < *min_fraction_covered_bases
                        {
                            0.0
                        } else {
                            let min_index: usize = (*min * total_bases as f32).floor() as usize;
                            let max_index: usize = (*max * total_bases as f32).ceil() as usize;
                            if *num_covered_bases == 0 {
                                return 0.0;
                            }
                            counts[0] += unobserved_contig_length;

                            let mut num_accounted_for: usize = 0;
                            let mut total: usize = 0;
                            let mut started = false;
                            let mut i = 0;
                            for num_covered in counts.iter() {
                                num_accounted_for += *num_covered as usize;
                                debug!(
                                    "start: i {}, num_accounted_for {}, total {}, min {}, max {}",
                                    i, num_accounted_for, total, min_index, max_index
                                );
                                if num_accounted_for >= min_index {
                                    debug!("inside");
                                    if started {
                                        if num_accounted_for > max_index {
                                            debug!(
                                                "num_accounted_for {}, *num_covered {}",
                                                num_accounted_for, *num_covered
                                            );
                                            let num_excess =
                                                num_accounted_for - *num_covered as usize;
                                            let num_wanted = match max_index >= num_excess {
                                                true => max_index - num_excess + 1,
                                                false => 0,
                                            };
                                            debug!("num wanted1: {}", num_wanted);
                                            total += num_wanted * i;
                                            break;
                                        } else {
                                            total += *num_covered as usize * i;
                                        }
                                    } else {
                                        if num_accounted_for > max_index {
                                            // all coverages are the same in the trimmed set
                                            total = (max_index - min_index + 1) * i;
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
                                debug!(
                                    "end i {}, num_accounted_for {}, total {}",
                                    i, num_accounted_for, total
                                );

                                i += 1;
                            }
                            total as f32 / (max_index - min_index) as f32
                        }
                    }
                };
                return answer;
            }
            CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                counts: _,
                observed_contig_length,
                num_covered_bases,
                num_mapped_reads: _,
                contig_end_exclusion,
                min_fraction_covered_bases,
            } => {
                // No need to actually calculate any kind of coverage, just return
                // whether any coverage was detected
                match observed_contig_length {
                    0 => 0.0,
                    _ => {
                        let total_bases = *observed_contig_length
                            + CoverageEstimator::calculate_unobserved_bases(
                                unobserved_contig_lengths,
                                *contig_end_exclusion,
                            );
                        if (*num_covered_bases as f32 / total_bases as f32)
                            < *min_fraction_covered_bases
                        {
                            0.0
                        } else {
                            // Hack: Return the number of zero coverage bases as the
                            // coverage, plus 1 so it is definitely non-zero, so that
                            // the print_genome function knows this info.
                            (total_bases - *num_covered_bases + 1) as f32
                        }
                    }
                }
            }
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
                total_bases,
                num_covered_bases,
                num_mapped_reads: _,
                min_fraction_covered_bases,
            } => {
                let final_total_bases: u64 =
                    *total_bases + unobserved_contig_lengths.iter().sum::<u64>();
                if final_total_bases == 0
                    || (*num_covered_bases as f32 / final_total_bases as f32)
                        < *min_fraction_covered_bases
                {
                    return 0.0;
                } else {
                    return *num_covered_bases as f32 / final_total_bases as f32;
                }
            }
            CoverageEstimator::NumCoveredBasesCoverageEstimator {
                total_bases,
                num_covered_bases,
                num_mapped_reads: _,
                min_fraction_covered_bases,
            } => {
                let final_total_bases: u64 =
                    *total_bases + unobserved_contig_lengths.iter().sum::<u64>();
                if final_total_bases == 0
                    || (*num_covered_bases as f32 / final_total_bases as f32)
                        < *min_fraction_covered_bases
                {
                    return 0.0;
                } else {
                    return *num_covered_bases as f32;
                }
            }
            CoverageEstimator::RPKMCoverageEstimator {
                total_bases,
                num_covered_bases,
                num_mapped_reads,
                min_fraction_covered_bases,
            } => {
                let final_total_bases: u64 =
                    *total_bases + unobserved_contig_lengths.iter().sum::<u64>();
                if final_total_bases == 0
                    || (*num_covered_bases as f32 / final_total_bases as f32)
                        < *min_fraction_covered_bases
                {
                    return 0.0;
                } else {
                    // Here we do not know the number of mapped reads total.
                    // Instead we divide by that later.
                    debug!("RPKM: {} {}", num_mapped_reads, final_total_bases);
                    return match final_total_bases == 0 {
                        true => 0.0,
                        false => {
                            (*num_mapped_reads * (10u64.pow(9))) as f32 / final_total_bases as f32
                        }
                    };
                }
            }
            CoverageEstimator::TPMCoverageEstimator {
                total_bases,
                num_covered_bases,
                num_mapped_reads,
                min_fraction_covered_bases,
            } => {
                let final_total_bases: u64 =
                    *total_bases + unobserved_contig_lengths.iter().sum::<u64>();
                if final_total_bases == 0
                    || (*num_covered_bases as f32 / final_total_bases as f32)
                        < *min_fraction_covered_bases
                {
                    return 0.0;
                } else {
                    // Here we do not know the number of mapped reads total.
                    // Instead we divide by that later.
                    debug!("TPM: {} {}", num_mapped_reads, final_total_bases);
                    return match final_total_bases == 0 {
                        true => 0.0,
                        false => {
                            ((*num_mapped_reads as f64).ln() - (final_total_bases as f64).ln())
                                .exp() as f32
                        }
                    };
                }
            }
            CoverageEstimator::VarianceGenomeCoverageEstimator {
                observed_contig_length,
                ref mut counts,
                num_covered_bases,
                num_mapped_reads: _,
                contig_end_exclusion,
                min_fraction_covered_bases,
            } => {
                let unobserved_contig_length = CoverageEstimator::calculate_unobserved_bases(
                    unobserved_contig_lengths,
                    *contig_end_exclusion,
                );
                let total_bases = *observed_contig_length + unobserved_contig_length;
                debug!("Calculating coverage with observed length {}, unobserved_length {:?} and counts {:?}", num_covered_bases, unobserved_contig_lengths, counts);
                match total_bases {
                    0 => 0.0,
                    _ => {
                        if (*num_covered_bases as f32 / total_bases as f32)
                            < *min_fraction_covered_bases
                        {
                            0.0
                        } else if total_bases < 3 {
                            0.0
                        } else if counts.len() == 0 {
                            0.0 // no mapped reads
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
                                if *num_covered == 0 {
                                    continue;
                                }
                                let nc = *num_covered as usize;
                                ex += (x - k) * nc;
                                ex2 += (x - k) * (x - k) * nc;
                            }
                            // Return sample variance not population variance since
                            // almost all MAGs are incomplete.
                            (ex2 as f32 - (ex * ex) as f32 / total_bases as f32)
                                / (total_bases - 1) as f32
                        }
                    }
                }
            }
            CoverageEstimator::ReferenceLengthCalculator {
                observed_contig_length,
                ..
            } => (*observed_contig_length + unobserved_contig_lengths.iter().sum::<u64>()) as f32,
            CoverageEstimator::ReadCountCalculator { num_mapped_reads } => *num_mapped_reads as f32,
            CoverageEstimator::ReadsPerBaseCalculator {
                observed_contig_length,
                num_mapped_reads,
            } => {
                *num_mapped_reads as f32
                    / (*observed_contig_length + unobserved_contig_lengths.iter().sum::<u64>())
                        as f32
            }
        }
    }

    fn copy(&self) -> CoverageEstimator {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                total_count: _,
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                total_mismatches: _,
                contig_end_exclusion,
                min_fraction_covered_bases,
                exclude_mismatches,
            } => CoverageEstimator::new_estimator_mean(
                *min_fraction_covered_bases,
                *contig_end_exclusion,
                *exclude_mismatches,
            ),
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
                counts: _,
                observed_contig_length: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                contig_end_exclusion,
                min_fraction_covered_bases,
                min,
                max,
            } => CoverageEstimator::new_estimator_trimmed_mean(
                *min,
                *max,
                *min_fraction_covered_bases,
                *contig_end_exclusion,
            ),
            CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                counts: _,
                observed_contig_length: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                contig_end_exclusion,
                min_fraction_covered_bases,
            } => CoverageEstimator::new_estimator_pileup_counts(
                *min_fraction_covered_bases,
                *contig_end_exclusion,
            ),
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                min_fraction_covered_bases,
            } => CoverageEstimator::new_estimator_covered_fraction(*min_fraction_covered_bases),
            CoverageEstimator::NumCoveredBasesCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                min_fraction_covered_bases,
            } => CoverageEstimator::new_estimator_covered_bases(*min_fraction_covered_bases),
            CoverageEstimator::RPKMCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                min_fraction_covered_bases,
            } => CoverageEstimator::new_estimator_rpkm(*min_fraction_covered_bases),
            CoverageEstimator::TPMCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                min_fraction_covered_bases,
            } => CoverageEstimator::new_estimator_tpm(*min_fraction_covered_bases),
            CoverageEstimator::VarianceGenomeCoverageEstimator {
                observed_contig_length: _,
                counts: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                contig_end_exclusion,
                min_fraction_covered_bases,
            } => CoverageEstimator::new_estimator_variance(
                *min_fraction_covered_bases,
                *contig_end_exclusion,
            ),
            CoverageEstimator::ReferenceLengthCalculator { .. } => {
                CoverageEstimator::new_estimator_length()
            }
            CoverageEstimator::ReadCountCalculator { .. } => {
                CoverageEstimator::new_estimator_read_count()
            }
            CoverageEstimator::ReadsPerBaseCalculator { .. } => {
                CoverageEstimator::new_estimator_reads_per_base()
            }
        }
    }

    fn print_coverage<T: CoverageTaker>(&self, coverage: &f32, coverage_taker: &mut T) {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator { .. }
            | CoverageEstimator::TrimmedMeanGenomeCoverageEstimator { .. }
            | CoverageEstimator::CoverageFractionGenomeCoverageEstimator { .. }
            | CoverageEstimator::NumCoveredBasesCoverageEstimator { .. }
            | CoverageEstimator::RPKMCoverageEstimator { .. }
            | CoverageEstimator::TPMCoverageEstimator { .. }
            | CoverageEstimator::VarianceGenomeCoverageEstimator { .. }
            | CoverageEstimator::ReferenceLengthCalculator { .. }
            | CoverageEstimator::ReadCountCalculator { .. }
            | CoverageEstimator::ReadsPerBaseCalculator { .. } => {
                coverage_taker.add_single_coverage(*coverage);
            }
            CoverageEstimator::PileupCountsGenomeCoverageEstimator { counts, .. } => {
                let mut i = 0;
                debug!("{:?}", counts);
                for num_covered in counts.iter() {
                    let cov: u64 = match i {
                        0 => {
                            let c = coverage.floor() as u64;
                            match c {
                                0 => 0,
                                _ => c - 1,
                            }
                        }
                        _ => *num_covered,
                    };
                    coverage_taker.add_coverage_entry(i, cov);
                    i += 1
                }
            }
        }
    }

    fn print_zero_coverage<T: CoverageTaker>(&self, coverage_taker: &mut T, entry_length: u64) {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator { .. }
            | CoverageEstimator::TrimmedMeanGenomeCoverageEstimator { .. }
            | CoverageEstimator::CoverageFractionGenomeCoverageEstimator { .. }
            | CoverageEstimator::NumCoveredBasesCoverageEstimator { .. }
            | CoverageEstimator::RPKMCoverageEstimator { .. }
            | CoverageEstimator::TPMCoverageEstimator { .. }
            | CoverageEstimator::VarianceGenomeCoverageEstimator { .. }
            | CoverageEstimator::ReadCountCalculator { .. }
            | CoverageEstimator::ReadsPerBaseCalculator { .. } => {
                coverage_taker.add_single_coverage(0.0);
            }
            CoverageEstimator::PileupCountsGenomeCoverageEstimator { .. } => {}
            CoverageEstimator::ReferenceLengthCalculator { .. } => {
                coverage_taker.add_single_coverage(entry_length as f32);
            }
        }
    }

    fn num_mapped_reads(&self) -> u64 {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                total_count: _,
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
                counts: _,
                observed_contig_length: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                counts: _,
                observed_contig_length: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::NumCoveredBasesCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::RPKMCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::TPMCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::VarianceGenomeCoverageEstimator {
                observed_contig_length: _,
                counts: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::ReferenceLengthCalculator {
                observed_contig_length: _,
                num_mapped_reads,
            }
            | CoverageEstimator::ReadCountCalculator { num_mapped_reads }
            | CoverageEstimator::ReadsPerBaseCalculator {
                observed_contig_length: _,
                num_mapped_reads,
            } => *num_mapped_reads,
        }
    }
}
