use std;

use rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read as bamRead;

use std::str;
use log;







pub trait GenomeCoverageEstimator {
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
impl GenomeCoverageEstimator for MeanGenomeCoverageEstimator {
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


pub fn genome_coverage2<T: GenomeCoverageEstimator>(
    bam_files: &Vec<&str>,
    split_char: u8,
    print_stream: &mut std::io::Write,
    coverage_estimator: &mut T,
    print_zero_coverage_genomes: bool) {

    for bam_file in bam_files {
        debug!("Working on BAM file {}", bam_file);
        let mut bam = bam::Reader::from_path(bam_file).expect(
            &format!("Unable to find BAM file {}", bam_file));
        let header = bam.header().clone();
        let target_names = header.target_names();

        let fill_genome_length_forwards = |current_tid, target_genome| {
            // Iterating reads skips over contigs with no mapped reads, but the
            // length of these contigs is required to calculate the average
            // across all contigs. This closure returns the number of bases in
            // contigs with tid > current_tid that are part of the current
            // genome.
            let mut extra: u32 = 0;
            let total_refs = header.target_count();
            let mut my_tid = current_tid + 1;
            while my_tid < total_refs {
                let my_genome = extract_genome(my_tid, &target_names, split_char);
                if my_genome == target_genome {
                    extra += header.target_len(my_tid).expect("Malformed bam header or programming error encountered");
                    my_tid += 1;
                } else {
                    break;
                }
            }
            return extra
        };
        let fill_genome_length_backwards = |current_tid, target_genome| {
            if current_tid == 0 {return 0}
            let mut extra: u32 = 0;
            let mut my_tid = current_tid - 1;
            loop {
                let my_genome = extract_genome(my_tid, &target_names, split_char);
                if my_genome == target_genome {
                    extra += header.target_len(my_tid).expect("Malformed bam header or programming error encountered");
                    if my_tid == 0 {
                        break
                    } else {
                        my_tid -= 1;
                    }
                } else {
                    break;
                }
            }
            return extra
        };
        let fill_genome_length_backwards_to_last = |current_tid, last_tid, target_genome| {
            if current_tid == 0 {return 0};
            let mut extra: u32 = 0;
            let mut my_tid = current_tid - 1;
            while my_tid > last_tid {
                let my_genome = extract_genome(my_tid, &target_names, split_char);
                if my_genome == target_genome {
                    extra += header.target_len(my_tid).expect("Malformed bam header or programming error encountered");
                    my_tid -= 1;
                } else {
                    break;
                }
            }
            return extra
        };


        let mut last_tid: u32 = 0;
        let mut doing_first = true;
        let mut last_genome: &[u8] = "error genome".as_bytes();
        let mut unobserved_contig_length: u32 = 0;
        let mut ups_and_downs: Vec<i32> = Vec::new();
        let stoit_name = std::path::Path::new(bam_file).file_stem().unwrap().to_str().expect(
            "failure to convert bam file name to stoit name - UTF8 error maybe?");
        let mut record: bam::record::Record = bam::record::Record::new();
        while bam.read(&mut record).is_ok() {
            // if reference has changed, finish a genome or not
            let tid = record.tid() as u32;
            let current_genome = extract_genome(tid as u32, &target_names, split_char);
            if tid != last_tid || doing_first {
                if doing_first == true {
                    coverage_estimator.setup();
                    last_genome = current_genome;
                    unobserved_contig_length = fill_genome_length_backwards(tid, current_genome);
                    doing_first = false;
                    if print_zero_coverage_genomes {
                        print_previous_zero_coverage_genomes2(
                            stoit_name, b"", current_genome, tid, coverage_estimator,
                            &target_names, split_char, print_stream);
                    }

                } else if current_genome == last_genome {
                    coverage_estimator.add_contig(&ups_and_downs);
                    // Collect the length of reference sequences from this
                    // genome that had no hits that were just skipped over.
                    unobserved_contig_length += fill_genome_length_backwards_to_last(
                        tid, last_tid as u32, current_genome);

                } else {
                    coverage_estimator.add_contig(&ups_and_downs);
                    // Collect the length of refs from the end of the last genome that had no hits
                    unobserved_contig_length += fill_genome_length_backwards_to_last(
                        tid, last_tid as u32, current_genome);
                    // Determine coverage of previous genome
                    let coverage = coverage_estimator.calculate_coverage(unobserved_contig_length);

                    // Print coverage of previous genome
                    if coverage > 0.0 {
                        coverage_estimator.print_genome(
                            &stoit_name,
                            &str::from_utf8(last_genome).unwrap(),
                            &coverage,
                            print_stream);
                    }
                    coverage_estimator.setup();
                    if print_zero_coverage_genomes {
                        print_previous_zero_coverage_genomes2(
                            stoit_name, last_genome, current_genome, tid, coverage_estimator,
                            &target_names, split_char, print_stream);
                    }
                    last_genome = current_genome;
                    unobserved_contig_length = fill_genome_length_backwards(tid, current_genome);
                }

                ups_and_downs = vec![0; header.target_len(tid as u32).expect("Corrupt BAM file?") as usize];
                last_tid = tid;
            }


            // Add coverage info for the current record
            // for each chunk of the cigar string
            debug!("read name {:?}", std::str::from_utf8(record.qname()).unwrap());
            let mut cursor: usize = record.pos() as usize;
            for cig in record.cigar().iter() {
                debug!("Found cigar {:} from {}", cig, cursor);
                match cig.char() {
                    'M' => {
                        // if M, increment start and decrement end index
                        debug!("Adding M at {} and {}", cursor, cursor + cig.len() as usize);
                        ups_and_downs[cursor] += 1;
                        let final_pos = cursor + cig.len() as usize;
                        if final_pos < ups_and_downs.len(){ // True unless the read hits the contig end.
                            ups_and_downs[final_pos] -= 1;
                        }
                        cursor += cig.len() as usize;
                    },
                    'D' => {
                        // if D, move the cursor
                        cursor += cig.len() as usize;
                    },
                    '=' => panic!("CIGAR '=' detected, but this case is not correctly handled for now"),
                    _ => {}
                }
            }
        }
        // Print the last genome
        coverage_estimator.add_contig(&ups_and_downs);
        // Collect the length of refs from the end of the last genome that had no hits
        unobserved_contig_length += fill_genome_length_forwards(last_tid, last_genome);
        // Determine coverage of previous genome
        let coverage = coverage_estimator.calculate_coverage(unobserved_contig_length);

        // Print coverage of previous genome
        if coverage > 0.0 {
            coverage_estimator.print_genome(
                &stoit_name,
                &str::from_utf8(last_genome).unwrap(),
                &coverage,
                print_stream);
        }
        if print_zero_coverage_genomes {
            print_previous_zero_coverage_genomes2(
                stoit_name, last_genome, b"", header.target_count()-1, coverage_estimator,
                &target_names, split_char, print_stream);
        }
    }
}




pub trait PileupGenomeCoverageEstimator {
    fn setup_genome(&mut self);

    #[inline]
    fn process_pileup(&mut self, pileup: rust_htslib::bam::pileup::Pileup);

    #[inline]
    fn finish_genome(&mut self, total_bases: u32) -> f32;

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

pub struct PileupMeanEstimator {
    total_count: u32,
    min_num_covered_bases: f32,
    num_covered_bases: u32
}
impl PileupMeanEstimator {
    pub fn new(min_num_covered_bases: f32) -> PileupMeanEstimator {
        PileupMeanEstimator {
            total_count: 0,
            min_num_covered_bases: min_num_covered_bases,
            num_covered_bases: 0}
    }
}
impl PileupGenomeCoverageEstimator for PileupMeanEstimator {
    fn setup_genome(&mut self) {
        self.total_count = 0;
        self.num_covered_bases = 0;
    }
    fn process_pileup(&mut self, pileup: rust_htslib::bam::pileup::Pileup) {
        self.total_count += pileup.depth();
        self.num_covered_bases += 1;
    }
    fn finish_genome(&mut self, total_bases: u32) -> f32 {
        let answer = match total_bases {
            0 => 0.0,
            _ => {
                if self.num_covered_bases as f32 / total_bases as f32 >= self.min_num_covered_bases {
                    self.total_count as f32 / total_bases as f32
                } else {
                    0.0
                }
            }
        };
        return answer
    }
}

pub struct PileupTrimmedMeanEstimator {
    counts: Vec<u32>,
    min: f32,
    max: f32,
    min_fraction_covered_bases: f32
}
impl PileupTrimmedMeanEstimator {
    pub fn new(min: f32, max: f32, min_fraction_covered_bases: f32) -> PileupTrimmedMeanEstimator {
        PileupTrimmedMeanEstimator {
            counts: vec!(),
            min: min,
            max: max,
            min_fraction_covered_bases: min_fraction_covered_bases}
    }
}
impl PileupGenomeCoverageEstimator for PileupTrimmedMeanEstimator {
    fn setup_genome(&mut self) {
        self.counts = vec!();
    }
    fn process_pileup(&mut self, pileup: rust_htslib::bam::pileup::Pileup) {
        self.counts.push(pileup.depth());
    }
    fn finish_genome(&mut self, total_bases: u32) -> f32 {
        let answer = match total_bases {
            0 => 0.0,
            _ => {
                if (self.counts.len() as f32 / total_bases as f32) < self.min_fraction_covered_bases {
                    0.0
                } else {
                    let num_zeroes = total_bases - self.counts.len() as u32;
                    self.counts.extend(vec![0; num_zeroes as usize]);
                    self.counts.sort_unstable();
                    let min_index: usize = (self.min * self.counts.len() as f32).floor() as usize;
                    let max_index: usize = (self.max * self.counts.len() as f32).ceil() as usize;
                    let total = self.counts[min_index..(max_index+1)].iter().fold(0, |acc, &x| acc + x);
                    debug!("Min {}, Max {} total {}", min_index, max_index, total);
                    total as f32 / (max_index-min_index) as f32
                }
            }
        };
        return answer
    }
}

/// Potentially optimised version of PileupTrimmedMeanEstimator, where the
/// internal storage is a vector of "count counts". However, it has now just
/// been co-opted to print out the number of bases at each coverage level, per
/// genome.
pub struct PileupTrimmedMeanEstimator2 {
    counts: Vec<u32>,
    min: f32,
    max: f32,
    min_fraction_covered_bases: f32
}
impl PileupTrimmedMeanEstimator2 {
    pub fn new(min: f32, max: f32, min_fraction_covered_bases: f32) -> PileupTrimmedMeanEstimator2 {
        PileupTrimmedMeanEstimator2 {
            counts: vec!(),
            min: min,
            max: max,
            min_fraction_covered_bases: min_fraction_covered_bases}
    }
}
impl PileupGenomeCoverageEstimator for PileupTrimmedMeanEstimator2 {
    fn setup_genome(&mut self) {
        self.counts = vec!();
    }
    fn process_pileup(&mut self, pileup: rust_htslib::bam::pileup::Pileup) {
        let depth = pileup.depth() as usize;
        if self.counts.len() <= depth {
            self.counts.resize(depth+1, 0);
        }
        self.counts[depth] += 1;
    }
    fn finish_genome(&mut self, total_bases: u32) -> f32 {
        let answer = match total_bases {
            0 => 0.0,
            _ => {
                let num_bases_covered = self.counts.iter().fold(0, |acc, &x| acc + x);
                if (num_bases_covered as f32 / total_bases as f32) < self.min_fraction_covered_bases {
                    0.0
                } else {
                    let min_index: usize = (self.min * total_bases as f32).floor() as usize;
                    let max_index: usize = (self.max * total_bases as f32).ceil() as usize;
                    if num_bases_covered == 0 {return 0.0;}
                    let num_zeroes = total_bases - num_bases_covered;
                    self.counts[0] = num_zeroes;
                    if log_enabled!(log::LogLevel::Debug){
                        let mut i = 0;
                        for num_covered in self.counts.iter() {
                            debug!("per_base_coverage: {:}\t{:}", i, *num_covered);
                            i += 1
                        }
                    }

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

    fn print_genome<'a >(&self, stoit_name: &str, genome: &str, _coverage: &f32,
                         print_stream: &'a mut std::io::Write) -> &'a mut std::io::Write {
        let mut i = 0;
        debug!("starting to print {}", genome);
        debug!("{:?}",self.counts);
        for num_covered in self.counts.iter() {
            writeln!(print_stream, "{}\t{}\t{:}\t{:}", stoit_name, genome, i, *num_covered).unwrap();
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

fn extract_genome<'a>(tid: u32, target_names: &'a Vec<&[u8]>, split_char: u8) -> &'a [u8] {
    let target_name = target_names[tid as usize];
    debug!("target name {:?}, separator {:?}", target_name, split_char);
    let offset = find_first(target_name, split_char).expect(
        &format!("Contig name {} does not contain split symbol, so cannot determine which genome it belongs to",
                 str::from_utf8(target_name).unwrap()));
    return &target_name[(0..offset)];
}

// Print zero coverage for genomes that have no reads mapped. Genomes are
// detected from the header, counting backwards from the current tid until the
// last seen genome is encountered, or we reach the beginning of the tid array.
fn print_previous_zero_coverage_genomes<'a>(
    stoit_name: &str,
    last_genome: &[u8],
    current_genome: &[u8],
    current_tid: u32,
    pileup_coverage_estimator: &'a PileupGenomeCoverageEstimator,
    target_names: &Vec<&[u8]>,
    split_char: u8,
    print_stream: &mut std::io::Write)
    -> &'a PileupGenomeCoverageEstimator {

    let mut my_current_genome = current_genome;
    let mut tid = current_tid;
    while tid > 0 {
        let genome = extract_genome(tid, &target_names, split_char);
        if genome == last_genome { break; }
        else if genome != my_current_genome {
            // In-between genome encountered for the first time.
            my_current_genome = genome;
            pileup_coverage_estimator.print_zero_coverage(
                &stoit_name, &str::from_utf8(genome).unwrap(), print_stream);
        }
        tid = tid - 1;
    };
    return pileup_coverage_estimator;
}

// Print zero coverage for genomes that have no reads mapped. Genomes are
// detected from the header, counting backwards from the current tid until the
// last seen genome is encountered, or we reach the beginning of the tid array.
fn print_previous_zero_coverage_genomes2<'a>(
    stoit_name: &str,
    last_genome: &[u8],
    current_genome: &[u8],
    current_tid: u32,
    pileup_coverage_estimator: &'a GenomeCoverageEstimator,
    target_names: &Vec<&[u8]>,
    split_char: u8,
    print_stream: &mut std::io::Write)
    -> &'a GenomeCoverageEstimator {

    let mut my_current_genome = current_genome;
    let mut tid = current_tid;
    while tid > 0 {
        let genome = extract_genome(tid, &target_names, split_char);
        if genome == last_genome { break; }
        else if genome != my_current_genome {
            // In-between genome encountered for the first time.
            my_current_genome = genome;
            pileup_coverage_estimator.print_zero_coverage(
                &stoit_name, &str::from_utf8(genome).unwrap(), print_stream);
        }
        tid = tid - 1;
    };
    return pileup_coverage_estimator;
}
pub fn genome_coverage<T: PileupGenomeCoverageEstimator>(
    bam_files: &Vec<&str>,
    split_char: u8,
    print_stream: &mut std::io::Write,
    pileup_coverage_estimator: &mut T,
    print_zero_coverage_genomes: bool) {

    for bam_file in bam_files {
        debug!("Working on BAM file {}", bam_file);
        let mut bam = bam::Reader::from_path(bam_file).expect(
            &format!("Unable to find BAM file {}", bam_file));
        let header = bam.header().clone();
        let target_names = header.target_names();

        let fill_genome_length_forwards = |current_tid, target_genome| {
            // pileup skips over contigs with no mapped reads, but the length of
            // these contigs is required to calculate the average across all
            // contigs. This closure returns the number of bases in contigs with
            // tid > current_tid that are part of the current genome.
            let mut extra: u32 = 0;
            let total_refs = header.target_count();
            let mut my_tid = current_tid + 1;
            while my_tid < total_refs {
                let my_genome = extract_genome(my_tid, &target_names, split_char);
                if my_genome == target_genome {
                    extra += header.target_len(my_tid).expect("Malformed bam header or programming error encountered");
                    my_tid += 1;
                } else {
                    break;
                }
            }
            return extra
        };
        let fill_genome_length_backwards = |current_tid, target_genome| {
            if current_tid == 0 {return 0}
            let mut extra: u32 = 0;
            let mut my_tid = current_tid - 1;
            loop { //my_tid >= 0 unnecessary due to type limits
                let my_genome = extract_genome(my_tid, &target_names, split_char);
                if my_genome == target_genome {
                    extra += header.target_len(my_tid).expect("Malformed bam header or programming error encountered");
                    if my_tid == 0 {
                        break
                    } else {
                        my_tid -= 1;
                    }
                } else {
                    break;
                }
            }
            return extra
        };
        let fill_genome_length_backwards_to_last = |current_tid, last_tid, target_genome| {
            if current_tid == 0 {return 0};
            let mut extra: u32 = 0;
            let mut my_tid = current_tid - 1;
            while my_tid > last_tid {
                let my_genome = extract_genome(my_tid, &target_names, split_char);
                if my_genome == target_genome {
                    extra += header.target_len(my_tid).expect("Malformed bam header or programming error encountered");
                    my_tid -= 1;
                } else {
                    break;
                }
            }
            return extra
        };

        let mut current_genome_length: u32 = 0;
        let mut last_tid: u32 = 0;
        let mut doing_first = true;
        let mut last_genome: &[u8] = "error genome".as_bytes();
        let stoit_name = std::path::Path::new(bam_file).file_stem().unwrap().to_str().expect(
            "failure to convert bam file name to stoit name - UTF8 error maybe?");
        pileup_coverage_estimator.setup_genome();
        for p in bam.pileup() {
            let pileup = p.unwrap();

            let tid = pileup.tid();
            if tid != last_tid || doing_first {
                // find the new name of the genome
                let current_genome = extract_genome(tid, &target_names, split_char);
                debug!("Found current genome {:?}", str::from_utf8(current_genome).unwrap());

                if doing_first == true {
                    last_genome = current_genome;
                    current_genome_length = fill_genome_length_backwards(tid, current_genome);
                    doing_first = false;
                    if print_zero_coverage_genomes {
                        print_previous_zero_coverage_genomes(
                            stoit_name, b"", current_genome, tid, pileup_coverage_estimator,
                            &target_names, split_char, print_stream);
                    }
                } else if current_genome == last_genome {
                    current_genome_length += fill_genome_length_backwards_to_last(tid, last_tid, current_genome);
                } else {
                    current_genome_length += fill_genome_length_forwards(last_tid, last_genome);
                    let coverage = pileup_coverage_estimator.finish_genome(current_genome_length);
                    if coverage > 0.0 {
                        pileup_coverage_estimator.print_genome(&stoit_name,
                                                               &str::from_utf8(last_genome).unwrap(),
                                                               &coverage,
                                                               print_stream);
                    }
                    pileup_coverage_estimator.setup_genome();
                    if print_zero_coverage_genomes {
                        print_previous_zero_coverage_genomes(
                            stoit_name, last_genome, current_genome, tid, pileup_coverage_estimator,
                            &target_names, split_char, print_stream);
                    }
                    last_genome = current_genome;
                    current_genome_length = fill_genome_length_backwards(tid, current_genome);
                }

                current_genome_length += header.target_len(tid).expect("malformed bam header");
                debug!("genome length now {}", current_genome_length);
                last_tid = tid;
            }
            pileup_coverage_estimator.process_pileup(pileup);
        }
        current_genome_length += fill_genome_length_forwards(last_tid, last_genome);
        let coverage = pileup_coverage_estimator.finish_genome(current_genome_length);
        if coverage > 0.0 {
            pileup_coverage_estimator.print_genome(&stoit_name, &str::from_utf8(last_genome).unwrap(),
                                                   &coverage, print_stream);
        };
        if print_zero_coverage_genomes {
            print_previous_zero_coverage_genomes(
                stoit_name, last_genome, b"", header.target_count()-1, pileup_coverage_estimator,
                &target_names, split_char, print_stream);
        }
    };
}

/// Finds the first occurence of element in a slice
fn find_first<T>(slice: &[T], element: T) -> Result<usize, &'static str>
    where T: std::cmp::PartialEq<T> {

    let mut index: usize = 0;
    for el in slice {
        if *el == element {
            return Ok(index)
            //let res: Result<usize, None> = Ok(index)
        }
        index += 1;
    }
    return Err("Element not found in slice")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use env_logger;

    #[test]
    fn initialize_logger() {
        env_logger::init().unwrap();
    }

    #[test]
    fn test_one_genome_two_contigs_first_covered(){
        let mut stream = Cursor::new(Vec::new());
        genome_coverage2(
            &vec!["test/data/2seqs.reads_for_seq1.bam"],
            'q' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true);
        assert_eq!(
            "2seqs.reads_for_seq1\tse\t0.6\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_two_contigs_second_covered(){
        let mut stream = Cursor::new(Vec::new());
        genome_coverage2(
            &vec!["test/data/2seqs.reads_for_seq2.bam"],
            'q' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true);
        assert_eq!(
            "2seqs.reads_for_seq2\tse\t0.6\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_two_contigs_both_covered(){
        let mut stream = Cursor::new(Vec::new());
        genome_coverage2(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            'e' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.0),
            true);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_min_fraction_covered_under_min(){
        let mut stream = Cursor::new(Vec::new());
        genome_coverage2(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            'e' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.76),
            true);
        assert_eq!(
            "",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_min_fraction_covered_just_ok(){
        let mut stream = Cursor::new(Vec::new());
        genome_coverage2(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            'e' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.759),
            true);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_two_contigs_trimmed_mean(){
        let mut stream = Cursor::new(Vec::new());
        genome_coverage(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            'e' as u8,
            &mut stream,
            &mut PileupTrimmedMeanEstimator::new(0.1, 0.9, 0.0),
            true);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.08875\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_two_contigs_trimmed_mean2(){
        let mut stream = Cursor::new(Vec::new());
        genome_coverage(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            'e' as u8,
            &mut stream,
            &mut PileupTrimmedMeanEstimator2::new(0.1, 0.9, 0.0),
            true);
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t0\t482\n2seqs.reads_for_seq1_and_seq2\ts\t1\t922\n2seqs.reads_for_seq1_and_seq2\ts\t2\t371\n2seqs.reads_for_seq1_and_seq2\ts\t3\t164\n2seqs.reads_for_seq1_and_seq2\ts\t4\t61\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_zero_coverage_genomes(){
        let mut stream = Cursor::new(Vec::new());
        genome_coverage2(
            &vec!["test/data/7seqs.reads_for_seq1_and_seq2.bam"],
            '~' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.1),
            true);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome1\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome4\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome3\t0.0\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome6\t0.0\n",
            str::from_utf8(stream.get_ref()).unwrap());

        stream = Cursor::new(Vec::new());
        genome_coverage2(
            &vec!["test/data/7seqs.reads_for_seq1_and_seq2.bam"],
            '~' as u8,
            &mut stream,
            &mut MeanGenomeCoverageEstimator::new(0.1),
            false);
        assert_eq!(
            "7seqs.reads_for_seq1_and_seq2\tgenome2\t1.2\n7seqs.reads_for_seq1_and_seq2\tgenome5\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }
}
