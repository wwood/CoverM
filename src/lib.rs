extern crate bio;
#[macro_use]
extern crate log;

extern crate rust_htslib;
use self::rust_htslib::bam;
use self::rust_htslib::bam::Read as bamRead;

use std::str;

extern crate env_logger;

pub trait PileupGenomeCoverageEstimator {
    #[inline]
    fn process_pileup(&mut self, pileup: rust_htslib::bam::pileup::Pileup);

    #[inline]
    fn finish_genome(&mut self, total_bases: u32) -> f32;
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
        self.total_count = 0;
        self.num_covered_bases = 0;
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
        self.counts = vec!();
        return answer
    }
}

/// Potentially optimised version of PileupTrimmedMeanEstimator, where the
/// internal storage is a vector of "count counts"
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
    fn process_pileup(&mut self, pileup: rust_htslib::bam::pileup::Pileup) {
        let depth = pileup.depth() as usize;
        if self.counts.len() <= depth {
            self.counts.resize(depth+1, 0);
        }
        self.counts[depth] += 1;
    }
    fn finish_genome(&mut self, total_bases: u32) -> f32 {
        debug!("{:?}", self.counts);
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
                    debug!("{:?}", self.counts);

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
        self.counts = vec!();
        return answer
    }
}

pub fn genome_coverage<T: PileupGenomeCoverageEstimator>(
    bam_files: &Vec<&str>,
    split_char: u8,
    print_stream: &mut std::io::Write,
    pileup_coverage_estimator: &mut T) {

    for bam_file in bam_files {
        debug!("Working on BAM file {}", bam_file);
        let mut bam = bam::Reader::from_path(bam_file).expect(
            &format!("Unable to find BAM file {}", bam_file));
        let header = bam.header().clone();
        let target_names = header.target_names();

        let mut print_genome = |stoit_name, genome, coverage| {
            debug!("{:?} {:?}", str::from_utf8(genome).unwrap(), coverage);
            if coverage > 0.0 {
                writeln!(print_stream, "{}\t{}\t{}",
                         stoit_name,
                         str::from_utf8(genome).unwrap(),
                         coverage).unwrap();
            }
        };
        let extract_genome = |tid| {
            let target_name = target_names[tid as usize];
            debug!("target name {:?}, separator {:?}", target_name, split_char);
            let offset = find_first(target_name, split_char).expect(
                &format!("Contig name {} does not contain split symbol, so cannot determine which genome it belongs to",
                         str::from_utf8(target_name).unwrap()));
            return &target_name[(0..offset)];
        };
        let fill_genome_length_forwards = |current_tid, target_genome| {
            // pileup skips over contigs with no mapped reads, but the length of
            // these contigs is required to calculate the average across all
            // contigs. This closure returns the number of bases in contigs with
            // tid > current_tid that are part of the current genome.
            let mut extra: u32 = 0;
            let total_refs = header.target_count();
            let mut my_tid = current_tid + 1;
            while my_tid < total_refs {
                let my_genome = extract_genome(my_tid);
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
                let my_genome = extract_genome(my_tid);
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
                let my_genome = extract_genome(my_tid);
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
        for p in bam.pileup() {
            let pileup = p.unwrap();

            let tid = pileup.tid();
            if tid != last_tid || doing_first {
                // find the new name of the genome
                let current_genome = extract_genome(tid);
                debug!("Found current genome {:?}", str::from_utf8(current_genome).unwrap());

                if doing_first == true {
                    last_genome = current_genome;
                    current_genome_length = fill_genome_length_backwards(tid, current_genome);
                    doing_first = false;
                } else if current_genome == last_genome {
                    current_genome_length += fill_genome_length_backwards_to_last(tid, last_tid, current_genome);
                } else {
                    current_genome_length += fill_genome_length_forwards(last_tid, last_genome);
                    print_genome(stoit_name, last_genome, pileup_coverage_estimator.finish_genome(current_genome_length));
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
        print_genome(stoit_name, last_genome, pileup_coverage_estimator.finish_genome(current_genome_length));
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

    #[test]
    fn initialize_logger() {
        env_logger::init().unwrap();
    }

    #[test]
    fn test_one_genome_two_contigs_first_covered(){
        let mut stream = Cursor::new(Vec::new());
        genome_coverage(
            &vec!["test/data/2seqs.reads_for_seq1.bam"],
            'q' as u8,
            &mut stream,
            &mut PileupMeanEstimator::new(0.0));
        assert_eq!(
            "2seqs.reads_for_seq1\tse\t0.6\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_two_contigs_second_covered(){
        let mut stream = Cursor::new(Vec::new());
        genome_coverage(
            &vec!["test/data/2seqs.reads_for_seq2.bam"],
            'q' as u8,
            &mut stream,
            &mut PileupMeanEstimator::new(0.0));
        assert_eq!(
            "2seqs.reads_for_seq2\tse\t0.6\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_two_contigs_both_covered(){
        let mut stream = Cursor::new(Vec::new());
        genome_coverage(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            'e' as u8,
            &mut stream,
            &mut PileupMeanEstimator::new(0.0));
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.2\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_min_fraction_covered_under_min(){
        let mut stream = Cursor::new(Vec::new());
        genome_coverage(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            'e' as u8,
            &mut stream,
            &mut PileupMeanEstimator::new(0.76));
        assert_eq!(
            "",
            str::from_utf8(stream.get_ref()).unwrap())
    }

    #[test]
    fn test_one_genome_min_fraction_covered_just_ok(){
        let mut stream = Cursor::new(Vec::new());
        genome_coverage(
            &vec!["test/data/2seqs.reads_for_seq1_and_seq2.bam"],
            'e' as u8,
            &mut stream,
            &mut PileupMeanEstimator::new(0.759));
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
            &mut PileupTrimmedMeanEstimator::new(0.1, 0.9, 0.0));
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
            &mut PileupTrimmedMeanEstimator2::new(0.1, 0.9, 0.0));
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.08875\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }
}
