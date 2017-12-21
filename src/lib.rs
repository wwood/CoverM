extern crate bio;
#[macro_use]
extern crate log;

extern crate rust_htslib;
use self::rust_htslib::bam;
use self::rust_htslib::bam::Read as bamRead;

use std::str;

extern crate env_logger;
use env_logger::LogBuilder;

pub trait PileupGenomeCoverageEstimator {
    #[inline]
    fn process_pileup(&mut self, pileup: rust_htslib::bam::pileup::Pileup);

    #[inline]
    fn finish_genome(&mut self, total_bases: u32) -> f32;
}

pub struct PileupMeanEstimator {
    total_count: u32
}
impl PileupMeanEstimator {
    pub fn new() -> PileupMeanEstimator {
        PileupMeanEstimator {total_count: 0}
    }
}
impl PileupGenomeCoverageEstimator for PileupMeanEstimator {
    fn process_pileup(&mut self, pileup: rust_htslib::bam::pileup::Pileup) {
        self.total_count += pileup.depth()
    }
    fn finish_genome(&mut self, total_bases: u32) -> f32 {
        let answer = match total_bases {
            0 => 0.0,
            _ => self.total_count as f32 / total_bases as f32
        };
        self.total_count = 0;
        return answer
    }
}

pub struct PileupTrimmedMeanEstimator {
    counts: Vec<u32>,
    min: f32,
    max: f32
}
impl PileupTrimmedMeanEstimator {
    pub fn new() -> PileupTrimmedMeanEstimator {
        PileupTrimmedMeanEstimator {
            counts: vec!(),
            max: 0.9,
            min: 0.1}
    }
}
impl PileupGenomeCoverageEstimator for PileupTrimmedMeanEstimator {
    fn process_pileup(&mut self, pileup: rust_htslib::bam::pileup::Pileup) {
        self.counts.push(pileup.depth())
    }
    fn finish_genome(&mut self, total_bases: u32) -> f32 {
        let answer = match total_bases {
            0 => 0.0,
            _ => {
                let num_zeroes = total_bases - self.counts.len() as u32;
                self.counts.extend(vec![0; num_zeroes as usize]);
                self.counts.sort_unstable();
                let min_index: usize = (self.min * self.counts.len() as f32).floor() as usize;
                let max_index: usize = (self.max * self.counts.len() as f32).ceil() as usize;
                let total = self.counts[min_index..(max_index+1)].iter().fold(0, |acc, &x| acc + x);
                debug!("Min {}, Max {} total {}", min_index, max_index, total);
                total as f32 / (max_index-min_index) as f32
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
        let bam = bam::Reader::from_path(bam_file).expect(
            &format!("Unable to find BAM file {}", bam_file));
        let header = bam.header();
        let target_names = header.target_names();

        let mut print_genome = |stoit_name, genome, coverage| {
            debug!("{:?} {:?}", str::from_utf8(genome).unwrap(), coverage);
            writeln!(print_stream, "{}\t{}\t{}",
                     stoit_name,
                     str::from_utf8(genome).unwrap(),
                     coverage).unwrap();
        };
        let extract_genome = |tid| {
            let target_name = target_names[tid as usize];
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
            &mut PileupMeanEstimator::new());
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
            &mut PileupMeanEstimator::new());
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
            &mut PileupMeanEstimator::new());
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
            &mut PileupTrimmedMeanEstimator::new());
        assert_eq!(
            "2seqs.reads_for_seq1_and_seq2\ts\t1.08875\n",
            str::from_utf8(stream.get_ref()).unwrap())
    }
}
