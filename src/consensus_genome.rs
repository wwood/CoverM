use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use std::collections::{HashMap, HashSet};
use std::fs::{create_dir_all, File};
use std::io::Write;
use std::path::Path;

#[derive(Clone)]
pub struct ConsensusGenomeOutputConfig {
    pub folder: String,
    pub min_coverage: u32,
}

#[derive(Clone, Copy, Default)]
pub struct ConsensusBaseCounts {
    a: u32,
    c: u32,
    g: u32,
    t: u32,
    skipped: u32,
    coverage: u32,
}

impl ConsensusBaseCounts {
    pub fn add_base(&mut self, base: u8) {
        self.coverage += 1;
        match base.to_ascii_uppercase() {
            b'A' => self.a += 1,
            b'C' => self.c += 1,
            b'G' => self.g += 1,
            b'T' => self.t += 1,
            _ => {}
        }
    }

    pub fn add_skipped(&mut self) {
        self.coverage += 1;
        self.skipped += 1;
    }

    pub fn consensus_base(&self, min_coverage: u32) -> u8 {
        if self.coverage < min_coverage {
            return b'N';
        }

        let mut best = 0;
        let mut best_base = b'N';
        if self.a > best {
            best = self.a;
            best_base = b'A';
        }
        if self.c > best {
            best = self.c;
            best_base = b'C';
        }
        if self.g > best {
            best = self.g;
            best_base = b'G';
        }
        if self.t > best {
            best = self.t;
            best_base = b'T';
        }

        if self.skipped > best {
            b'-'
        } else {
            best_base
        }
    }
}

fn write_wrapped_fasta_record(file: &mut File, name: &str, sequence: &[u8]) {
    writeln!(file, ">{name}").expect("Unable to write FASTA header");
    for chunk in sequence.chunks(80) {
        writeln!(file, "{}", std::str::from_utf8(chunk).unwrap())
            .expect("Unable to write FASTA sequence");
    }
}

pub fn safe_genome_file_name(genome_name: &str) -> String {
    genome_name
        .chars()
        .map(|c| match c {
            '/' | '\\' | ':' | '*' | '?' | '"' | '<' | '>' | '|' => '_',
            _ => c,
        })
        .collect()
}

pub fn update_consensus_counts_for_record(
    counts: &mut [ConsensusBaseCounts],
    record: &bam::record::Record,
) {
    let seq = record.seq().as_bytes();
    let mut read_cursor: usize = 0;
    let mut ref_cursor: usize = record.pos() as usize;
    for cig in record.cigar().iter() {
        match cig {
            Cigar::Match(l) | Cigar::Diff(l) | Cigar::Equal(l) => {
                for _ in 0..*l as usize {
                    if ref_cursor < counts.len() && read_cursor < seq.len() {
                        counts[ref_cursor].add_base(seq[read_cursor]);
                    }
                    ref_cursor += 1;
                    read_cursor += 1;
                }
            }
            Cigar::Del(l) | Cigar::RefSkip(l) => {
                for _ in 0..*l as usize {
                    if ref_cursor < counts.len() {
                        counts[ref_cursor].add_skipped();
                    }
                    ref_cursor += 1;
                }
            }
            Cigar::Ins(l) | Cigar::SoftClip(l) => {
                read_cursor += *l as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
}

pub fn write_consensus_genomes_from_counts(
    stoit_name: &str,
    output_config: &ConsensusGenomeOutputConfig,
    genomes_to_write: &[usize],
    genome_names: &[String],
    genome_to_tids: &[Vec<u32>],
    tid_to_counts: &HashMap<u32, Vec<ConsensusBaseCounts>>,
) {
    let mut seen_sanitized = HashMap::new();
    let mut seen_paths = HashSet::new();
    for genome_index in genomes_to_write {
        let genome_name = &genome_names[*genome_index];
        let sanitized = safe_genome_file_name(genome_name);
        if let Some(previous) = seen_sanitized.insert(sanitized.clone(), genome_name.clone()) {
            panic!(
                "Genome names '{}' and '{}' resolve to the same consensus output filename '{}.fna'",
                previous, genome_name, sanitized
            );
        }
        let relative_path = format!("{}/{}.fna", stoit_name, sanitized);
        if !seen_paths.insert(relative_path.clone()) {
            panic!(
                "Duplicate consensus output path detected: {}",
                relative_path
            );
        }
    }

    let sample_dir = Path::new(&output_config.folder).join(stoit_name);
    create_dir_all(&sample_dir).expect("Unable to create consensus genome output directory");
    for genome_index in genomes_to_write {
        let genome_name = &genome_names[*genome_index];
        let output_path = sample_dir.join(format!("{}.fna", safe_genome_file_name(genome_name)));
        let mut output_file = File::create(&output_path).unwrap_or_else(|_| {
            panic!(
                "Unable to create consensus FASTA file for genome '{}' at {}",
                genome_name,
                output_path.display()
            )
        });
        let mut consensus = Vec::new();
        for tid in &genome_to_tids[*genome_index] {
            if let Some(counts) = tid_to_counts.get(tid) {
                for c in counts {
                    consensus.push(c.consensus_base(output_config.min_coverage));
                }
            }
        }
        write_wrapped_fasta_record(&mut output_file, genome_name, &consensus);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_consensus_base_counts_threshold_and_skipped() {
        let mut counts = ConsensusBaseCounts::default();
        for _ in 0..4 {
            counts.add_base(b'A');
        }
        assert_eq!(counts.consensus_base(5), b'N');
        for _ in 0..5 {
            counts.add_skipped();
        }
        assert_eq!(counts.consensus_base(5), b'-');
    }

    #[test]
    fn test_consensus_no_nucleotide_evidence_is_n() {
        let counts = ConsensusBaseCounts::default();
        assert_eq!(counts.consensus_base(0), b'N');

        let mut counts2 = ConsensusBaseCounts::default();
        counts2.add_base(b'N');
        assert_eq!(counts2.consensus_base(0), b'N');
    }

    #[test]
    fn test_safe_genome_file_name() {
        assert_eq!(safe_genome_file_name("a/b:c*?\"<d>|"), "a_b_c____d__");
    }

    #[test]
    #[should_panic(expected = "resolve to the same consensus output filename")]
    fn test_panics_on_sanitized_filename_collision() {
        let cfg = ConsensusGenomeOutputConfig {
            folder: tempfile::tempdir()
                .unwrap()
                .path()
                .to_str()
                .unwrap()
                .to_string(),
            min_coverage: 1,
        };
        let genomes = vec!["a/b".to_string(), "a:b".to_string()];
        let genome_to_tids = vec![vec![0], vec![1]];
        let mut tid_to_counts = HashMap::new();
        tid_to_counts.insert(0, vec![ConsensusBaseCounts::default()]);
        tid_to_counts.insert(1, vec![ConsensusBaseCounts::default()]);
        write_consensus_genomes_from_counts(
            "sample1",
            &cfg,
            &[0, 1],
            &genomes,
            &genome_to_tids,
            &tid_to_counts,
        );
    }
}
