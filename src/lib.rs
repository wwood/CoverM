pub mod bam_generator;
pub mod cli;
pub mod contig;
pub mod coverage_printer;
pub mod coverage_takers;
pub mod external_command_checker;
pub mod filter;
pub mod genome;
pub mod genome_exclusion;
pub mod genome_parsing;
pub mod genomes_and_contigs;
pub mod mapping_index_maintenance;
pub mod mapping_parameters;
pub mod mosdepth_genome_coverage_estimators;
pub mod shard_bam_reader;

use rust_htslib::bam::record::Record;
use std::sync::Arc;

extern crate bio;
#[macro_use]
extern crate log;

extern crate env_logger;
extern crate nix;
extern crate rand;
extern crate rust_htslib;
extern crate tempdir;
extern crate tempfile;
#[macro_use]
extern crate serde;
extern crate clap;
#[macro_use]
extern crate lazy_static;
extern crate bird_tool_utils;
extern crate bird_tool_utils_man;
extern crate galah;
extern crate roff;
extern crate version_compare;

pub const CONCATENATED_FASTA_FILE_SEPARATOR: &str = "~";

pub const AUTHOR: &str =
    "Ben J. Woodcroft, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology";
pub const AUTHOR_AND_EMAIL: &str =
    "Ben J. Woodcroft, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology <benjwoodcroft near gmail.com>";

#[derive(PartialEq, Debug, Clone)]
pub struct ReadsMapped {
    num_mapped_reads: u64,
    num_reads: u64,
}

#[derive(Clone, Debug)]
pub struct FlagFilter {
    pub include_improper_pairs: bool,
    pub include_supplementary: bool,
    pub include_secondary: bool,
}

impl FlagFilter {
    pub fn passes(&self, record: &Record) -> bool {
        if !self.include_secondary && record.is_secondary() {
            return false;
        }
        if !self.include_supplementary && record.is_supplementary() {
            return false;
        }
        if !self.include_improper_pairs && !record.is_proper_pair() {
            return false;
        }
        return true;
    }
}

pub struct OutputWriter {
    pub output_file: Option<Arc<std::sync::Mutex<std::fs::File>>>,
}

impl OutputWriter {
    pub fn generate(file_path_opt: Option<&str>) -> OutputWriter {
        debug!("Generating output writer from path {:?}", file_path_opt);
        match file_path_opt {
            Some(file_path) => {
                if file_path == "-" {
                    info!("Outputing to STDOUT");
                    OutputWriter { output_file: None }
                } else {
                    let path = std::path::Path::new(&file_path);
                    info!("Writing output to file: {}", file_path);
                    OutputWriter {
                        output_file: Some(Arc::new(std::sync::Mutex::new(
                            std::fs::File::create(path)
                                .expect(&format!("Failed to create output file: {}", file_path)),
                        ))),
                    }
                }
            }
            None => {
                debug!("Outputing to STDOUT");
                OutputWriter { output_file: None }
            }
        }
    }
}

impl std::io::Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match &mut self.output_file {
            None => std::io::stdout().write(buf),
            Some(f) => f.lock().expect("failed to unlock output file").write(buf),
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        match &mut self.output_file {
            None => std::io::stdout().flush(),
            Some(f) => f.lock().expect("failed to unlock output file").flush(),
        }
    }
}

impl Clone for OutputWriter {
    fn clone(&self) -> OutputWriter {
        OutputWriter {
            output_file: match &self.output_file {
                Some(f) => Some(f.clone()),
                None => None,
            },
        }
    }
}
