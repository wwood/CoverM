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
extern crate version_compare;

pub const CONCATENATED_FASTA_FILE_SEPARATOR: &str = "~";

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
