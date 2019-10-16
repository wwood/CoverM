pub mod contig;
pub mod genome;
pub mod mosdepth_genome_coverage_estimators;
pub mod genomes_and_contigs;
pub mod bam_generator;
pub mod filter;
pub mod external_command_checker;
pub mod mapping_index_maintenance;
pub mod coverage_takers;
pub mod mapping_parameters;
pub mod coverage_printer;
pub mod shard_bam_reader;
pub mod genome_exclusion;
pub mod cli;
pub mod genome_parsing;

extern crate bio;
#[macro_use]
extern crate log;

extern crate rust_htslib;
extern crate env_logger;
extern crate nix;
extern crate tempdir;
extern crate tempfile;
extern crate rand;
#[macro_use]
extern crate serde;
extern crate clap;
#[macro_use]
extern crate lazy_static;

pub const CONCATENATED_FASTA_FILE_SEPARATOR: &str = "~";

#[derive(PartialEq, Debug)]
pub struct ReadsMapped {
    num_mapped_reads: u64,
    num_reads: u64
}

#[derive(Clone, Debug)]
pub struct FlagFilter {
    pub include_improper_pairs: bool,
    pub include_supplementary: bool,
    pub include_secondary: bool,
}
