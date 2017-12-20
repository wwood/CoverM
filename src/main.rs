extern crate coverm;

use std::env;
use std::str;

extern crate clap;
use clap::*;

extern crate log;
use log::LogLevelFilter;
extern crate env_logger;
use env_logger::LogBuilder;

fn main(){
    let mut app = build_cli();
    let matches = app.clone().get_matches();

    match matches.subcommand_name() {
        Some("genome") => {
            let m = matches.subcommand_matches("genome").unwrap();
            let separator: u8 = m.value_of("separator").unwrap().as_bytes()[0];
            let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
            set_log_level(m);
            coverm::genome_coverage(
                &bam_files,
                separator,
                &mut std::io::stdout(),
                &mut coverm::PileupMeanEstimator::new());
        },
        _ => {
            app.print_help().unwrap();
            println!();
        }
    }
}

fn set_log_level(matches: &clap::ArgMatches) {
    let mut log_level = LogLevelFilter::Info;
    if matches.is_present("verbose") {
        log_level = LogLevelFilter::Debug;
    }
    if matches.is_present("quiet") {
        log_level = LogLevelFilter::Error;
    }
    let mut builder = LogBuilder::new();
    builder.filter(None, log_level);
    if env::var("RUST_LOG").is_ok() {
        builder.parse(&env::var("RUST_LOG").unwrap());
    }
    builder.init().unwrap();
}

fn build_cli() -> App<'static, 'static> {
    //-f, --fasta-files=<FILE>...         'Read contig to genome mapping from these fasta files'
    let genome_args: &'static str = "-b, --bam-files=<BAM>...      'Sorted BAM files contain reads mapped to target contigs'
                      -s, --separator=<CHARACTER>         'Character used in contig name to separate genome (first) from contig (second) e.g. ~ for BAM reference names like genome1~contig2'

                      -v, --verbose       'Print extra debug logging information'
                      -q, --quiet         'Unless there is an error, do not print logging information'";

    return App::new("coverm")
        .version("0.1.0-pre")
        .author("Ben J. Woodcroft <benjwoodcroft near gmail.com>")
        .about("Mapping coverage analysis for metagenomics")
        .args_from_usage("-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'")
        .subcommand(
            SubCommand::with_name("genome")
                .about("Calculate coverage of genomes")
                .args_from_usage(&genome_args));
}
