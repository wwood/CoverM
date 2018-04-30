extern crate coverm;

use std::env;
use std::str;
use std::process;

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
            let separator_str = m.value_of("separator").unwrap().as_bytes();
            if separator_str.len() != 1 {
                eprintln!(
                    "error: Separator can only be a single character, found {} ({}).",
                    separator_str.len(),
                    str::from_utf8(separator_str).unwrap());
                process::exit(1);
            }
            let separator: u8 = separator_str[0];
            let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
            set_log_level(m);
            let method = m.value_of("method").unwrap();
            let min_fraction_covered = value_t!(m.value_of("min-covered-fraction"), f32).unwrap();
            if min_fraction_covered > 1.0 || min_fraction_covered < 0.0 {
                eprintln!("Minimum fraction covered parameter cannot be < 0 or > 1, found {}", min_fraction_covered);
                process::exit(1)
            }
            let print_zeros = !m.is_present("no-zeros");
            match method {
                "mean" => coverm::genome_coverage(
                    &bam_files,
                    separator,
                    &mut std::io::stdout(),
                    &mut coverm::PileupMeanEstimator::new(min_fraction_covered),
                    print_zeros),
                _ => {
                    let min = value_t!(m.value_of("trim-min"), f32).unwrap();
                    let max = value_t!(m.value_of("trim-max"), f32).unwrap();
                    if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                        eprintln!("error: Trim bounds must be between 0 and 1, and min must be less than max, found {} and {}", min, max);
                        process::exit(1)
                    }
                    match method {
                        "trimmed_mean" => coverm::genome_coverage(
                            &bam_files,
                            separator,
                            &mut std::io::stdout(),
                            &mut coverm::PileupTrimmedMeanEstimator::new(
                                min, max, min_fraction_covered),
                            print_zeros),
                        "coverage_histogram" => coverm::genome_coverage(
                            &bam_files,
                            separator,
                            &mut std::io::stdout(),
                            &mut coverm::PileupTrimmedMeanEstimator2::new(
                                min, max, min_fraction_covered),
                            print_zeros),
                        _ => panic!("programming error")
                    }
                }
            }
        },
        Some("contig") => {
            let m = matches.subcommand_matches("contig").unwrap();
            set_log_level(m);
            let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
            coverm::contig::contig_coverage(
                &bam_files,
                &mut std::io::stdout());
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
                      -s, --separator=<CHARACTER>         'Character used in contig name to separate genome (first) from contig (second) e.g. '~' for BAM reference names like genome1~contig2'

                      -v, --verbose       'Print extra debug logging information'
                      -q, --quiet         'Unless there is an error, do not print logging information'";
    let contig_args: &'static str = "-b, --bam-files=<BAM>...      'Sorted BAM files contain reads mapped to target contigs'

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
                .args_from_usage(&genome_args)
                .arg(Arg::with_name("method")
                     .short("m")
                     .long("method")
                     .help("Method for calculating coverage")
                     .takes_value(true)
                     .possible_values(&["mean", "trimmed_mean", "coverage_histogram"])
                     .default_value("mean"))
                .arg(Arg::with_name("trim-min")
                     .long("trim-min")
                     .help("Minimum for trimmed mean calculations")
                     .default_value("0.05"))
                .arg(Arg::with_name("trim-max")
                     .long("trim-max")
                     .help("Maximum for trimmed mean calculations")
                     .default_value("0.95"))
                .arg(Arg::with_name("min-covered-fraction")
                     .long("min-covered-fraction")
                     .help("Minimum fraction of the genome covered (when less than this, coverage is set to zero)")
                     .default_value("0.02"))
                .arg(Arg::with_name("no-zeros")
                     .long("no-zeros")
                     .help("Omit printing of genomes that have insufficient coverage")))
        .subcommand(
            SubCommand::with_name("contig")
                .about("Calculate coverage of contigs")
                .args_from_usage(&contig_args));
}
