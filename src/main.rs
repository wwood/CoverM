extern crate coverm;
use coverm::mosdepth_genome_coverage_estimators::*;

use std::env;
use std::str;
use std::process;

extern crate clap;
use clap::*;

#[macro_use]
extern crate log;
use log::LogLevelFilter;
extern crate env_logger;
use env_logger::LogBuilder;

fn main(){
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    let get_trimmed_mean_estimator = |m: &clap::ArgMatches, min_fraction_covered| {
        let min = value_t!(m.value_of("trim-min"), f32).unwrap();
        let max = value_t!(m.value_of("trim-max"), f32).unwrap();
        if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
            eprintln!("error: Trim bounds must be between 0 and 1, and min must be less than max, found {} and {}", min, max);
            process::exit(1)
        }
        TrimmedMeanGenomeCoverageEstimator::new(
            min, max, min_fraction_covered)
    };

    match matches.subcommand_name() {

        Some("genome") => {
            let m = matches.subcommand_matches("genome").unwrap();

            let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
            set_log_level(m);
            let method = m.value_of("method").unwrap();
            let min_fraction_covered = value_t!(m.value_of("min-covered-fraction"), f32).unwrap();
            if min_fraction_covered > 1.0 || min_fraction_covered < 0.0 {
                eprintln!("Minimum fraction covered parameter cannot be < 0 or > 1, found {}", min_fraction_covered);
                process::exit(1)
            }
            let print_zeros = !m.is_present("no-zeros");

            if m.is_present("separator") {
                let separator_str = m.value_of("separator").unwrap().as_bytes();
                if separator_str.len() != 1 {
                    eprintln!(
                        "error: Separator can only be a single character, found {} ({}).",
                        separator_str.len(),
                        str::from_utf8(separator_str).unwrap());
                    process::exit(1);
                }
                let separator: u8 = separator_str[0];
                match method {
                    "mean" => coverm::genome::mosdepth_genome_coverage(
                        &bam_files,
                        separator,
                        &mut std::io::stdout(),
                        &mut MeanGenomeCoverageEstimator::new(min_fraction_covered),
                        print_zeros),
                    "coverage_histogram" => coverm::genome::mosdepth_genome_coverage(
                        &bam_files,
                        separator,
                        &mut std::io::stdout(),
                        &mut PileupCountsGenomeCoverageEstimator::new(
                            min_fraction_covered),
                        print_zeros),
                    "trimmed_mean" => {
                        coverm::genome::mosdepth_genome_coverage(
                            &bam_files,
                            separator,
                            &mut std::io::stdout(),
                            &mut get_trimmed_mean_estimator(m, min_fraction_covered),
                            print_zeros)},
                    _ => panic!("programming error")
                }
            } else {
                let genomes_and_contigs;
                if m.is_present("genome-fasta-files"){
                    let genome_fasta_files: Vec<&str> = m.values_of("genome-fasta-files").unwrap().collect();
                    genomes_and_contigs = coverm::read_genome_fasta_files(&genome_fasta_files);
                } else if m.is_present("genome-fasta-directory") {
                    let dir = m.value_of("genome-fasta-directory").unwrap();
                    let paths = std::fs::read_dir(dir).unwrap();
                    let mut genome_fasta_files: Vec<String> = vec!();
                    for path in paths {
                        let file = path.unwrap().path();
                        match file.extension() {
                            Some(ext) => {
                                if ext == "fna" {
                                    let s = String::from(file.to_string_lossy());
                                    genome_fasta_files.push(s);
                                } else {
                                    info!("Not using directory entry '{}' as a genome FASTA file",
                                           file.to_str().expect("UTF8 error in filename"));
                                }
                            },
                            None => {
                                info!("Not using directory entry '{}' as a genome FASTA file",
                                       file.to_str().expect("UTF8 error in filename"));
                            }
                        }
                    }
                    let mut strs: Vec<&str> = vec!();
                    for f in &genome_fasta_files {
                        strs.push(f);
                    }
                    genomes_and_contigs = coverm::read_genome_fasta_files(&strs);
                } else {
                    eprintln!("Either a separator (-s) or path(s) to genome FASTA files (with -d or -f) must be given");
                    process::exit(1);
                }
                match method {
                    "mean" => coverm::genome::mosdepth_genome_coverage_with_contig_names(
                        &bam_files,
                        &genomes_and_contigs,
                        &mut std::io::stdout(),
                        &mut MeanGenomeCoverageEstimator::new(min_fraction_covered),
                        print_zeros),
                    "coverage_histogram" => coverm::genome::mosdepth_genome_coverage_with_contig_names(
                        &bam_files,
                        &genomes_and_contigs,
                        &mut std::io::stdout(),
                        &mut PileupCountsGenomeCoverageEstimator::new(
                            min_fraction_covered),
                        print_zeros),
                    "trimmed_mean" => {
                        coverm::genome::mosdepth_genome_coverage_with_contig_names(
                            &bam_files,
                            &genomes_and_contigs,
                            &mut std::io::stdout(),
                            &mut get_trimmed_mean_estimator(m, min_fraction_covered),
                            print_zeros)},
                    _ => panic!("programming error")
                }
            }
        },
        Some("contig") => {
            let m = matches.subcommand_matches("contig").unwrap();
            set_log_level(m);
            let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
            let method = m.value_of("method").unwrap();
            let min_fraction_covered = value_t!(m.value_of("min-covered-fraction"), f32).unwrap();
            let print_zeros = !m.is_present("no-zeros");
            match method {
                "mean" => coverm::contig::contig_coverage(
                    &bam_files,
                    &mut std::io::stdout(),
                    &mut MeanGenomeCoverageEstimator::new(min_fraction_covered),
                    print_zeros),
                "coverage_histogram" => coverm::contig::contig_coverage(
                    &bam_files,
                    &mut std::io::stdout(),
                    &mut PileupCountsGenomeCoverageEstimator::new(
                        min_fraction_covered),
                    print_zeros),
                "trimmed_mean" => {
                    coverm::contig::contig_coverage(
                        &bam_files,
                        &mut std::io::stdout(),
                        &mut get_trimmed_mean_estimator(m, min_fraction_covered),
                        print_zeros)},
                _ => panic!("programming error")
            }
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
    let genome_help: &'static str =
        "coverm genome: Calculate read coverage per-genome

Define the contigs in each genome (required):
   -s, --separator <CHARACTER>           This character separates genome names
                                         from contig names
   -f, --genome-fasta-files <PATH> ..    Path to FASTA files of each genome
   -d, --genome-fasta-directory <PATH>   Directory containing FASTA files of each
                                         genome

Define mapping(s) (required):
   -b, --bam-files <PATH> ..             Path to reference-sorted BAM files

Other arguments (optional):
   -m, --method METHOD                   Method for calculating coverage (mean,
                                         trimmed_mean or coverage_histogram)
                                         [default: mean]
   --min-covered-fraction FRACTION       Genomes with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0.02]
   --trim-min FRACTION                   Remove this smallest fraction of positions
                                         when calculating trimmed_mean
                                         [default: 0.05]
   --trim-max FRACTION                   Maximum fraction for trimmed_mean
                                         calculations [default: 0.95]
   --no-zeros                            Omit printing of genomes that have zero
                                         coverage [default: false]
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Ben J. Woodcroft <benjwoodcroft near gmail.com>
";

    let contig_help: &'static str =
        "coverm contig: Calculate read coverage per-contig

Define mapping(s) (required):
   -b, --bam-files <PATH> ..             Path to reference-sorted BAM files

Other arguments (optional):
   -m, --method METHOD                   Method for calculating coverage (mean,
                                         trimmed_mean or coverage_histogram)
                                         [default: mean]
   --min-covered-fraction FRACTION       Genomes with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0.02]
   --trim-min FRACTION                   Remove this smallest fraction of positions
                                         when calculating trimmed_mean
                                         [default: 0.05]
   --trim-max FRACTION                   Maximum fraction for trimmed_mean
                                         calculations [default: 0.95]
   --no-zeros                            Omit printing of genomes that have zero
                                         coverage [default: false]
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Ben J. Woodcroft <benjwoodcroft near gmail.com>";

    return App::new("coverm")
        .version("0.1.0-pre")
        .author("Ben J. Woodcroft <benjwoodcroft near gmail.com>")
        .about("Mapping coverage analysis for metagenomics")
        .args_from_usage("-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'")
        .help("
Mapping coverage analysis for metagenomics

Usage: coverm <subcommand> ...

Subcommands:
\tcontig\tCalculate coverage of contigs
\tgenome\tCalculate coverage of genomes

Other options:
\t-V, --version\tPrint version information

Ben J. Woodcroft <benjwoodcroft near gmail.com>
")
        .subcommand(
            SubCommand::with_name("genome")
                .about("Calculate coverage of genomes")
                .help(genome_help)
                .arg(Arg::with_name("bam-files")
                     .short("b")
                     .long("bam-files")
                     .multiple(true)
                     .takes_value(true)
                     .required(true))

                .arg(Arg::with_name("separator")
                     .short("s")
                     .long("separator")
                     .takes_value(true))
                .arg(Arg::with_name("genome-fasta-files")
                     .short("f")
                     .long("genome-fasta-files")
                     .multiple(true)
                     .takes_value(true))
                .arg(Arg::with_name("genome-fasta-directory")
                     .short("d")
                     .long("genome-fasta-directory")
                     .takes_value(true))

                .arg(Arg::with_name("method")
                     .short("m")
                     .long("method")
                     .takes_value(true)
                     .possible_values(&["mean", "trimmed_mean", "coverage_histogram"])
                     .default_value("mean"))
                .arg(Arg::with_name("trim-min")
                     .long("trim-min")
                     .default_value("0.05"))
                .arg(Arg::with_name("trim-max")
                     .long("trim-max")
                     .default_value("0.95"))
                .arg(Arg::with_name("min-covered-fraction")
                     .long("min-covered-fraction")
                     .default_value("0.02"))
                .arg(Arg::with_name("no-zeros")
                     .long("no-zeros"))

                .arg(Arg::with_name("verbose")
                     .short("v")
                     .long("verbose"))
                .arg(Arg::with_name("quiet")
                     .short("q")
                     .long("quiet")))
        .subcommand(
            SubCommand::with_name("contig")
                .about("Calculate coverage of contigs")
                .help(contig_help)

                .arg(Arg::with_name("bam-files")
                     .short("b")
                     .long("bam-files")
                     .multiple(true)
                     .takes_value(true)
                     .required(true))

                .arg(Arg::with_name("method")
                     .short("m")
                     .long("method")
                     .takes_value(true)
                     .possible_values(&["mean", "trimmed_mean", "coverage_histogram"])
                     .default_value("mean"))
                .arg(Arg::with_name("min-covered-fraction")
                     .long("min-covered-fraction")
                     .default_value("0.02"))
                .arg(Arg::with_name("trim-min")
                     .long("trim-min")
                     .default_value("0.05"))
                .arg(Arg::with_name("trim-max")
                     .long("trim-max")
                     .default_value("0.95"))
                .arg(Arg::with_name("no-zeros")
                     .long("no-zeros"))
                .arg(Arg::with_name("verbose")
                     .short("v")
                     .long("verbose"))
                .arg(Arg::with_name("quiet")
                     .short("q")
                     .long("quiet")));
}
