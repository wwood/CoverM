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

    let mut htype = &mut HeaderTypes::created();
    let mut htype_hist = &mut HeaderTypes::created();
    pub fn print_header(htype: &mut HeaderTypes){
        for h in htype.headers.iter(){
            print!("{}\t", h)
        }
        println!("")
    }
    let mut output_stream = Vec::new();
    let mut limit_stream = false;
    pub fn update_outputs(output: &mut Vec<OutputStream>, mut input: Vec<OutputStream>) -> Vec<OutputStream> {
        // let it = output.iter().zip(input.iter());
        if output.len()==0{
            return input
        } else {
            // let cnt = 0;
            let mut output_st: Vec<OutputStream> = Vec::new();
            if input.len() > 0 {
                for (i, v) in input.iter().enumerate() {
                    for m in v.methods.iter(){
                        output[i].methods.push(*m);
                    }
                    output_st.push(output[i].clone());
                    // let mut cnt = cnt + 1;
                }
            } else{
                output[0].methods.append(&mut input[0].methods);
                output_st.push(output[0].clone());
            }
        return output_st
        }
    }
    pub fn print_output_stream(output_stream: Vec<OutputStream>){
        for mut out in output_stream{
            out.print_output();
        }
    }
    match matches.subcommand_name() {

        Some("genome") => {
            let m = matches.subcommand_matches("genome").unwrap();

            let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
            set_log_level(m);
            let method: Vec<&str> = m.values_of("method").unwrap().collect();
            let min_fraction_covered = value_t!(m.value_of("min-covered-fraction"), f32).unwrap();
            if min_fraction_covered > 1.0 || min_fraction_covered < 0.0 {
                eprintln!("Minimum fraction covered parameter cannot be < 0 or > 1, found {}", min_fraction_covered);
                process::exit(1)
            }
            let print_zeros = !m.is_present("no-zeros");
            let flag_filter = !m.is_present("no-flag-filter");
            let single_genome = m.is_present("single-genome");
            let headers = !m.is_present("remove-headers");

            if m.is_present("separator") || single_genome {
                let separator: u8 = match single_genome {
                    true => "0".as_bytes()[0],
                    false => {
                        let separator_str = m.value_of("separator").unwrap().as_bytes();
                        if separator_str.len() != 1 {
                            eprintln!(
                                "error: Separator can only be a single character, found {} ({}).",
                                separator_str.len(),
                                str::from_utf8(separator_str).unwrap());
                            process::exit(1);
                        }
                        separator_str[0]
                    }
                };
                for cover in method{
                    match cover {
                        "mean" => {
                            if headers{
                                let ref mut htype = &mut MeanGenomeCoverageEstimator::add_to_header(&mut htype);
                                print_header(htype)
                            }
                            let mut input_stream = coverm::genome::mosdepth_genome_coverage(
                            &bam_files,
                            separator,
                            limit_stream,
                            &mut MeanGenomeCoverageEstimator::new(min_fraction_covered),
                            print_zeros,
                            flag_filter,
                            single_genome);
                            let mut out = update_outputs(&mut output_stream, input_stream);
                            output_stream = out.clone();
                        },
                        "coverage_histogram" => {
                            limit_stream = true;
                            if headers{
                                let ref mut htype_hist = &mut PileupCountsGenomeCoverageEstimator::add_to_header(&mut htype_hist);
                                print_header(htype_hist)
                            }
                            coverm::genome::mosdepth_genome_coverage(
                            &bam_files,
                            separator,
                            limit_stream,
                            &mut PileupCountsGenomeCoverageEstimator::new(
                                min_fraction_covered),
                            print_zeros,
                            flag_filter,
                            single_genome);
                            // output_stream = out.clone();
                            },
                        "trimmed_mean" => {
                            if headers{
                                let ref mut htype = &mut TrimmedMeanGenomeCoverageEstimator::add_to_header(&mut htype);
                                print_header(htype)
                            }
                            let mut input_stream = coverm::genome::mosdepth_genome_coverage(
                                &bam_files,
                                separator,
                                limit_stream,
                                &mut get_trimmed_mean_estimator(m, min_fraction_covered),
                                print_zeros,
                                flag_filter,
                                single_genome);
                                let mut out = update_outputs(&mut output_stream, input_stream);
                                output_stream = out.clone();
                            },
                        "covered_fraction" => {
                            if headers{
                                let ref mut htype = &mut CoverageFractionGenomeCoverageEstimator::add_to_header(&mut htype);
                                print_header(htype)
                            }
                            let mut input_stream = coverm::genome::mosdepth_genome_coverage(
                            &bam_files,
                            separator,
                            limit_stream,
                            &mut CoverageFractionGenomeCoverageEstimator::new(
                                min_fraction_covered),
                            print_zeros,
                            flag_filter,
                            single_genome);
                            let mut out = update_outputs(&mut output_stream, input_stream);
                            output_stream = out.clone();
                        },
                        "variance" => {
                            if headers{
                                let ref mut htype = &mut VarianceGenomeCoverageEstimator::add_to_header(&mut htype);
                                print_header(htype)
                            }
                            let mut input_stream = coverm::genome::mosdepth_genome_coverage(
                            &bam_files,
                            separator,
                            limit_stream,
                            &mut VarianceGenomeCoverageEstimator::new(
                                min_fraction_covered),
                            print_zeros,
                            flag_filter,
                            single_genome);
                            let mut out = update_outputs(&mut output_stream, input_stream);
                            output_stream = out.clone();
                            },
                        _ => panic!("programming error")
                    }
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
                for cover in method{
                    match cover {

                        "mean" => {
                            if headers{
                                let ref mut htype = &mut MeanGenomeCoverageEstimator::add_to_header(&mut htype);
                                print_header(htype)
                            }
                            let mut input_stream = coverm::genome::mosdepth_genome_coverage_with_contig_names(
                            &bam_files,
                            &genomes_and_contigs,
                            limit_stream,
                            &mut MeanGenomeCoverageEstimator::new(min_fraction_covered),
                            print_zeros,
                            flag_filter);
                            let mut out = update_outputs(&mut output_stream, input_stream);
                            output_stream = out.clone();
                            },
                        "coverage_histogram" => {
                            limit_stream = true;
                            if headers{
                                let ref mut htype_hist = &mut PileupCountsGenomeCoverageEstimator::add_to_header(&mut htype_hist);
                                print_header(htype_hist)
                            }
                            coverm::genome::mosdepth_genome_coverage_with_contig_names(
                            &bam_files,
                            &genomes_and_contigs,
                            limit_stream,
                            &mut PileupCountsGenomeCoverageEstimator::new(
                                min_fraction_covered),
                            print_zeros,
                            flag_filter);
                            },
                        "trimmed_mean" => {
                            if headers{
                                let ref mut htype = &mut TrimmedMeanGenomeCoverageEstimator::add_to_header(&mut htype);
                                print_header(htype)
                            }
                            let mut input_stream = coverm::genome::mosdepth_genome_coverage_with_contig_names(
                                &bam_files,
                                &genomes_and_contigs,
                                limit_stream,
                                &mut get_trimmed_mean_estimator(m, min_fraction_covered),
                                print_zeros,
                                flag_filter);
                            let mut out = update_outputs(&mut output_stream, input_stream);
                            output_stream = out.clone();
                            },
                        "covered_fraction" => {
                            if headers{
                                let ref mut htype = &mut CoverageFractionGenomeCoverageEstimator::add_to_header(&mut htype);
                                print_header(htype)
                            }
                            let mut input_stream = coverm::genome::mosdepth_genome_coverage_with_contig_names(
                            &bam_files,
                            &genomes_and_contigs,
                            limit_stream,
                            &mut CoverageFractionGenomeCoverageEstimator::new(
                                min_fraction_covered),
                            print_zeros,
                            flag_filter);
                            let mut out = update_outputs(&mut output_stream, input_stream);
                            output_stream = out.clone();
                            },
                        "variance" => {
                            if headers{
                                let ref mut htype = &mut VarianceGenomeCoverageEstimator::add_to_header(&mut htype);
                                print_header(htype)
                            }
                            let mut input_stream = coverm::genome::mosdepth_genome_coverage_with_contig_names(
                            &bam_files,
                            &genomes_and_contigs,
                            limit_stream,
                            &mut VarianceGenomeCoverageEstimator::new(
                                min_fraction_covered),
                            print_zeros,
                            flag_filter);
                            let mut out = update_outputs(&mut output_stream, input_stream);
                            output_stream = out.clone();
                            },

                        _ => panic!("programming error")
                    }
                }
            }
        },
        Some("contig") => {
            let m = matches.subcommand_matches("contig").unwrap();
            set_log_level(m);
            let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
            let method: Vec<&str> = m.values_of("method").unwrap().collect();
            let min_fraction_covered = value_t!(m.value_of("min-covered-fraction"), f32).unwrap();
            let print_zeros = !m.is_present("no-zeros");
            let flag_filter = !m.is_present("no-flag-filter");
            let headers = !m.is_present("remove-headers");
            for cover in method {
                match cover {
                    "mean" => {
                        if headers{
                            let ref mut htype = &mut MeanGenomeCoverageEstimator::add_to_header(&mut htype);
                            print_header(htype)
                        }
                        let mut input_stream = coverm::contig::contig_coverage(
                        &bam_files,
                        limit_stream,
                        &mut MeanGenomeCoverageEstimator::new(min_fraction_covered),
                        print_zeros,
                        flag_filter);
                        let mut out = update_outputs(&mut output_stream, input_stream);
                        output_stream = out.clone();
                    },
                    "coverage_histogram" => {
                        limit_stream = true;
                        if headers{
                            let ref mut htype_hist = &mut PileupCountsGenomeCoverageEstimator::add_to_header(&mut htype_hist);
                            print_header(htype_hist)
                        }
                        coverm::contig::contig_coverage(
                        &bam_files,
                        limit_stream,
                        &mut PileupCountsGenomeCoverageEstimator::new(
                            min_fraction_covered),
                        print_zeros,
                        flag_filter);},
                    "trimmed_mean" => {
                        if headers{
                            let ref mut htype = &mut TrimmedMeanGenomeCoverageEstimator::add_to_header(&mut htype);
                            print_header(htype)
                        }
                        let mut input_stream = coverm::contig::contig_coverage(
                            &bam_files,
                            limit_stream,
                            &mut get_trimmed_mean_estimator(m, min_fraction_covered),
                            print_zeros,
                            flag_filter);
                        let mut out = update_outputs(&mut output_stream, input_stream);
                        output_stream = out.clone();},
                    "covered_fraction" => {
                        if headers{
                            let ref mut htype = &mut CoverageFractionGenomeCoverageEstimator::add_to_header(&mut htype);
                            print_header(htype)
                        }
                        let mut input_stream = coverm::contig::contig_coverage(
                        &bam_files,
                        limit_stream,
                        &mut CoverageFractionGenomeCoverageEstimator::new(
                            min_fraction_covered),
                        print_zeros,
                        flag_filter);
                        let mut out = update_outputs(&mut output_stream, input_stream);
                        output_stream = out.clone();
                        },
                    "variance" => {
                        if headers{
                            let ref mut htype = &mut VarianceGenomeCoverageEstimator::add_to_header(&mut htype);
                            print_header(htype)
                        }
                        let mut input_stream = coverm::contig::contig_coverage(
                        &bam_files,
                        limit_stream,
                        &mut VarianceGenomeCoverageEstimator::new(
                            min_fraction_covered),
                        print_zeros,
                        flag_filter);
                        let mut out = update_outputs(&mut output_stream, input_stream);
                        output_stream = out.clone();
                        },
                    _ => panic!("programming error")
                }
            }
        },
        _ => {
            app.print_help().unwrap();
            println!();
        }

    }
    if !limit_stream{
        print_header(htype);
        print_output_stream(output_stream);
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

Define the contigs in each genome (one of the following required):
   -s, --separator <CHARACTER>           This character separates genome names
                                         from contig names
   -f, --genome-fasta-files <PATH> ..    Path to FASTA files of each genome
   -d, --genome-fasta-directory <PATH>   Directory containing FASTA files of each
                                         genome
   --single-genome                       All contigs are from the same genome

Define mapping(s) (required):
   -b, --bam-files <PATH> ..             Path to reference-sorted BAM file(s)

Other arguments (optional):
   -m, --method METHOD                   Method for calculating coverage. One of:
                                              mean (default)
                                              trimmed_mean
                                              coverage_histogram
                                              covered_fraction
                                              variance
   --min-covered-fraction FRACTION       Genomes with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0.10]
   --trim-min FRACTION                   Remove this smallest fraction of positions
                                         when calculating trimmed_mean
                                         [default: 0.05]
   --trim-max FRACTION                   Maximum fraction for trimmed_mean
                                         calculations [default: 0.95]
   --no-zeros                            Omit printing of genomes that have zero
                                         coverage
   --no-flag-filter                      Do not ignore secondary and supplementary
                                         alignments, and improperly paired reads
   --headers                             Omit headers from output
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Ben J. Woodcroft <benjwoodcroft near gmail.com>
";

    let contig_help: &'static str =
        "coverm contig: Calculate read coverage per-contig

Define mapping(s) (required):
   -b, --bam-files <PATH> ..             Path to reference-sorted BAM file(s)

Other arguments (optional):
   -m, --method METHOD                   Method for calculating coverage. One of:
                                              mean (default)
                                              trimmed_mean
                                              coverage_histogram
                                              covered_fraction
                                              variance
   --min-covered-fraction FRACTION       Genomes with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0]
   --trim-min FRACTION                   Remove this smallest fraction of positions
                                         when calculating trimmed_mean
                                         [default: 0.05]
   --trim-max FRACTION                   Maximum fraction for trimmed_mean
                                         calculations [default: 0.95]
   --no-zeros                            Omit printing of genomes that have zero
                                         coverage
   --no-flag-filter                      Do not ignore secondary and supplementary
                                         alignments, and improperly paired reads
   --headers                             Omit headers from output
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Ben J. Woodcroft <benjwoodcroft near gmail.com>";

    return App::new("coverm")
        .version(crate_version!())
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
                     .conflicts_with("genome-fasta-files")
                     .conflicts_with("genome-fasta-directory")
                     .conflicts_with("single-genome")
                     .takes_value(true))
                .arg(Arg::with_name("genome-fasta-files")
                     .short("f")
                     .long("genome-fasta-files")
                     .multiple(true)
                     .conflicts_with("separator")
                     .conflicts_with("genome-fasta-directory")
                     .conflicts_with("single-genome")
                     .takes_value(true))
                .arg(Arg::with_name("genome-fasta-directory")
                     .short("d")
                     .long("genome-fasta-directory")
                     .conflicts_with("separator")
                     .conflicts_with("genome-fasta-files")
                     .conflicts_with("single-genome")
                     .takes_value(true))
                .arg(Arg::with_name("single-genome")
                     .long("single-genome")
                     .conflicts_with("separator")
                     .conflicts_with("genome-fasta-files")
                     .conflicts_with("genome-fasta-directory"))

                .arg(Arg::with_name("method")
                     .short("m")
                     .long("method")
                     .takes_value(true)
                     .multiple(true)
                     .possible_values(&[
                         "mean",
                         "trimmed_mean",
                         "coverage_histogram",
                         "covered_fraction",
                         "variance"])
                     .default_value("mean"))
                .arg(Arg::with_name("trim-min")
                     .long("trim-min")
                     .default_value("0.05"))
                .arg(Arg::with_name("trim-max")
                     .long("trim-max")
                     .default_value("0.95"))
                .arg(Arg::with_name("min-covered-fraction")
                     .long("min-covered-fraction")
                     .default_value("0.10"))
                .arg(Arg::with_name("no-zeros")
                     .long("no-zeros"))
                .arg(Arg::with_name("no-flag-filter")
                     .long("no-flag-filter"))
                .arg(Arg::with_name("remove-headers")
                     .long("remove-headers"))

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
                     .multiple(true)
                     .possible_values(&[
                         "mean",
                         "trimmed_mean",
                         "coverage_histogram",
                         "covered_fraction",
                         "variance"])
                     .default_value("mean"))
                .arg(Arg::with_name("min-covered-fraction")
                     .long("min-covered-fraction")
                     .default_value("0.0"))
                .arg(Arg::with_name("trim-min")
                     .long("trim-min")
                     .default_value("0.05"))
                .arg(Arg::with_name("trim-max")
                     .long("trim-max")
                     .default_value("0.95"))
                .arg(Arg::with_name("no-zeros")
                     .long("no-zeros"))
                .arg(Arg::with_name("no-flag-filter")
                     .long("no-flag-filter"))
                .arg(Arg::with_name("remove-headers")
                     .long("remove-headers"))
                .arg(Arg::with_name("verbose")
                     .short("v")
                     .long("verbose"))
                .arg(Arg::with_name("quiet")
                     .short("q")
                     .long("quiet")));
}
