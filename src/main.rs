extern crate coverm;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::bam_generator::*;
use coverm::filter;
use coverm::external_command_checker;

extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;

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

    match matches.subcommand_name() {

        Some("genome") => {
            let m = matches.subcommand_matches("genome").unwrap();
            set_log_level(m);
            let filtering = doing_filtering(m);
            debug!("Doing filtering? {}", filtering);

            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
                if filtering {
                    let filter_params = FilterParameters::generate_from_clap(m);
                    run_genome(
                        coverm::bam_generator::generate_filtered_bam_readers_from_bam_files(
                            bam_files,
                            filter_params.min_aligned_length,
                            filter_params.min_percent_identity),
                        m);
                } else {
                    run_genome(
                        coverm::bam_generator::generate_named_bam_readers_from_bam_files(
                            bam_files),
                        m);
                }
            } else {
                external_command_checker::check_for_bwa();
                external_command_checker::check_for_samtools();
                if filtering {
                    debug!("Mapping and filtering..");
                    run_genome(
                        get_streamed_filtered_bam_readers(m),
                        m);
                } else {
                    run_genome(
                        get_streamed_bam_readers(m),
                        m);
                }
            }
        },
        Some("contig") => {
            let m = matches.subcommand_matches("contig").unwrap();
            set_log_level(m);
            let method = m.value_of("method").unwrap();
            let min_fraction_covered = value_t!(m.value_of("min-covered-fraction"), f32).unwrap();
            let print_zeros = !m.is_present("no-zeros");
            let flag_filter = !m.is_present("no-flag-filter");
            let filtering = doing_filtering(m);

            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
                if filtering {
                    let filter_params = FilterParameters::generate_from_clap(m);
                    let mut bam_readers = coverm::bam_generator::generate_filtered_bam_readers_from_bam_files(
                        bam_files,
                        filter_params.min_aligned_length,
                        filter_params.min_percent_identity);
                    run_contig(method, bam_readers, min_fraction_covered, print_zeros, flag_filter, m);
                } else {
                    let mut bam_readers = coverm::bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    run_contig(method, bam_readers, min_fraction_covered, print_zeros, flag_filter, m);
                }
            } else {
                external_command_checker::check_for_bwa();
                external_command_checker::check_for_samtools();
                if filtering {
                    debug!("Filtering..");
                    run_contig(method,
                               get_streamed_filtered_bam_readers(m),
                               min_fraction_covered, print_zeros, flag_filter, m);
                } else {
                    debug!("Not filtering..");
                    run_contig(method,
                               get_streamed_bam_readers(m),
                               min_fraction_covered, print_zeros, flag_filter, m);
                }
            }
        },
        Some("filter") => {
            let m = matches.subcommand_matches("filter").unwrap();
            set_log_level(m);

            let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
            let output_bam_files: Vec<&str> = m.values_of("output-bam-files").unwrap().collect();
            if bam_files.len() != output_bam_files.len() {
                panic!("The number of input BAM files must be the same as the number output")
            }

            let min_aligned_length = value_t!(m.value_of("min-aligned-length"), u32).unwrap();
            let min_percent_identity = value_t!(m.value_of("min-percent-identity"), f32).unwrap();

            let num_threads = value_t!(m.value_of("threads"), u16).unwrap();


            for (bam, output) in bam_files.iter().zip(output_bam_files.iter()) {
                let mut reader = bam::Reader::from_path(bam).expect(
                    &format!("Unable to find BAM file {}", bam));
                let header = bam::header::Header::from_template(reader.header());
                let mut writer = bam::Writer::from_path(
                    output,
                    &header
                ).expect(&format!("Failed to write BAM file {}", output));
                writer.set_threads(num_threads as usize).expect("Failed to set num threads in writer");
                let mut filtered = filter::ReferenceSortedBamFilter::new(
                    reader, min_aligned_length, min_percent_identity);

                let mut record = bam::record::Record::new();
                while filtered.read(&mut record).is_ok() {
                    debug!("Writing.. {:?}", record.qname());
                    writer.write(&record).expect("Failed to write BAM record");
                }
            }
        },
        _ => {
            app.print_help().unwrap();
            println!();
        }
    }
}

fn doing_filtering(m: &clap::ArgMatches) -> bool {
    m.is_present("min-aligned-length") ||
        m.is_present("min-percent-identity")
}

fn run_genome<R: coverm::bam_generator::NamedBamReader,
              T: coverm::bam_generator::NamedBamReaderGenerator<R>>(
    bam_generators: Vec<T>,
    m: &clap::ArgMatches) {

    let method = m.value_of("method").unwrap();
    let min_fraction_covered = value_t!(m.value_of("min-covered-fraction"), f32).unwrap();
    if min_fraction_covered > 1.0 || min_fraction_covered < 0.0 {
        eprintln!("Minimum fraction covered parameter cannot be < 0 or > 1, found {}", min_fraction_covered);
        process::exit(1)
    }
    let print_zeros = !m.is_present("no-zeros");
    let flag_filter = !m.is_present("no-flag-filter");
    let single_genome = m.is_present("single-genome");

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

        match method {
            "mean" => coverm::genome::mosdepth_genome_coverage(
                bam_generators,
                separator,
                &mut std::io::stdout(),
                &mut MeanGenomeCoverageEstimator::new(min_fraction_covered),
                print_zeros,
                flag_filter,
                single_genome),
            "coverage_histogram" => coverm::genome::mosdepth_genome_coverage(
                bam_generators,
                separator,
                &mut std::io::stdout(),
                &mut PileupCountsGenomeCoverageEstimator::new(
                    min_fraction_covered),
                print_zeros,
                flag_filter,
                single_genome),
            "trimmed_mean" => {
                coverm::genome::mosdepth_genome_coverage(
                    bam_generators,
                    separator,
                    &mut std::io::stdout(),
                    &mut get_trimmed_mean_estimator(m, min_fraction_covered),
                    print_zeros,
                    flag_filter,
                    single_genome)},
            "covered_fraction" => coverm::genome::mosdepth_genome_coverage(
                bam_generators,
                separator,
                &mut std::io::stdout(),
                &mut CoverageFractionGenomeCoverageEstimator::new(
                    min_fraction_covered),
                print_zeros,
                flag_filter,
                single_genome),
            "variance" => coverm::genome::mosdepth_genome_coverage(
                bam_generators,
                separator,
                &mut std::io::stdout(),
                &mut VarianceGenomeCoverageEstimator::new(
                    min_fraction_covered),
                print_zeros,
                flag_filter,
                single_genome),
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
                bam_generators,
                &genomes_and_contigs,
                &mut std::io::stdout(),
                &mut MeanGenomeCoverageEstimator::new(min_fraction_covered),
                print_zeros,
                flag_filter),
            "coverage_histogram" => coverm::genome::mosdepth_genome_coverage_with_contig_names(
                bam_generators,
                &genomes_and_contigs,
                &mut std::io::stdout(),
                &mut PileupCountsGenomeCoverageEstimator::new(
                    min_fraction_covered),
                print_zeros,
                flag_filter),
            "trimmed_mean" => {
                coverm::genome::mosdepth_genome_coverage_with_contig_names(
                    bam_generators,
                    &genomes_and_contigs,
                    &mut std::io::stdout(),
                    &mut get_trimmed_mean_estimator(m, min_fraction_covered),
                    print_zeros,
                    flag_filter)},
            "covered_fraction" => coverm::genome::mosdepth_genome_coverage_with_contig_names(
                bam_generators,
                &genomes_and_contigs,
                &mut std::io::stdout(),
                &mut CoverageFractionGenomeCoverageEstimator::new(
                    min_fraction_covered),
                print_zeros,
                flag_filter),
            "variance" => coverm::genome::mosdepth_genome_coverage_with_contig_names(
                bam_generators,
                &genomes_and_contigs,
                &mut std::io::stdout(),
                &mut VarianceGenomeCoverageEstimator::new(
                    min_fraction_covered),
                print_zeros,
                flag_filter),

            _ => panic!("programming error")
        }
    }
}

struct FilterParameters {
    min_aligned_length: u32,
    min_percent_identity: f32
}
impl FilterParameters {
    pub fn generate_from_clap(m: &clap::ArgMatches) -> FilterParameters {
        FilterParameters {
            min_aligned_length: match m.is_present("min-aligned-length") {
                true => value_t!(m.value_of("min-aligned-length"), u32).unwrap(),
                false => 0
            },
            min_percent_identity: match m.is_present("min-percent-identity") {
                true => value_t!(m.value_of("min-percent-identity"), f32).unwrap(),
                false => 0.0
            }
        }
    }
}

fn get_streamed_bam_readers(m: &clap::ArgMatches) -> Vec<StreamingNamedBamReaderGenerator> {
    let reference = m.value_of("reference").unwrap();
    let read1: Vec<&str> = m.values_of("read1").unwrap().collect();
    let read2: Vec<&str> = m.values_of("read2").unwrap().collect();
    let threads: u16 = m.value_of("threads").unwrap().parse::<u16>()
        .expect("Failed to convert threads argument into integer");
    if read1.len() != read2.len() {
        panic!("The number of forward read files ({}) was not the same as the number of reverse read files ({})",
               read1.len(), read2.len())
    }
    let mut bam_readers = vec![];
    coverm::bam_generator::generate_bwa_index(&reference);
    for (i, _) in read1.iter().enumerate() {
        bam_readers.push(
            coverm::bam_generator::generate_named_bam_readers_from_read_couple(
                reference, read1[i], read2[i], threads));
        debug!("Back");
    }
    debug!("Finished BAM setup");
    return bam_readers
}

fn get_streamed_filtered_bam_readers(m: &clap::ArgMatches) -> Vec<StreamingFilteredNamedBamReaderGenerator> {
    let reference = m.value_of("reference").unwrap();
    let read1: Vec<&str> = m.values_of("read1").unwrap().collect();
    let read2: Vec<&str> = m.values_of("read2").unwrap().collect();
    let threads: u16 = m.value_of("threads").unwrap().parse::<u16>()
        .expect("Failed to convert threads argument into integer");
    if read1.len() != read2.len() {
        panic!("The number of forward read files ({}) was not the same as the number of reverse read files ({})",
               read1.len(), read2.len())
    }
    let mut bam_readers = vec![];
    let filter_params = FilterParameters::generate_from_clap(m);
    coverm::bam_generator::generate_bwa_index(&reference);
    for (i, _) in read1.iter().enumerate() {
        bam_readers.push(
            coverm::bam_generator::generate_filtered_named_bam_readers_from_read_couple(
                reference, read1[i], read2[i], threads,
                filter_params.min_aligned_length,
                filter_params.min_percent_identity
            ));
        debug!("Back");
    }
    debug!("Finished BAM setup");
    return bam_readers
}

fn get_trimmed_mean_estimator(
    m: &clap::ArgMatches,
    min_fraction_covered: f32) -> TrimmedMeanGenomeCoverageEstimator {
    let min = value_t!(m.value_of("trim-min"), f32).unwrap();
    let max = value_t!(m.value_of("trim-max"), f32).unwrap();
    if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
        eprintln!("error: Trim bounds must be between 0 and 1, and min must be less than max, found {} and {}", min, max);
        process::exit(1)
    }
    TrimmedMeanGenomeCoverageEstimator::new(
        min, max, min_fraction_covered)
}


fn run_contig<R: coverm::bam_generator::NamedBamReader,
        T: coverm::bam_generator::NamedBamReaderGenerator<R>>(
    method: &str,
    bam_readers: Vec<T>,
    min_fraction_covered: f32,
    print_zeros: bool,
    flag_filter: bool,
    m: &clap::ArgMatches) {
    match method {
        "mean" => coverm::contig::contig_coverage(
            bam_readers,
            &mut std::io::stdout(),
            &mut MeanGenomeCoverageEstimator::new(min_fraction_covered),
            print_zeros,
            flag_filter),
        "coverage_histogram" => coverm::contig::contig_coverage(
            bam_readers,
            &mut std::io::stdout(),
            &mut PileupCountsGenomeCoverageEstimator::new(
                min_fraction_covered),
            print_zeros,
            flag_filter),
        "trimmed_mean" => {
            coverm::contig::contig_coverage(
                bam_readers,
                &mut std::io::stdout(),
                &mut get_trimmed_mean_estimator(m, min_fraction_covered),
                print_zeros,
                flag_filter)},
        "covered_fraction" => coverm::contig::contig_coverage(
            bam_readers,
            &mut std::io::stdout(),
            &mut CoverageFractionGenomeCoverageEstimator::new(
                min_fraction_covered),
            print_zeros,
            flag_filter),
        "variance" => coverm::contig::contig_coverage(
            bam_readers,
            &mut std::io::stdout(),
            &mut VarianceGenomeCoverageEstimator::new(
                min_fraction_covered),
            print_zeros,
            flag_filter),
        _ => panic!("programming error")
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

Define the contigs in each genome (exactly one of the following is required):
   -s, --separator <CHARACTER>           This character separates genome names
                                         from contig names
   -f, --genome-fasta-files <PATH> ..    Path to FASTA files of each genome
   -d, --genome-fasta-directory <PATH>   Directory containing FASTA files of each
                                         genome
   --single-genome                       All contigs are from the same genome

Define mapping(s) (required):
  Either define BAM:
    -b, --bam-files <PATH> ..            Path to reference-sorted BAM file(s)

  Or do mapping:
   -1 <PATH> ..                          Forward FASTA/Q files for mapping
   -2 <PATH> ..                          Reverse FASTA/Q files for mapping
   -r, --reference <PATH>                BWA indexed FASTA file of contigs
   -t, --threads <INT>                   Number of threads to use for mapping

Alignment filtering (optional):
   --min-aligned-length <INT>            Exclude pairs with smaller numbers of
                                         aligned bases [default: 0]
   --min-percent-identity <FLOAT>        Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95% [default 0.0]

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
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Ben J. Woodcroft <benjwoodcroft near gmail.com>
";

    let contig_help: &'static str =
        "coverm contig: Calculate read coverage per-contig

Define mapping(s) (required):
  Either define BAM:
    -b, --bam-files <PATH> ..            Path to reference-sorted BAM file(s)

  Or do mapping:
   -1 <PATH> ..                          Forward FASTA/Q files for mapping
   -2 <PATH> ..                          Reverse FASTA/Q files for mapping
   -r, --reference <PATH>                BWA indexed FASTA file of contigs
   -t, --threads <INT>                   Number of threads to use for mapping

Alignment filtering (optional):
   --min-aligned-length <INT>            Exclude pairs with smaller numbers of
                                         aligned bases [default: none]
   --min-percent-identity <FLOAT>        Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95% [default: none]

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
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Ben J. Woodcroft <benjwoodcroft near gmail.com>";

    let filter_help: &'static str =
        "coverm filter: Remove alignments with insufficient identity.

Only primary, non-supplementary alignments are considered, and output files
are grouped by reference, but not sorted by position.

Files (both required):
   -b, --bam-files <PATH> ..           Path to reference-sorted BAM file(s)
   -o, --output-bam-files <PATH> ..    Path to corresponding output file(s)

Thresholds:
   --min-aligned-length <INT>          Exclude pairs with smaller numbers of
                                       aligned bases [default: 0]
   --min-percent-identity <FLOAT>      Exclude pairs by overall percent
                                       identity e.g. 0.95 for 95% [default 0.0]

Other:
   -t, --threads <INT>                 Number of threads for output compression
                                       [default 1]
   -v, --verbose                       Print extra debugging information
   -q, --quiet                         Unless there is an error, do not print
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

Main modes:
\tcontig\tCalculate coverage of contigs
\tgenome\tCalculate coverage of genomes

Utilities:
\tfilter\tRemove alignments with insufficient identity

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
                     .takes_value(true))
                .arg(Arg::with_name("read1")
                     .short("-1")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless("bam-files")
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("read2")
                     .short("-2")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless("bam-files")
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("reference")
                     .short("-r")
                     .long("reference")
                     .takes_value(true)
                     .required_unless("bam-files")
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("threads")
                     .short("-t")
                     .long("threads")
                     .default_value("1")
                     .takes_value(true))

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

                .arg(Arg::with_name("min-aligned-length")
                     .long("min-aligned-length")
                     .conflicts_with("no-flag-filter"))
                .arg(Arg::with_name("min-percent-identity")
                     .long("min-percent-identity")
                     .conflicts_with("no-flag-filter"))

                .arg(Arg::with_name("method")
                     .short("m")
                     .long("method")
                     .takes_value(true)
                     .possible_values(&[
                         "mean",
                         "trimmed_mean",
                         "coverage_histogram",
                         "covered_fraction"])
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
                     .required_unless("read1"))
                .arg(Arg::with_name("read1")
                     .short("-1")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless("bam-files")
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("read2")
                     .short("-2")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless("bam-files")
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("reference")
                     .short("-r")
                     .long("reference")
                     .takes_value(true)
                     .required_unless("bam-files")
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("threads")
                     .short("-t")
                     .long("threads")
                     .default_value("1")
                     .takes_value(true))

                .arg(Arg::with_name("min-aligned-length")
                     .long("min-aligned-length")
                     .takes_value(true)
                     .conflicts_with("no-flag-filter"))
                .arg(Arg::with_name("min-percent-identity")
                     .long("min-percent-identity")
                     .takes_value(true)
                     .conflicts_with("no-flag-filter"))

                .arg(Arg::with_name("method")
                     .short("m")
                     .long("method")
                     .takes_value(true)
                     .possible_values(&[
                         "mean",
                         "trimmed_mean",
                         "coverage_histogram",
                         "covered_fraction"])
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
                .arg(Arg::with_name("verbose")
                     .short("v")
                     .long("verbose"))
                .arg(Arg::with_name("quiet")
                     .short("q")
                     .long("quiet")))
        .subcommand(
            SubCommand::with_name("filter")
                .about("Remove alignments with insufficient identity")
                .help(filter_help)

                .arg(Arg::with_name("bam-files")
                     .short("b")
                     .long("bam-files")
                     .multiple(true)
                     .takes_value(true)
                     .required(true))
                .arg(Arg::with_name("output-bam-files")
                     .short("o")
                     .long("output-bam-files")
                     .multiple(true)
                     .takes_value(true)
                     .required(true))

                .arg(Arg::with_name("min-aligned-length")
                     .long("min-aligned-length")
                     .default_value("0"))
                .arg(Arg::with_name("min-percent-identity")
                     .long("min-percent-identity")
                     .default_value("0.0"))
                .arg(Arg::with_name("threads")
                     .long("threads")
                     .short("t")
                     .default_value("1"))

                .arg(Arg::with_name("verbose")
                     .short("v")
                     .long("verbose"))
                .arg(Arg::with_name("quiet")
                     .short("q")
                     .long("quiet")))
}
