extern crate coverm;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::bam_generator::*;
use coverm::filter;
use coverm::external_command_checker;
use coverm::coverage_takers::*;
use coverm::mapping_parameters::*;

extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;

use std::env;
use std::str;
use std::process;
use std::io::Write;

extern crate clap;
use clap::*;

#[macro_use]
extern crate log;
extern crate env_logger;
use log::LevelFilter;
use env_logger::Builder;

extern crate tempfile;


fn main(){
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    let print_stream = &mut std::io::stdout();
    set_log_level(&matches, false);

    match matches.subcommand_name() {

        Some("genome") => {
            let m = matches.subcommand_matches("genome").unwrap();
            set_log_level(m, true);
            let filtering = doing_filtering(m);
            debug!("Doing filtering? {}", filtering);

            let mut stream = std::io::stdout();
            let mut estimators_and_taker = EstimatorsAndTaker::generate_from_clap(
                m, &mut stream);
            write!(print_stream, "Sample\tGenome");
            for (i, estimator) in estimators_and_taker.estimators.iter().enumerate() {
                for header in estimator.column_headers() {
                    if estimators_and_taker.columns_to_normalise.contains(&i) {
                        write!(print_stream, "\tRelative Abundance (%)");
                    } else {
                        write!(print_stream, "\t{}", header);
                    }
                }
            }
            writeln!(print_stream);

            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
                if filtering {
                    let filter_params = FilterParameters::generate_from_clap(m);
                    run_genome(
                        coverm::bam_generator::generate_filtered_bam_readers_from_bam_files(
                            bam_files,
                            filter_params.min_aligned_length,
                            filter_params.min_percent_identity,
                            filter_params.min_aligned_percent),
                        m,
                        &mut estimators_and_taker);
                } else {
                    run_genome(
                        coverm::bam_generator::generate_named_bam_readers_from_bam_files(
                            bam_files),
                        m,
                        &mut estimators_and_taker);
                }
            } else {
                external_command_checker::check_for_bwa();
                external_command_checker::check_for_samtools();
                if filtering {
                    debug!("Mapping and filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(m);
                    let mut all_generators = vec!();
                    for set in generator_sets {
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    run_genome(
                        all_generators,
                        m,
                        &mut estimators_and_taker);
                } else {
                    let generator_sets = get_streamed_bam_readers(m);
                    let mut all_generators = vec!();
                    for set in generator_sets {
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    run_genome(
                        all_generators,
                        m,
                        &mut estimators_and_taker);
                };
            }
        },
        Some("contig") => {
            let m = matches.subcommand_matches("contig").unwrap();
            set_log_level(m, true);
            let print_zeros = !m.is_present("no-zeros");
            let flag_filter = !m.is_present("no-flag-filter");
            let filtering = doing_filtering(m);

            let mut stream = std::io::stdout();
            let mut estimators_and_taker = EstimatorsAndTaker::generate_from_clap(
                m, &mut stream);
            write!(print_stream, "Sample\tContig");
            for estimator in &(estimators_and_taker.estimators) {
                for ref header in estimator.column_headers() {
                    write!(print_stream, "\t{}", header);
                }
            }
            writeln!(print_stream);

            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
                if filtering {
                    let filter_params = FilterParameters::generate_from_clap(m);
                    let mut bam_readers = coverm::bam_generator::generate_filtered_bam_readers_from_bam_files(
                        bam_files,
                        filter_params.min_aligned_length,
                        filter_params.min_percent_identity,
                        filter_params.min_aligned_percent);
                    run_contig(
                        &mut estimators_and_taker.estimators,
                        bam_readers,
                        print_zeros,
                        flag_filter,
                        &mut estimators_and_taker.taker);
                } else {
                    let mut bam_readers = coverm::bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    run_contig(
                        &mut estimators_and_taker.estimators,
                        bam_readers,
                        print_zeros,
                        flag_filter,
                        &mut estimators_and_taker.taker);
                }
            } else {
                external_command_checker::check_for_bwa();
                external_command_checker::check_for_samtools();
                if filtering {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(m);
                    for generator_set in generator_sets {
                        run_contig(
                            &mut estimators_and_taker.estimators,
                            generator_set.generators,
                            print_zeros,
                            flag_filter,
                            &mut estimators_and_taker.taker);
                    }
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m);
                    for generator_set in generator_sets {
                        run_contig(
                            &mut estimators_and_taker.estimators,
                            generator_set.generators,
                            print_zeros,
                            flag_filter,
                            &mut estimators_and_taker.taker);
                    }
                }
            }
        },
        Some("filter") => {
            let m = matches.subcommand_matches("filter").unwrap();
            set_log_level(m, true);

            let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
            let output_bam_files: Vec<&str> = m.values_of("output-bam-files").unwrap().collect();
            if bam_files.len() != output_bam_files.len() {
                panic!("The number of input BAM files must be the same as the number output")
            }

            let min_aligned_length = value_t!(m.value_of("min-aligned-length"), u32).unwrap();
            let min_percent_identity = value_t!(m.value_of("min-percent-identity"), f32).unwrap();
            let min_aligned_percent = value_t!(m.value_of("min-aligned-percent"), f32).unwrap();

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
                    reader, min_aligned_length, min_percent_identity, min_aligned_percent);

                let mut record = bam::record::Record::new();
                while filtered.read(&mut record).is_ok() {
                    debug!("Writing.. {:?}", record.qname());
                    writer.write(&record).expect("Failed to write BAM record");
                }
            }
        },
        Some("make") => {
            let m = matches.subcommand_matches("make").unwrap();
            set_log_level(m, true);
            external_command_checker::check_for_bwa();
            external_command_checker::check_for_samtools();

            let output_directory = m.value_of("output-directory").unwrap();
            let params = MappingParameters::generate_from_clap(&m);
            let mut generator_sets = vec!();

            for reference_wise_params in params {
                let mut bam_readers = vec![];
                let index = coverm::bwa_index_maintenance::generate_bwa_index(
                    reference_wise_params.reference);

                for p in reference_wise_params {
                    bam_readers.push(
                        coverm::bam_generator::generate_bam_maker_generator_from_reads(
                            index.index_path(),
                            p.read1,
                            p.read2,
                            p.read_format.clone(),
                            p.threads,
                            &generate_cached_bam_file_name(
                                output_directory, p.reference, p.read1)));
                }

                debug!("Finished BAM setup");
                let to_return = BamGeneratorSet {
                    generators: bam_readers,
                    index: index
                };
                generator_sets.push(to_return);
            }

            let mut i = 1;
            for generator_set in generator_sets {
                for generator in generator_set.generators {
                    info!("Running mapping number {} ..", i);
                    generator.start().finish();
                    i += 1;
                }
            }
        },
        _ => {
            app.print_help().unwrap();
            println!();
        }
    }
}

struct EstimatorsAndTaker<'a> {
    estimators: Vec<CoverageEstimator>,
    taker: CoverageTakerType<'a>,
    columns_to_normalise: Vec<usize>,
}

impl<'a> EstimatorsAndTaker<'a> {
    pub fn generate_from_clap(
        m: &clap::ArgMatches, stream: &'a mut std::io::Stdout)
        -> EstimatorsAndTaker<'a> {
        let mut estimators = vec!();
        let min_fraction_covered = value_t!(m.value_of("min-covered-fraction"), f32).unwrap();

        if min_fraction_covered > 1.0 || min_fraction_covered < 0.0 {
            eprintln!("Minimum fraction covered parameter cannot be < 0 or > 1, found {}", min_fraction_covered);
            process::exit(1)
        }
        let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u32).unwrap();

        let methods: Vec<&str> = m.values_of("methods").unwrap().collect();
        let mut columns_to_normalise: Vec<usize> = vec!();

        for (i, method) in methods.iter().enumerate() {
            let estimator = match method {
                &"mean" => {
                    CoverageEstimator::new_estimator_mean(
                        min_fraction_covered, contig_end_exclusion)
                },
                &"coverage_histogram" => {
                    CoverageEstimator::new_estimator_pileup_counts(
                        min_fraction_covered, contig_end_exclusion)
                },
                &"trimmed_mean" => {
                    let min = value_t!(m.value_of("trim-min"), f32).unwrap();
                    let max = value_t!(m.value_of("trim-max"), f32).unwrap();
                    if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                        eprintln!("error: Trim bounds must be between 0 and 1, and min must be less than max, found {} and {}", min, max);
                        process::exit(1)
                    }
                    CoverageEstimator::new_estimator_trimmed_mean(
                        min, max, min_fraction_covered, contig_end_exclusion)
                },
                &"covered_fraction" => {
                    CoverageEstimator::new_estimator_covered_fraction(
                        min_fraction_covered, contig_end_exclusion)
                },
                &"variance" =>{
                    CoverageEstimator::new_estimator_variance(
                        min_fraction_covered, contig_end_exclusion)
                },
                &"length" =>{
                    CoverageEstimator::new_estimator_length()
                },
                &"relative_abundance" => {
                    columns_to_normalise.push(i);
                    CoverageEstimator::new_estimator_mean(
                        min_fraction_covered, contig_end_exclusion)
                }
                _ => panic!("programming error")
            };
            estimators.push(estimator);
        }

        let taker;
        if methods.contains(&"coverage_histogram") {
            if methods.len() > 1 {
                panic!("Cannot specify the coverage_histogram method with any other coverage methods")
            } else {
                debug!("Coverage histogram type coverage taker being used");
                taker = CoverageTakerType::new_pileup_coverage_coverage_printer(stream)
            }
        } else if columns_to_normalise.len() == 0 {
            debug!("Streaming regular coverage output");
            taker = CoverageTakerType::new_single_float_coverage_streaming_coverage_printer(
                stream)
        } else {
            debug!("Cached regular coverage taker from columns: {:?}",
                   columns_to_normalise);
            taker = CoverageTakerType::new_cached_single_float_coverage_taker(
                estimators.len())
        }

        return EstimatorsAndTaker {
            estimators: estimators,
            taker: taker,
            columns_to_normalise: columns_to_normalise
        }
    }
}


fn doing_filtering(m: &clap::ArgMatches) -> bool {
    m.is_present("min-aligned-length") ||
        m.is_present("min-percent-identity") ||
        m.is_present("min-aligned-percent")
}

fn run_genome<'a,
              R: coverm::bam_generator::NamedBamReader,
              T: coverm::bam_generator::NamedBamReaderGenerator<R>>(
    bam_generators: Vec<T>,
    m: &clap::ArgMatches,
    estimators_and_taker: &'a mut EstimatorsAndTaker<'a>) {

    let print_zeros = !m.is_present("no-zeros");
    let flag_filter = !m.is_present("no-flag-filter");
    let single_genome = m.is_present("single-genome");
    let reads_mapped = match m.is_present("separator") || single_genome {
        true => {
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
            coverm::genome::mosdepth_genome_coverage(
                bam_generators,
                separator,
	              &mut estimators_and_taker.taker,
                print_zeros,
                &mut estimators_and_taker.estimators,
                flag_filter,
                single_genome)
        },

        false => {
            let genomes_and_contigs;
            if m.is_present("genome-fasta-files"){
                let genome_fasta_files: Vec<&str> = m.values_of("genome-fasta-files").unwrap().collect();
                genomes_and_contigs = coverm::read_genome_fasta_files(&genome_fasta_files);
            } else if m.is_present("genome-fasta-directory") {
                let dir = m.value_of("genome-fasta-directory").unwrap();
                let paths = std::fs::read_dir(dir).unwrap();
                let mut genome_fasta_files: Vec<String> = vec!();
                let extension = m.value_of("genome-fasta-extension").unwrap();
                for path in paths {
                    let file = path.unwrap().path();
                    match file.extension() {
                        Some(ext) => {
                            if ext == extension {
                                let s = String::from(file.to_string_lossy());
                                genome_fasta_files.push(s);
                            } else {
                                info!(
                                    "Not using directory entry '{}' as a genome FASTA file, as it does not end with the extension '{}'",
                                    file.to_str().expect("UTF8 error in filename"),
                                    extension);
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
                if strs.len() == 0 {
                    panic!("Found 0 genomes from the genome-fasta-directory, cannot continue.")
                }
                info!("Calculating coverage for {} genomes ..", strs.len());
                genomes_and_contigs = coverm::read_genome_fasta_files(&strs);
            } else {
                panic!("Either a separator (-s) or path(s) to genome FASTA files (with -d or -f) must be given");
            }
            coverm::genome::mosdepth_genome_coverage_with_contig_names(
                bam_generators,
                &genomes_and_contigs,
                &mut estimators_and_taker.taker,
                print_zeros,
                flag_filter,
                &mut estimators_and_taker.estimators)
        }
    };

    match estimators_and_taker.taker {
        CoverageTakerType::CachedSingleFloatCoverageTaker{..} => {
            coverm::coverage_takers::print_sparse_cached_coverage_taker(
                &estimators_and_taker.taker, &mut std::io::stdout(), &reads_mapped,
                &estimators_and_taker.columns_to_normalise);
        }, _ => {}
    }
}

struct FilterParameters {
    min_aligned_length: u32,
    min_percent_identity: f32,
    min_aligned_percent: f32
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
            },
            min_aligned_percent: match m.is_present("min-aligned-percent") {
                true => value_t!(m.value_of("min-aligned-percent"), f32).unwrap(),
                false => 0.0
            }
        }
    }
}



fn get_streamed_bam_readers<'a>(
    m: &clap::ArgMatches) -> Vec<BamGeneratorSet<StreamingNamedBamReaderGenerator>> {

    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }

    let params = MappingParameters::generate_from_clap(&m);
    let mut generator_set = vec!();
    for reference_wise_params in params {
        let mut bam_readers = vec![];
        let index = coverm::bwa_index_maintenance::generate_bwa_index(
            reference_wise_params.reference);

        let reference = reference_wise_params.reference;
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.is_present("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.value_of("bam-file-cache-directory").unwrap(),
                        reference, naming_readset);
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        for p in reference_wise_params {
            bam_readers.push(
                coverm::bam_generator::generate_named_bam_readers_from_reads(
                    index.index_path(),
                    p.read1,
                    p.read2,
                    p.read_format.clone(),
                    p.threads,
                    bam_file_cache(p.read1).as_ref().map(String::as_ref)));
        }

        debug!("Finished BAM setup");
        let to_return = BamGeneratorSet {
            generators: bam_readers,
            index: index
        };
        generator_set.push(to_return);
    };
    return generator_set;
}

fn generate_cached_bam_file_name(directory: &str, reference: &str, read1_path: &str) -> String {
    std::path::Path::new(directory).to_str()
        .expect("Unable to covert bam-file-cache-directory name into str").to_string()+"/"+
        &std::path::Path::new(reference).file_name()
        .expect("Unable to convert reference to file name").to_str()
        .expect("Unable to covert file name into str").to_string()+"."+
        &std::path::Path::new(read1_path).file_name()
        .expect("Unable to convert read1 name to file name").to_str()
        .expect("Unable to covert file name into str").to_string()+".bam"
}

fn setup_bam_cache_directory(cache_directory: &str) {
    let path = std::path::Path::new(cache_directory);
    if path.is_dir() {
        if path.metadata().expect("Unable to read metadata for cache directory")
            .permissions().readonly() {
                panic!("Cache directory {} does not appear to be writeable, not continuing",
                       cache_directory);
            } else {
                info!("Writing BAM files to already existing directory {}", cache_directory)
            }
    } else {
        match path.parent() {
            Some(parent) => {
                println!("parent: {:?}",parent);
                let parent2 = match parent == std::path::Path::new("") {
                    true => std::path::Path::new("."),
                    false => parent
                };
                if parent2.canonicalize().expect(
                    &format!("Unable to canonicalize parent of cache directory {}", cache_directory)).is_dir() {
                    if parent2.metadata().expect(
                        &format!("Unable to get metadata for parent of cache directory {}",
                                 cache_directory))
                        .permissions().readonly() {
                            panic!(
                                "The parent directory of the (currently non-existent) \
                                 cache directory {} is not writeable, not continuing",
                                cache_directory);
                        } else {
                            info!("Creating cache directory {}", cache_directory);
                            std::fs::create_dir(path).expect("Unable to create cache directory");
                        }
                } else {
                    panic!("The parent directory of the cache directory {} does not \
                            yet exist, so not creating that cache directory, and not continuing.",
                           cache_directory)
                }
            },
            None => {
                panic!("Cannot create root directory {}", cache_directory)
            }
        }
    }
    // Test writing a tempfile to the directory, to test it actually is
    // writeable.
    let tf_result = tempfile::tempfile_in(path);
    if tf_result.is_err() {
        panic!("Failed to create test file in bam cache directory: {}", tf_result.err().unwrap())
    }
}

fn get_streamed_filtered_bam_readers(
    m: &clap::ArgMatches) -> Vec<BamGeneratorSet<StreamingFilteredNamedBamReaderGenerator>> {

    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }

    let params = MappingParameters::generate_from_clap(&m);
    let mut generator_set = vec!();
    for reference_wise_params in params {
        let mut bam_readers = vec![];
        let filter_params = FilterParameters::generate_from_clap(m);
        let index = coverm::bwa_index_maintenance::generate_bwa_index(
            reference_wise_params.reference);

        let reference = reference_wise_params.reference;
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.is_present("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.value_of("bam-file-cache-directory").unwrap(),
                        reference, naming_readset);
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        for p in reference_wise_params {
            bam_readers.push(
                coverm::bam_generator::generate_filtered_named_bam_readers_from_reads(
                    index.index_path(),
                    p.read1,
                    p.read2,
                    p.read_format.clone(),
                    p.threads,
                    bam_file_cache(p.read1).as_ref().map(String::as_ref),
                    filter_params.min_aligned_length,
                    filter_params.min_percent_identity,
                    filter_params.min_aligned_percent));
        }

        debug!("Finished BAM setup");
        let to_return = BamGeneratorSet {
            generators: bam_readers,
            index: index
        };
        generator_set.push(to_return);
    }
    return generator_set;
}



fn run_contig<R: coverm::bam_generator::NamedBamReader,
              T: coverm::bam_generator::NamedBamReaderGenerator<R>,
              C: coverm::coverage_takers::CoverageTaker>(
    coverage_estimators: &mut Vec<CoverageEstimator>,
    bam_readers: Vec<T>,
    print_zeros: bool,
    flag_filter: bool,
    coverage_taker: &mut C) {

    coverm::contig::contig_coverage(
        bam_readers,
        coverage_taker,
        coverage_estimators,
        print_zeros,
        flag_filter);
}


fn set_log_level(matches: &clap::ArgMatches, is_last: bool) {
    let mut log_level = LevelFilter::Info;
    let mut specified = false;
    if matches.is_present("verbose") {
        specified = true;
        log_level = LevelFilter::Debug;
    }
    if matches.is_present("quiet") {
        specified = true;
        log_level = LevelFilter::Error;
    }
    if specified || is_last {
        let mut builder = Builder::new();
        builder.filter_level(log_level);
        if env::var("RUST_LOG").is_ok() {
            builder.parse(&env::var("RUST_LOG").unwrap());
        }
        if builder.try_init().is_err() {
            panic!("Failed to set log level - has it been specified multiple times?")
        }
    }
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
   -x, --genome-fasta-extension <EXT>    File extension of genomes in the directory
                                         specified with -d/--genome-fasta-directory
                                         [default \"fna\"]
   --single-genome                       All contigs are from the same genome

Define mapping(s) (required):
  Either define BAM:
    -b, --bam-files <PATH> ..            Path to reference-sorted BAM file(s)

  Or do mapping:
   -r, --reference <PATH>                FASTA file of contigs or BWA index stem
   -t, --threads <INT>                   Number of threads to use for mapping
   -1 <PATH> ..                          Forward FASTA/Q file(s) for mapping
   -2 <PATH> ..                          Reverse FASTA/Q file(s) for mapping
   -c, --coupled <PATH> <PATH> ..        Forward and reverse pairs of FASTA/Q files(s)
                                         for mapping.
   --interleaved <PATH> ..               Interleaved FASTA/Q files(s) for mapping.
   --single <PATH> ..                    Unpaired FASTA/Q files(s) for mapping.

Alignment filtering (optional):
   --min-aligned-length <INT>            Exclude pairs with smaller numbers of
                                         aligned bases [default: 0]
   --min-percent-identity <FLOAT>        Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95% [default 0.0]
   --min-aligned-percent <FLOAT>         Exclude pairs by percent aligned
                                         identity e.g. 0.95 for 95% [default 0.0]

Other arguments (optional):
   -m, --methods <METHOD> [METHOD ..]    Method(s) for calculating coverage.
                                         One of:
                                              relative_abundance (default)
                                              mean
                                              trimmed_mean
                                              coverage_histogram
                                              covered_fraction
                                              variance
                                              length
   --min-covered-fraction FRACTION       Genomes with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0.10]
   --contig-end-exclusion                Exclude bases at the ends of reference
                                         sequences from calculation [default: 75]
   --trim-min FRACTION                   Remove this smallest fraction of positions
                                         when calculating trimmed_mean
                                         [default: 0.05]
   --trim-max FRACTION                   Maximum fraction for trimmed_mean
                                         calculations [default: 0.95]
   --no-zeros                            Omit printing of genomes that have zero
                                         coverage
   --no-flag-filter                      Do not ignore secondary and supplementary
                                         alignments, and improperly paired reads
   --bam-file-cache-directory            Output BAM files generated during
                                         alignment to this directory
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
   -r, --reference <PATH>                FASTA file of contigs or BWA index stem
   -t, --threads <INT>                   Number of threads to use for mapping
   -1 <PATH> ..                          Forward FASTA/Q file(s) for mapping
   -2 <PATH> ..                          Reverse FASTA/Q file(s) for mapping
   -c, --coupled <PATH> <PATH> ..        Forward and reverse pairs of FASTA/Q files(s)
                                         for mapping.
   --interleaved <PATH> ..               Interleaved FASTA/Q files(s) for mapping.
   --single <PATH> ..                    Unpaired FASTA/Q files(s) for mapping.

Alignment filtering (optional):
   --min-aligned-length <INT>            Exclude pairs with smaller numbers of
                                         aligned bases [default: none]
   --min-percent-identity <FLOAT>        Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95% [default: none]
   --min-aligned-percent <FLOAT>         Exclude pairs by percent aligned
                                         identity e.g. 0.95 for 95% [default 0.0]

Other arguments (optional):
   -m, --methods <METHOD> [METHOD ..]    Method(s) for calculating coverage.
                                         One of:
                                              mean (default)
                                              trimmed_mean
                                              coverage_histogram
                                              covered_fraction
                                              variance
                                              length
   --min-covered-fraction FRACTION       Genom
es with less coverage than this
                                         reported as having zero coverage.
                                         [default: 0]
   --contig-end-exclusion                Exclude bases at the ends of reference
                                         sequences from calculation [default: 75]
   --trim-min FRACTION                   Remove this smallest fraction of positions
                                         when calculating trimmed_mean
                                         [default: 0.05]
   --trim-max FRACTION                   Maximum fraction for trimmed_mean
                                         calculations [default: 0.95]
   --no-zeros                            Omit printing of genomes that have zero
                                         coverage
   --no-flag-filter                      Do not ignore secondary and supplementary
                                         alignments, and improperly paired reads
   --bam-file-cache-directory            Output BAM files generated during
                                         alignment to this directory
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
   --min-aligned-percent <FLOAT>       Exclude pairs by percent aligned
                                       identity e.g. 0.95 for 95% [default 0.0]

Other:
   -t, --threads <INT>                 Number of threads for output compression
                                       [default 1]
   -v, --verbose                       Print extra debugging information
   -q, --quiet                         Unless there is an error, do not print
                                       log messages

Ben J. Woodcroft <benjwoodcroft near gmail.com>";

    let make_help: &'static str =
        "coverm make: Generate BAM files through mapping.

Output (required):
   -o, --output-directory <DIR>          Where generated BAM files will go

Mapping parameters:
   -r, --reference <PATH>                FASTA file of contig(s) or BWA index stem
   -t, --threads <INT>                   Number of threads to use for mapping
   -1 <PATH> ..                          Forward FASTA/Q file(s) for mapping
   -2 <PATH> ..                          Reverse FASTA/Q file(s) for mapping
   -c, --coupled <PATH> <PATH> ..        Forward and reverse pairs of FASTA/Q files(s)
                                         for mapping.
   --interleaved <PATH> ..               Interleaved FASTA/Q files(s) for mapping.
   --single <PATH> ..                    Unpaired FASTA/Q files(s) for mapping.

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
\tmake\tGenerate BAM files through alignment
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
                     .requires("read2")
                     .required_unless_one(
                         &["bam-files","coupled","interleaved","single"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("read2")
                     .short("-2")
                     .multiple(true)
                     .takes_value(true)
                     .requires("read1")
                     .required_unless_one(
                         &["bam-files","coupled","interleaved","single"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("coupled")
                     .short("-c")
                     .long("coupled")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["bam-files","read1","interleaved","single"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("interleaved")
                     .long("interleaved")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["bam-files","read1","coupled","single"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("single")
                     .long("single")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["bam-files","read1","coupled","interleaved"])
                     .requires("no-flag-filter")
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("reference")
                     .short("-r")
                     .long("reference")
                     .takes_value(true)
                     .required_unless("bam-files")
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("bam-file-cache-directory")
                     .long("bam-file-cache-directory")
                     .takes_value(true)
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
                     .required_unless_one(
                         &["genome-fasta-files","genome-fasta-directory","single-genome"])
                     .takes_value(true))
                .arg(Arg::with_name("genome-fasta-files")
                     .short("f")
                     .long("genome-fasta-files")
                     .multiple(true)
                     .conflicts_with("separator")
                     .conflicts_with("genome-fasta-directory")
                     .conflicts_with("single-genome")
                     .required_unless_one(
                         &["separator","genome-fasta-directory","single-genome"])
                     .takes_value(true))
                .arg(Arg::with_name("genome-fasta-directory")
                     .short("d")
                     .long("genome-fasta-directory")
                     .conflicts_with("separator")
                     .conflicts_with("genome-fasta-files")
                     .conflicts_with("single-genome")
                     .required_unless_one(
                         &["genome-fasta-files","separator","single-genome"])
                     .takes_value(true))
                .arg(Arg::with_name("genome-fasta-extension")
                     .short("x")
                     .long("genome-fasta-extension")
                     // Unsure why, but uncommenting causes test failure - clap
                     // bug?
                     //.requires("genome-fasta-directory")
                     .default_value("fna")
                     .takes_value(true))
                .arg(Arg::with_name("single-genome")
                     .long("single-genome")
                     .conflicts_with("separator")
                     .conflicts_with("genome-fasta-files")
                     .conflicts_with("genome-fasta-directory"))

                .arg(Arg::with_name("min-aligned-length")
                     .long("min-aligned-length")
                     .takes_value(true)
                     .conflicts_with("no-flag-filter"))
                .arg(Arg::with_name("min-percent-identity")
                     .long("min-percent-identity")
                     .takes_value(true)
                     .conflicts_with("no-flag-filter"))
                .arg(Arg::with_name("min-aligned-percent")
                     .long("min-aligned-percent")
                     .takes_value(true)
                     .conflicts_with("no-flag-filter"))

                .arg(Arg::with_name("methods")
                     .short("m")
                     .long("method")
                     .long("methods")
                     .takes_value(true)
                     .multiple(true)
                     .possible_values(&[
                         "relative_abundance",
                         "mean",
                         "trimmed_mean",
                         "coverage_histogram",
                         "covered_fraction",
                         "variance",
                         "length"])
                     .default_value("relative_abundance"))
                .arg(Arg::with_name("trim-min")
                     .long("trim-min")
                     .default_value("0.05"))
                .arg(Arg::with_name("trim-max")
                     .long("trim-max")
                     .default_value("0.95"))
                .arg(Arg::with_name("min-covered-fraction")
                     .long("min-covered-fraction")
                     .default_value("0.10"))
                .arg(Arg::with_name("contig-end-exclusion")
                     .long("contig-end-exclusion")
                     .default_value("75"))
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
                     .takes_value(true))
                .arg(Arg::with_name("read1")
                     .short("-1")
                     .multiple(true)
                     .takes_value(true)
                     .requires("read2")
                     .required_unless_one(
                         &["bam-files","coupled","interleaved","single"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("read2")
                     .short("-2")
                     .multiple(true)
                     .takes_value(true)
                     .requires("read1")
                     .required_unless_one(
                         &["bam-files","coupled","interleaved","single"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("coupled")
                     .short("-c")
                     .long("coupled")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["bam-files","read1","interleaved","single"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("interleaved")
                     .long("interleaved")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["bam-files","read1","coupled","single"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("single")
                     .long("single")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["bam-files","read1","coupled","interleaved"])
                     .requires("no-flag-filter")
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("reference")
                     .short("-r")
                     .long("reference")
                     .takes_value(true)
                     .required_unless("bam-files")
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("bam-file-cache-directory")
                     .long("bam-file-cache-directory")
                     .takes_value(true)
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
                .arg(Arg::with_name("min-aligned-percent")
                     .long("min-aligned-percent")
                     .takes_value(true)
                     .conflicts_with("no-flag-filter"))

                .arg(Arg::with_name("methods")
                     .short("m")
                     .long("method")
                     .long("methods")
                     .takes_value(true)
                     .multiple(true)
                     .possible_values(&[
                         "mean",
                         "trimmed_mean",
                         "coverage_histogram",
                         "covered_fraction",
                         "variance",
                         "length"])
                     .default_value("mean"))
                .arg(Arg::with_name("min-covered-fraction")
                     .long("min-covered-fraction")
                     .default_value("0.0"))
                .arg(Arg::with_name("contig-end-exclusion")
                     .long("contig-end-exclusion")
                     .default_value("75"))
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
                .arg(Arg::with_name("min-aligned-percent")
                     .long("min-aligned-percent")
                     .takes_value(true)
                     .conflicts_with("no-flag-filter")
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
        .subcommand(
            SubCommand::with_name("make")
                .about("Generate BAM files through mapping")
                .help(make_help)
                .arg(Arg::with_name("output-directory")
                     .short("-o")
                     .long("output-directory")
                     .takes_value(true)
                     .required(true))
                .arg(Arg::with_name("read1")
                     .short("-1")
                     .multiple(true)
                     .takes_value(true)
                     .requires("read2")
                     .required_unless_one(
                         &["coupled","interleaved","single"]))
                .arg(Arg::with_name("read2")
                     .short("-2")
                     .multiple(true)
                     .takes_value(true)
                     .requires("read1")
                     .required_unless_one(
                         &["coupled","interleaved","single"]))
                .arg(Arg::with_name("coupled")
                     .short("-c")
                     .long("coupled")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["read1","interleaved","single"]))
                .arg(Arg::with_name("interleaved")
                     .long("interleaved")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["read1","coupled","single"]))
                .arg(Arg::with_name("single")
                     .long("single")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["read1","coupled","interleaved"]))
                .arg(Arg::with_name("reference")
                     .short("-r")
                     .long("reference")
                     .takes_value(true)
                     .required(true))
                .arg(Arg::with_name("threads")
                     .short("-t")
                     .long("threads")
                     .default_value("1")
                     .takes_value(true))
        )
}
