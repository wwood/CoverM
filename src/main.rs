extern crate coverm;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::bam_generator::*;
use coverm::filter;
use coverm::external_command_checker;
use coverm::coverage_takers::*;
use coverm::mapping_parameters::*;
use coverm::coverage_printer::*;
use coverm::FlagFilter;

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
extern crate env_logger;
use log::LevelFilter;
use env_logger::Builder;

extern crate tempfile;


fn main(){
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    let mut print_stream = &mut std::io::stdout();
    set_log_level(&matches, false);

    match matches.subcommand_name() {

        Some("genome") => {
            let m = matches.subcommand_matches("genome").unwrap();
            set_log_level(m, true);

            let mut estimators_and_taker = EstimatorsAndTaker::generate_from_clap(
                m, print_stream);
            estimators_and_taker = estimators_and_taker.print_headers(
                &"Genome", &mut std::io::stdout());
            let filter_params = FilterParameters::generate_from_clap(m);

            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
                if filter_params.doing_filtering() {
                    run_genome(
                        coverm::bam_generator::generate_filtered_bam_readers_from_bam_files(
                            bam_files,
                            filter_params.flag_filters,
                            filter_params.min_aligned_length_single,
                            filter_params.min_percent_identity_single,
                            filter_params.min_aligned_percent_single,
                            filter_params.min_aligned_length_pair,
                            filter_params.min_percent_identity_pair,
                            filter_params.min_aligned_percent_pair),
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
                if filter_params.doing_filtering() {
                    debug!("Mapping and filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(m);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
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
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
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
            let filter_params = FilterParameters::generate_from_clap(m);

            let mut estimators_and_taker = EstimatorsAndTaker::generate_from_clap(
                m, &mut print_stream);
            estimators_and_taker = estimators_and_taker.print_headers(
                &"Contig", &mut std::io::stdout());

            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
                if filter_params.doing_filtering() {
                    let mut bam_readers = coverm::bam_generator::generate_filtered_bam_readers_from_bam_files(
                        bam_files,
                        filter_params.flag_filters.clone(),
                        filter_params.min_aligned_length_single,
                        filter_params.min_percent_identity_single,
                        filter_params.min_aligned_percent_single,
                        filter_params.min_aligned_length_pair,
                        filter_params.min_percent_identity_pair,
                        filter_params.min_aligned_percent_pair);
                    run_contig(
                        &mut estimators_and_taker,
                        bam_readers,
                        print_zeros,
                        filter_params.flag_filters);
                } else {
                    let mut bam_readers = coverm::bam_generator::generate_named_bam_readers_from_bam_files(
                        bam_files);
                    run_contig(
                        &mut estimators_and_taker,
                        bam_readers,
                        print_zeros,
                        filter_params.flag_filters);
                }
            } else {
                external_command_checker::check_for_bwa();
                external_command_checker::check_for_samtools();
                if filter_params.doing_filtering() {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(m);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    debug!("Finished collecting generators.");
                    run_contig(
                        &mut estimators_and_taker,
                        all_generators,
                        print_zeros,
                        filter_params.flag_filters);
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m);
                    let mut all_generators = vec!();
                    let mut indices = vec!(); // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    run_contig(
                        &mut estimators_and_taker,
                        all_generators,
                        print_zeros,
                        filter_params.flag_filters.clone());
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

            let filter_params = FilterParameters::generate_from_clap(m);

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
                    reader,
                    filter_params.flag_filters.clone(),
                    filter_params.min_aligned_length_single,
                    filter_params.min_percent_identity_single,
                    filter_params.min_aligned_percent_single,
                    filter_params.min_aligned_length_pair,
                    filter_params.min_percent_identity_pair,
                    filter_params.min_aligned_percent_pair);

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
            setup_bam_cache_directory(output_directory);
            let params = MappingParameters::generate_from_clap(&m);
            let mut generator_sets = vec!();
            let discard_unmapped_reads = m.is_present("discard-unmapped");

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
                                output_directory, p.reference, p.read1),
                            discard_unmapped_reads,
                            p.bwa_options));
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
    printer: CoveragePrinter,
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

        let taker;
        let output_format = m.value_of("output-format").unwrap();
        let printer;

        if doing_metabat(&m) {
            estimators.push(CoverageEstimator::new_estimator_length());
            estimators.push(CoverageEstimator::new_estimator_mean(
                min_fraction_covered, contig_end_exclusion, false));
            estimators.push(CoverageEstimator::new_estimator_variance(
                min_fraction_covered, contig_end_exclusion));

            debug!("Cached regular coverage taker for metabat mode being used");
            taker = CoverageTakerType::new_cached_single_float_coverage_taker(
                estimators.len());
            printer = CoveragePrinter::MetabatAdjustedCoveragePrinter;
        } else {
            for (i, method) in methods.iter().enumerate() {
                match method {
                    &"mean" => {
                        estimators.push(CoverageEstimator::new_estimator_mean(
                            min_fraction_covered,
                            contig_end_exclusion,
                            false)); // TODO: Parameterise exclude_mismatches
                    },
                    &"coverage_histogram" => {
                        estimators.push(CoverageEstimator::new_estimator_pileup_counts(
                            min_fraction_covered, contig_end_exclusion));
                    },
                    &"trimmed_mean" => {
                        let min = value_t!(m.value_of("trim-min"), f32).unwrap();
                        let max = value_t!(m.value_of("trim-max"), f32).unwrap();
                        if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                            panic!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
                        }
                        estimators.push(CoverageEstimator::new_estimator_trimmed_mean(
                            min, max, min_fraction_covered, contig_end_exclusion));
                    },
                    &"covered_fraction" => {
                        estimators.push(CoverageEstimator::new_estimator_covered_fraction(
                            min_fraction_covered, contig_end_exclusion));
                    },
                    &"variance" =>{
                        estimators.push(CoverageEstimator::new_estimator_variance(
                            min_fraction_covered, contig_end_exclusion));
                    },
                    &"length" =>{
                        estimators.push(CoverageEstimator::new_estimator_length());
                    },
                    &"relative_abundance" => {
                        columns_to_normalise.push(i);
                        estimators.push(CoverageEstimator::new_estimator_mean(
                            min_fraction_covered, contig_end_exclusion, false));
                        // TODO: Parameterise exclude_mismatches
                    },
                    &"count" => {
                        estimators.push(CoverageEstimator::new_estimator_read_count());
                    },
                    _ => panic!("programming error")
                };
            }

            if methods.contains(&"coverage_histogram") {
                if methods.len() > 1 {
                    panic!("Cannot specify the coverage_histogram method with any other coverage methods")
                } else {
                    debug!("Coverage histogram type coverage taker being used");
                    taker = CoverageTakerType::new_pileup_coverage_coverage_printer(stream);
                    printer = CoveragePrinter::StreamedCoveragePrinter;
                }
            } else if columns_to_normalise.len() == 0 && output_format == "sparse" {
                debug!("Streaming regular coverage output");
                taker = CoverageTakerType::new_single_float_coverage_streaming_coverage_printer(
                    stream);
                printer = CoveragePrinter::StreamedCoveragePrinter;
            } else {
                debug!("Cached regular coverage taker with columns to normlise: {:?}",
                       columns_to_normalise);
                taker = CoverageTakerType::new_cached_single_float_coverage_taker(
                    estimators.len());
                printer = match output_format {
                    "sparse" => CoveragePrinter::SparseCachedCoveragePrinter,
                    "dense" => CoveragePrinter::DenseCachedCoveragePrinter {
                        entry_type: None, estimator_headers: None},
                    _ => panic!("Unexpected output format seen. Programming error")
                }
            }

        }

        return EstimatorsAndTaker {
            estimators: estimators,
            taker: taker,
            columns_to_normalise: columns_to_normalise,
            printer: printer,
        }
    }

    pub fn print_headers(
        mut self,
        entry_type: &str,
        print_stream: &mut std::io::Write) -> Self {

        let mut headers: Vec<String> = vec!();
        for e in self.estimators.iter() {
            for h in e.column_headers() {
                headers.push(h.to_string())
            }
        }
        for i in self.columns_to_normalise.iter() {
            headers[*i] = "Relative Abundance (%)".to_string();
        }
        self.printer.print_headers(&entry_type, headers, print_stream);
        return self;
    }
}


fn run_genome<'a,
              R: coverm::bam_generator::NamedBamReader,
              T: coverm::bam_generator::NamedBamReaderGenerator<R>>(
    bam_generators: Vec<T>,
    m: &clap::ArgMatches,
    estimators_and_taker: &'a mut EstimatorsAndTaker<'a>) {

    let print_zeros = !m.is_present("no-zeros");
    let proper_pairs_only = m.is_present("proper-pairs-only");
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
                proper_pairs_only,
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
                                    "Not using directory entry '{}' as a genome FASTA file, as \
                                     it does not end with the extension '{}'",
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
                proper_pairs_only,
                &mut estimators_and_taker.estimators)
        }
    };

    debug!("Finalising printing ..");
    estimators_and_taker.printer.finalise_printing(
        &estimators_and_taker.taker, &mut std::io::stdout(), Some(&reads_mapped),
        &estimators_and_taker.columns_to_normalise);
}

fn doing_metabat(m: &clap::ArgMatches) -> bool {
    match m.subcommand_name() {
        Some("contig") | None => {
            if !m.is_present("methods") { return false; }
            let methods: Vec<&str> = m.values_of("methods").unwrap().collect();
            if methods.contains(&"metabat") {
                if methods.len() > 1 {
                    panic!("Cannot specify the metabat method with any other coverage methods")
                } else {
                    return true
                }
            }
            return false
        },
        _ => {
            debug!("Not running in contig mode so cannot be in metabat mode");
            return false
        }
    }
}

#[derive(Debug)]
struct FilterParameters {
    flag_filters: FlagFilter,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
}
impl FilterParameters {
    pub fn generate_from_clap(m: &clap::ArgMatches) -> FilterParameters {
        let mut f = FilterParameters {
            flag_filters: FlagFilter {
                include_improper_pairs: !m.is_present("proper-pairs-only"),
                include_secondary: false,
                include_supplementary: false,
            },
            min_aligned_length_single: match m.is_present("min-aligned-length") {
                true => value_t!(m.value_of("min-aligned-length"), u32).unwrap(),
                false => 0
            },
            min_percent_identity_single: match m.is_present("min-percent-identity") {
                true => value_t!(m.value_of("min-percent-identity"), f32).unwrap(),
                false => 0.0
            },
            min_aligned_percent_single: match m.is_present("min-aligned-percent") {
                true => value_t!(m.value_of("min-aligned-percent"), f32).unwrap(),
                false => 0.0
            },
            min_aligned_length_pair: match m.is_present("min-aligned-length-pair") {
                true => value_t!(m.value_of("min-aligned-length-pair"), u32).unwrap(),
                false => 0
            },
            min_percent_identity_pair: match m.is_present("min-percent-identity-pair") {
                true => value_t!(m.value_of("min-percent-identity-pair"), f32).unwrap(),
                false => 0.0
            },
            min_aligned_percent_pair: match m.is_present("min-aligned-percent-pair") {
                true => value_t!(m.value_of("min-aligned-percent-pair"), f32).unwrap(),
                false => 0.0
            }
        };
        if doing_metabat(&m) {
            debug!("Setting single read percent identity threshold at 0.97 for \
                    MetaBAT adjusted coverage.");
            // we use >= where metabat uses >. Gah.
            f.min_percent_identity_single = 0.97001;
            f.flag_filters.include_improper_pairs = true;
            f.flag_filters.include_supplementary = true;
            f.flag_filters.include_secondary = true;
        }
        debug!("Filter parameters set as {:?}", f);
        return f
    }

    pub fn doing_filtering(&self) -> bool {
        return self.min_percent_identity_single > 0.0 || self.min_percent_identity_pair > 0.0 ||
            self.min_aligned_percent_single > 0.0 || self.min_aligned_percent_pair > 0.0 ||
            self.min_aligned_length_single > 0 || self.min_aligned_length_pair > 0
    }
}



fn get_streamed_bam_readers<'a>(
    m: &clap::ArgMatches) -> Vec<BamGeneratorSet<StreamingNamedBamReaderGenerator>> {

    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");

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
                    bam_file_cache(p.read1).as_ref().map(String::as_ref),
                    discard_unmapped,
                    p.bwa_options));
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
    let discard_unmapped = m.is_present("discard-unmapped");

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
                    filter_params.flag_filters.clone(),
                    filter_params.min_aligned_length_single,
                    filter_params.min_percent_identity_single,
                    filter_params.min_aligned_percent_single,
                    filter_params.min_aligned_length_pair,
                    filter_params.min_percent_identity_pair,
                    filter_params.min_aligned_percent_pair,
                    p.bwa_options,
                    discard_unmapped));
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



fn run_contig<'a,
              R: coverm::bam_generator::NamedBamReader,
              T: coverm::bam_generator::NamedBamReaderGenerator<R>>(
    estimators_and_taker: &'a mut EstimatorsAndTaker<'a>,
    bam_readers: Vec<T>,
    print_zeros: bool,
    flag_filters: FlagFilter) {

    coverm::contig::contig_coverage(
        bam_readers,
        &mut estimators_and_taker.taker,
        &mut estimators_and_taker.estimators,
        print_zeros,
        flag_filters);

    debug!("Finalising printing ..");

    estimators_and_taker.printer.finalise_printing(
        &estimators_and_taker.taker,
        &mut std::io::stdout(),
        None,
        &estimators_and_taker.columns_to_normalise);
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
   -f, --genome-fasta-files <PATH> ..    Path to FASTA files of each genome e.g.
                                         'pathA/genome1.fna pathB/genome2.fa'
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
                                         e.g. assembly output
   -t, --threads <INT>                   Number of threads to use for mapping
   -1 <PATH> ..                          Forward FASTA/Q file(s) for mapping
   -2 <PATH> ..                          Reverse FASTA/Q file(s) for mapping
   -c, --coupled <PATH> <PATH> ..        One or more pairs of forward and reverse
                                         FASTA/Q files for mapping in order
                                         <sample1_R1.fq.gz> <sample1_R2.fq.gz>
                                         <sample2_R1.fq.gz> <sample2_R2.fq.gz> ..
   --interleaved <PATH> ..               Interleaved FASTA/Q files(s) for mapping.
   --single <PATH> ..                    Unpaired FASTA/Q files(s) for mapping.
   --bwa-params PARAMS                   Extra parameters to provide to BWA. Note
                                         that usage of this parameter has security
                                         implications if untrusted input is specified.
                                         [default \"\"]

Alignment filtering (optional):
   --min-aligned-length <INT>            Exclude reads with smaller numbers of
                                         aligned bases [default: 0]
   --min-percent-identity <FLOAT>        Exclude reads by overall percent
                                         identity e.g. 0.95 for 95% [default 0.0]
   --min-aligned-percent <FLOAT>         Exclude reads by percent aligned
                                         identity e.g. 0.95 for 95% [default 0.0]
   --min-aligned-length-pair <INT>       Exclude pairs with smaller numbers of
                                         aligned bases [default: 0]
   --min-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95% [default 0.0]
   --min-aligned-percent-pair <FLOAT>    Exclude pairs by percent aligned
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
                                              count
   --output-format FORMAT                Shape of output: 'sparse' for long format,
                                         'dense' for species-by-site.
                                         [default: sparse]
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
   --proper-pairs-only                   Require reads to be mapped as proper pairs
   --bam-file-cache-directory            Output BAM files generated during
                                         alignment to this directory
   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Example usage:

  coverm genome -f bins/*fa -c read1.fq.gz read2.fq.gz -r concatenated_genomes.fna

Ben J. Woodcroft <benjwoodcroft near gmail.com>
";

    let contig_help: &'static str =
        "coverm contig: Calculate read coverage per-contig

Define mapping(s) (required):
  Either define BAM:
    -b, --bam-files <PATH> ..            Path to reference-sorted BAM file(s)

  Or do mapping:
   -r, --reference <PATH>                FASTA file of contigs or BWA index stem
                                         e.g. concatenated genomes or assembly
   -t, --threads <INT>                   Number of threads to use for mapping
   -1 <PATH> ..                          Forward FASTA/Q file(s) for mapping
   -2 <PATH> ..                          Reverse FASTA/Q file(s) for mapping
   -c, --coupled <PATH> <PATH> ..        One or more pairs of forward and reverse
                                         FASTA/Q files for mapping in order
                                         <sample1_R1.fq.gz> <sample1_R2.fq.gz>
                                         <sample2_R1.fq.gz> <sample2_R2.fq.gz> ..
   --interleaved <PATH> ..               Interleaved FASTA/Q files(s) for mapping.
   --single <PATH> ..                    Unpaired FASTA/Q files(s) for mapping.
   --bwa-params PARAMS                   Extra parameters to provide to BWA. Note
                                         that usage of this parameter has security
                                         implications if untrusted input is specified.
                                         [default \"\"]

Alignment filtering (optional):
   --min-aligned-length <INT>            Exclude reads with smaller numbers of
                                         aligned bases [default: 0]
   --min-percent-identity <FLOAT>        Exclude reads by overall percent
                                         identity e.g. 0.95 for 95% [default 0.0]
   --min-aligned-percent <FLOAT>         Exclude reads by percent aligned
                                         identity e.g. 0.95 for 95% [default 0.0]
   --min-aligned-length-pair <INT>       Exclude pairs with smaller numbers of
                                         aligned bases.
                                         Implies --proper-pairs-only.[default: 0]
   --min-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95%.
                                         Implies --proper-pairs-only. [default 0.0]
   --min-aligned-percent-pair <FLOAT>    Exclude pairs by percent aligned
                                         identity e.g. 0.95 for 95%.
                                         Implies --proper-pairs-only. [default 0.0]

Other arguments (optional):
   -m, --methods <METHOD> [METHOD ..]    Method(s) for calculating coverage.
                                         One of:
                                           mean (default)
                                           trimmed_mean
                                           coverage_histogram
                                           covered_fraction
                                           variance
                                           length
                                           count
                                           metabat (\"MetaBAT adjusted coverage\")
   --output-format FORMAT                Shape of output: 'sparse' for long format,
                                         'dense' for species-by-site.
                                         [default: sparse]
   --min-covered-fraction FRACTION       Contigs with less coverage than this
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
   --proper-pairs-only                   Require reads to be mapped as proper pairs
   --bam-file-cache-directory            Output BAM files generated during
                                         alignment to this directory
   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Example usage:

  coverm contig -b mapping.bam

Ben J. Woodcroft <benjwoodcroft near gmail.com>";

    let filter_help: &'static str =
        "coverm filter: Remove alignments with insufficient identity.

Only primary, non-supplementary alignments are considered, and output files
are grouped by reference, but not sorted by position.

Files (both required):
   -b, --bam-files <PATH> ..             Path to reference-sorted BAM file(s)
   -o, --output-bam-files <PATH> ..      Path to corresponding output file(s)

Thresholds:
   --min-aligned-length <INT>            Exclude reads with smaller numbers of
                                         aligned bases [default: 0]
   --min-percent-identity <FLOAT>        Exclude reads by overall percent
                                         identity e.g. 0.95 for 95% [default 0.0]
   --min-aligned-percent <FLOAT>         Exclude reads by percent aligned
                                         identity e.g. 0.95 for 95% [default 0.0]
   --min-aligned-length-pair <INT>       Exclude pairs with smaller numbers of
                                         aligned bases.
                                         Implies --proper-pairs-only.[default: 0]
   --min-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95%.
                                         Implies --proper-pairs-only. [default 0.0]
   --min-aligned-percent-pair <FLOAT>    Exclude pairs by percent aligned
                                         identity e.g. 0.95 for 95%.
                                         Implies --proper-pairs-only. [default 0.0]
   --proper-pairs-only                   Require reads to be mapped as proper pairs

Other:
   -t, --threads <INT>                   Number of threads for output compression
                                         [default 1]
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Example usage:

  coverm filter -b in.bam -o out.bam --min-aligned-length 75

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
   -c, --coupled <PATH> <PATH> ..        One or more pairs of forward and reverse
                                         FASTA/Q files for mapping in order
                                         <sample1_R1.fq.gz> <sample1_R2.fq.gz>
                                         <sample2_R1.fq.gz> <sample2_R2.fq.gz> ..
   --interleaved <PATH> ..               Interleaved FASTA/Q files(s) for mapping.
   --single <PATH> ..                    Unpaired FASTA/Q files(s) for mapping.
   --bwa-params PARAMS                   Extra parameters to provide to BWA. Note
                                         that usage of this parameter has security
                                         implications if untrusted input is specified.
                                         [default \"\"]
   --discard-unmapped                    Exclude unmapped reads from generated BAM files.

Example usage:

  coverm make -r combined_genomes.fna -1 read1.fq -2 read2.fq

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
                .arg(Arg::with_name("bwa-params")
                     .long("bwa-params")
                     .long("bwa-parameters")
                     .takes_value(true)
                     .allow_hyphen_values(true)
                     .requires("reference"))
                .arg(Arg::with_name("discard-unmapped")
                     .long("discard-unmapped")
                     .requires("bam-file-cache-directory"))

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
                     .takes_value(true))
                .arg(Arg::with_name("min-percent-identity")
                     .long("min-percent-identity")
                     .takes_value(true))
                .arg(Arg::with_name("min-aligned-percent")
                     .long("min-aligned-percent")
                     .takes_value(true))
                .arg(Arg::with_name("min-aligned-length-pair")
                     .long("min-aligned-length-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))
                .arg(Arg::with_name("min-percent-identity-pair")
                     .long("min-percent-identity-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))
                .arg(Arg::with_name("min-aligned-percent-pair")
                     .long("min-aligned-percent-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))

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
                         "length",
                         "count"])
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
                .arg(Arg::with_name("proper-pairs-only")
                     .long("proper-pairs-only"))
                .arg(Arg::with_name("output-format")
                     .long("output-format")
                     .possible_values(&["sparse","dense"])
                     .default_value("sparse"))

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
                .arg(Arg::with_name("bwa-params")
                     .long("bwa-params")
                     .long("bwa-parameters")
                     .takes_value(true)
                     .allow_hyphen_values(true)
                     .requires("reference"))
                .arg(Arg::with_name("discard-unmapped")
                     .long("discard-unmapped")
                     .requires("bam-file-cache-directory"))

                .arg(Arg::with_name("min-aligned-length")
                     .long("min-aligned-length")
                     .takes_value(true))
                .arg(Arg::with_name("min-percent-identity")
                     .long("min-percent-identity")
                     .takes_value(true))
                .arg(Arg::with_name("min-aligned-percent")
                     .long("min-aligned-percent")
                     .takes_value(true))
                .arg(Arg::with_name("min-aligned-length-pair")
                     .long("min-aligned-length-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))
                .arg(Arg::with_name("min-percent-identity-pair")
                     .long("min-percent-identity-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))
                .arg(Arg::with_name("min-aligned-percent-pair")
                     .long("min-aligned-percent-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))

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
                         "length",
                         "count",
                         "metabat"])
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
                .arg(Arg::with_name("proper-pairs-only")
                     .long("proper-pairs-only"))
                .arg(Arg::with_name("output-format")
                     .long("output-format")
                     .possible_values(&["sparse","dense"])
                     .default_value("sparse"))
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
                     .takes_value(true))
                .arg(Arg::with_name("min-percent-identity")
                     .long("min-percent-identity")
                     .takes_value(true))
                .arg(Arg::with_name("min-aligned-percent")
                     .long("min-aligned-percent")
                     .takes_value(true))
                .arg(Arg::with_name("min-aligned-length-pair")
                     .long("min-aligned-length-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))
                .arg(Arg::with_name("min-percent-identity-pair")
                     .long("min-percent-identity-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))
                .arg(Arg::with_name("min-aligned-percent-pair")
                     .long("min-aligned-percent-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))

                .arg(Arg::with_name("proper-pairs-only")
                     .long("proper-pairs-only"))
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
                .arg(Arg::with_name("discard-unmapped")
                     .long("discard-unmapped"))
                .arg(Arg::with_name("bwa-params")
                     .long("bwa-params")
                     .long("bwa-parameters")
                     .takes_value(true)
                     .allow_hyphen_values(true)
                     .requires("reference"))
        )
}
