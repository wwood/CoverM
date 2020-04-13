extern crate coverm;
use coverm::bam_generator::*;
use coverm::cli::*;
use coverm::coverage_printer::*;
use coverm::coverage_takers::*;
use coverm::external_command_checker;
use coverm::filter;
use coverm::genome_exclusion::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::mapping_parameters::*;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::shard_bam_reader::*;
use coverm::FlagFilter;
use coverm::CONCATENATED_FASTA_FILE_SEPARATOR;

extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;

use std::collections::HashSet;
use std::env;
use std::io::Write;
use std::process;
use std::str;

extern crate clap;
use clap::*;

#[macro_use]
extern crate log;
extern crate env_logger;
use env_logger::Builder;
use log::LevelFilter;

extern crate tempfile;
use tempfile::NamedTempFile;

const CONCATENATED_REFERENCE_CACHE_STEM: &str = "coverm-genome";
const DEFAULT_MAPPING_SOFTWARE_ENUM: MappingProgram = MappingProgram::MINIMAP2_SR;

fn galah_command_line_definition(
) -> galah::cluster_argument_parsing::GalahClustererCommandDefinition {
    galah::cluster_argument_parsing::GalahClustererCommandDefinition {
        dereplication_ani_argument: "dereplication-ani".to_string(),
        dereplication_prethreshold_ani_argument: "dereplication-prethreshold-ani".to_string(),
        dereplication_quality_formula_argument: "dereplication-quality-formula".to_string(),
        dereplication_precluster_method_argument: "dereplication-precluster-method".to_string(),
    }
}

fn main() {
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    let mut print_stream = &mut std::io::stdout();
    set_log_level(&matches, false);

    match matches.subcommand_name() {
        Some("genome") => {
            let m = matches.subcommand_matches("genome").unwrap();
            if m.is_present("full-help") {
                println!("{}", genome_full_help());
                process::exit(1);
            }
            set_log_level(m, true);

            let genome_names_content: Vec<u8>;

            let mut estimators_and_taker = EstimatorsAndTaker::generate_from_clap(m, print_stream);
            estimators_and_taker =
                estimators_and_taker.print_headers(&"Genome", &mut std::io::stdout());
            let filter_params = FilterParameters::generate_from_clap(m);
            let separator = parse_separator(m);

            let genomes_and_contigs_option_predereplication = if !m.is_present("separator")
                && !m.is_present("dereplicate")
                && !m.is_present("single-genome")
            {
                parse_all_genome_definitions(&m)
            } else {
                None
            };

            // This would be better as a separate function to make this function
            // smaller, but I find this hard because functions cannot return a
            // trait.
            let mut genome_exclusion_filter_separator_type: Option<SeparatorGenomeExclusionFilter> =
                None;
            let mut genome_exclusion_filter_non_type: Option<NoExclusionGenomeFilter> = None;
            let mut genome_exclusion_genomes_and_contigs: Option<GenomesAndContigsExclusionFilter> =
                None;
            enum GenomeExclusionTypes {
                SeparatorType,
                NoneType,
                GenomesAndContigsType,
            }
            let genome_exclusion_type = {
                if m.is_present("sharded") {
                    if m.is_present("exclude-genomes-from-deshard") {
                        let filename = m.value_of("exclude-genomes-from-deshard").unwrap();
                        genome_names_content = std::fs::read(filename).expect(&format!(
                            "Failed to open file '{}' containing list of excluded genomes",
                            filename
                        ));
                        let mut genome_names_hash: HashSet<&[u8]> = HashSet::new();
                        for n in genome_names_content.split(|s| *s == b"\n"[0]) {
                            if n != b"" {
                                genome_names_hash.insert(n);
                            }
                        }
                        if genome_names_hash.is_empty() {
                            warn!("No genomes read in that are to be excluded from desharding process");
                            genome_exclusion_filter_non_type = Some(NoExclusionGenomeFilter {});
                            GenomeExclusionTypes::NoneType
                        } else {
                            info!("Read in {} distinct genomes to exclude from desharding process e.g. '{}'",
                                  genome_names_hash.len(),
                                  std::str::from_utf8(genome_names_hash.iter().next().unwrap())
                                  .unwrap());
                            if separator.is_some() {
                                genome_exclusion_filter_separator_type =
                                    Some(SeparatorGenomeExclusionFilter {
                                        split_char: separator.unwrap(),
                                        excluded_genomes: genome_names_hash,
                                    });
                                GenomeExclusionTypes::SeparatorType
                            } else {
                                match genomes_and_contigs_option_predereplication {
                                    Some(ref gc) => {
                                        genome_exclusion_genomes_and_contigs =
                                            Some(GenomesAndContigsExclusionFilter {
                                                genomes_and_contigs: gc,
                                                excluded_genomes: genome_names_hash,
                                            });
                                        GenomeExclusionTypes::GenomesAndContigsType
                                    }
                                    None => unreachable!(),
                                }
                            }
                        }
                    } else {
                        debug!("Not excluding any genomes during the deshard process");
                        genome_exclusion_filter_non_type = Some(NoExclusionGenomeFilter {});
                        GenomeExclusionTypes::NoneType
                    }
                } else {
                    genome_exclusion_filter_non_type = Some(NoExclusionGenomeFilter {});
                    GenomeExclusionTypes::NoneType
                }
            };

            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();

                // Associate genomes and contig names, if required
                let genomes_and_contigs_option = parse_all_genome_definitions(&m);

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
                            filter_params.min_aligned_percent_pair,
                        ),
                        m,
                        &mut estimators_and_taker,
                        separator,
                        &genomes_and_contigs_option,
                    );
                } else if m.is_present("sharded") {
                    external_command_checker::check_for_samtools();
                    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
                    // Seems crazy, but I cannot work out how to make this more
                    // DRY, without making GenomeExclusion into an enum.
                    match genome_exclusion_type {
                        GenomeExclusionTypes::NoneType => {
                            run_genome(
                                coverm::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                                    bam_files,
                                    sort_threads,
                                    &genome_exclusion_filter_non_type.unwrap()),
                                m,
                                &mut estimators_and_taker,
                                separator,
                                &genomes_and_contigs_option);
                        }
                        GenomeExclusionTypes::SeparatorType => {
                            run_genome(
                                coverm::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                                    bam_files,
                                    sort_threads,
                                    &genome_exclusion_filter_separator_type.unwrap()),
                                m,
                                &mut estimators_and_taker,
                                separator,
                                &genomes_and_contigs_option);
                        }
                        GenomeExclusionTypes::GenomesAndContigsType => {
                            run_genome(
                                coverm::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                                    bam_files,
                                    sort_threads,
                                    &genome_exclusion_genomes_and_contigs.unwrap()),
                                m,
                                &mut estimators_and_taker,
                                separator,
                                &genomes_and_contigs_option);
                        }
                    }
                } else {
                    run_genome(
                        coverm::bam_generator::generate_named_bam_readers_from_bam_files(bam_files),
                        m,
                        &mut estimators_and_taker,
                        separator,
                        &genomes_and_contigs_option,
                    );
                }
            } else {
                let mapping_program = parse_mapping_program(&m);
                external_command_checker::check_for_samtools();

                // If genomes defined by file, then potentially dereplicate. If
                // no reference is specified, make a new concatenated one.

                // Get the GenomesAndContigs from the dereplicated genomes or
                // list of genome fasta files.

                let genome_fasta_files_opt = {
                    match bird_tool_utils::clap_utils::parse_list_of_genome_fasta_files(&m, false) {
                        Ok(paths) => {
                            if paths.len() == 0 {
                                error!(
                                    "Genome paths were described, but ultimately none were found"
                                );
                                process::exit(1);
                            }
                            if m.is_present("checkm-tab-table") || m.is_present("genome-info") {
                                let genomes_after_filtering =
                                    galah::cluster_argument_parsing::filter_genomes_through_checkm(
                                        &paths,
                                        &m,
                                        &galah_command_line_definition(),
                                    )
                                    .expect("Error parsing CheckM-related options");
                                info!(
                                    "After filtering by CheckM, {} genomes remained",
                                    genomes_after_filtering.len()
                                );
                                if genomes_after_filtering.len() == 0 {
                                    error!("All genomes were filtered out, so none remain to be mapped to");
                                    process::exit(1);
                                }
                                Some(
                                    genomes_after_filtering
                                        .iter()
                                        .map(|s| s.to_string())
                                        .collect(),
                                )
                            } else {
                                Some(paths)
                            }
                        }
                        Err(_) => None,
                    }
                };

                let (concatenated_genomes, genomes_and_contigs_option) = match m
                    .is_present("reference")
                {
                    true => match genome_fasta_files_opt {
                        Some(genome_paths) => (
                            None,
                            extract_genomes_and_contigs_option(
                                &m,
                                &genome_paths.iter().map(|s| s.as_str()).collect(),
                            ),
                        ),
                        None => (None, None),
                    },
                    false => {
                        // Dereplicate if required
                        let dereplicated_genomes: Vec<String> = if m.is_present("dereplicate") {
                            dereplicate(&m, &genome_fasta_files_opt.unwrap())
                        } else {
                            genome_fasta_files_opt.unwrap()
                        };
                        info!("Profiling {} genomes", dereplicated_genomes.len());

                        let list_of_genome_fasta_files = &dereplicated_genomes;
                        info!(
                            "Generating concatenated reference FASTA file of {} genomes ..",
                            list_of_genome_fasta_files.len()
                        );

                        (
                            Some(
                                coverm::mapping_index_maintenance::generate_concatenated_fasta_file(
                                    list_of_genome_fasta_files,
                                ),
                            ),
                            extract_genomes_and_contigs_option(
                                &m,
                                &dereplicated_genomes
                                    .clone()
                                    .iter()
                                    .map(|s| s.as_str())
                                    .collect(),
                            ),
                        )
                    }
                };

                if filter_params.doing_filtering() {
                    debug!("Mapping and filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(
                        m,
                        mapping_program,
                        &concatenated_genomes,
                        &filter_params,
                    );
                    let mut all_generators = vec![];
                    let mut indices = vec![]; // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    run_genome(
                        all_generators,
                        m,
                        &mut estimators_and_taker,
                        separator,
                        &genomes_and_contigs_option,
                    );
                } else if m.is_present("sharded") {
                    match genome_exclusion_type {
                        GenomeExclusionTypes::NoneType => {
                            run_genome(
                                get_sharded_bam_readers(
                                    m,
                                    mapping_program,
                                    &concatenated_genomes,
                                    &genome_exclusion_filter_non_type.unwrap(),
                                ),
                                m,
                                &mut estimators_and_taker,
                                separator,
                                &genomes_and_contigs_option,
                            );
                        }
                        GenomeExclusionTypes::SeparatorType => {
                            run_genome(
                                get_sharded_bam_readers(
                                    m,
                                    mapping_program,
                                    &concatenated_genomes,
                                    &genome_exclusion_filter_separator_type.unwrap(),
                                ),
                                m,
                                &mut estimators_and_taker,
                                separator,
                                &genomes_and_contigs_option,
                            );
                        }
                        GenomeExclusionTypes::GenomesAndContigsType => {
                            run_genome(
                                get_sharded_bam_readers(
                                    m,
                                    mapping_program,
                                    &concatenated_genomes,
                                    &genome_exclusion_genomes_and_contigs.unwrap(),
                                ),
                                m,
                                &mut estimators_and_taker,
                                separator,
                                &genomes_and_contigs_option,
                            );
                        }
                    }
                } else {
                    let generator_sets =
                        get_streamed_bam_readers(m, mapping_program, &concatenated_genomes);
                    let mut all_generators = vec![];
                    let mut indices = vec![]; // Prevent indices from being dropped
                    for set in generator_sets {
                        indices.push(set.index);
                        for g in set.generators {
                            all_generators.push(g)
                        }
                    }
                    run_genome(
                        all_generators,
                        m,
                        &mut estimators_and_taker,
                        separator,
                        &genomes_and_contigs_option,
                    );
                };
            }
        }
        Some("filter") => {
            let m = matches.subcommand_matches("filter").unwrap();
            if m.is_present("full-help") {
                println!("{}", filter_full_help());
                process::exit(1);
            }
            set_log_level(m, true);

            let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
            let output_bam_files: Vec<&str> = m.values_of("output-bam-files").unwrap().collect();
            if bam_files.len() != output_bam_files.len() {
                error!("The number of input BAM files must be the same as the number output");
                process::exit(1);
            }

            let filter_params = FilterParameters::generate_from_clap(m);

            let num_threads = value_t!(m.value_of("threads"), u16).unwrap();

            for (bam, output) in bam_files.iter().zip(output_bam_files.iter()) {
                let reader =
                    bam::Reader::from_path(bam).expect(&format!("Unable to find BAM file {}", bam));
                let header = bam::header::Header::from_template(reader.header());
                let mut writer =
                    bam::Writer::from_path(output, &header, rust_htslib::bam::Format::BAM)
                        .expect(&format!("Failed to write BAM file {}", output));
                writer
                    .set_threads(num_threads as usize)
                    .expect("Failed to set num threads in writer");
                let mut filtered = filter::ReferenceSortedBamFilter::new(
                    reader,
                    filter_params.flag_filters.clone(),
                    filter_params.min_aligned_length_single,
                    filter_params.min_percent_identity_single,
                    filter_params.min_aligned_percent_single,
                    filter_params.min_aligned_length_pair,
                    filter_params.min_percent_identity_pair,
                    filter_params.min_aligned_percent_pair,
                    !m.is_present("inverse"),
                );

                let mut record = bam::record::Record::new();
                while filtered
                    .read(&mut record)
                    .expect("Failure to read filtered BAM record")
                    == true
                {
                    debug!("Writing.. {:?}", record.qname());
                    writer.write(&record).expect("Failed to write BAM record");
                }
            }
        }
        Some("contig") => {
            let m = matches.subcommand_matches("contig").unwrap();
            if m.is_present("full-help") {
                println!("{}", contig_full_help());
                process::exit(1);
            }
            set_log_level(m, true);
            let print_zeros = !m.is_present("no-zeros");
            let filter_params = FilterParameters::generate_from_clap(m);
            let threads = m.value_of("threads").unwrap().parse().unwrap();

            let mut estimators_and_taker =
                EstimatorsAndTaker::generate_from_clap(m, &mut print_stream);
            estimators_and_taker =
                estimators_and_taker.print_headers(&"Contig", &mut std::io::stdout());

            if m.is_present("bam-files") {
                let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
                if filter_params.doing_filtering() {
                    let bam_readers =
                        coverm::bam_generator::generate_filtered_bam_readers_from_bam_files(
                            bam_files,
                            filter_params.flag_filters.clone(),
                            filter_params.min_aligned_length_single,
                            filter_params.min_percent_identity_single,
                            filter_params.min_aligned_percent_single,
                            filter_params.min_aligned_length_pair,
                            filter_params.min_percent_identity_pair,
                            filter_params.min_aligned_percent_pair,
                        );
                    run_contig(
                        &mut estimators_and_taker,
                        bam_readers,
                        print_zeros,
                        filter_params.flag_filters,
                        threads,
                    );
                } else if m.is_present("sharded") {
                    external_command_checker::check_for_samtools();
                    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
                    let bam_readers =
                        coverm::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                            bam_files,
                            sort_threads,
                            &NoExclusionGenomeFilter {},
                        );
                    run_contig(
                        &mut estimators_and_taker,
                        bam_readers,
                        print_zeros,
                        filter_params.flag_filters,
                        threads,
                    );
                } else {
                    let bam_readers =
                        coverm::bam_generator::generate_named_bam_readers_from_bam_files(bam_files);
                    run_contig(
                        &mut estimators_and_taker,
                        bam_readers,
                        print_zeros,
                        filter_params.flag_filters,
                        threads,
                    );
                }
            } else {
                let mapping_program = parse_mapping_program(&m);
                external_command_checker::check_for_samtools();

                if filter_params.doing_filtering() {
                    debug!("Filtering..");
                    let generator_sets = get_streamed_filtered_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &filter_params,
                    );
                    let mut all_generators = vec![];
                    let mut indices = vec![]; // Prevent indices from being dropped
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
                        filter_params.flag_filters,
                        threads,
                    );
                } else if m.is_present("sharded") {
                    let generator_sets = get_sharded_bam_readers(
                        m,
                        mapping_program,
                        &None,
                        &NoExclusionGenomeFilter {},
                    );
                    run_contig(
                        &mut estimators_and_taker,
                        generator_sets,
                        print_zeros,
                        filter_params.flag_filters,
                        threads,
                    );
                } else {
                    debug!("Not filtering..");
                    let generator_sets = get_streamed_bam_readers(m, mapping_program, &None);
                    let mut all_generators = vec![];
                    let mut indices = vec![]; // Prevent indices from being dropped
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
                        filter_params.flag_filters.clone(),
                        threads,
                    );
                }
            }
        }
        Some("make") => {
            let m = matches.subcommand_matches("make").unwrap();
            set_log_level(m, true);

            let mapping_program = parse_mapping_program(&m);
            external_command_checker::check_for_samtools();

            let output_directory = m.value_of("output-directory").unwrap();
            setup_bam_cache_directory(output_directory);
            let params = MappingParameters::generate_from_clap(&m, mapping_program, &None);
            let mut generator_sets = vec![];
            let discard_unmapped_reads = m.is_present("discard-unmapped");

            for reference_wise_params in params {
                let mut bam_readers = vec![];
                let index = setup_mapping_index(&reference_wise_params, &m, mapping_program);
                let ref_string = reference_wise_params.reference;

                for p in reference_wise_params {
                    bam_readers.push(
                        coverm::bam_generator::generate_bam_maker_generator_from_reads(
                            mapping_program,
                            match index {
                                Some(ref index) => index.index_path(),
                                None => ref_string,
                            },
                            p.read1,
                            p.read2,
                            p.read_format.clone(),
                            p.threads,
                            &generate_cached_bam_file_name(output_directory, p.reference, p.read1),
                            discard_unmapped_reads,
                            p.mapping_options,
                        ),
                    );
                }

                debug!("Finished BAM setup");
                let to_return = BamGeneratorSet {
                    generators: bam_readers,
                    index: index,
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
        }
        Some("shell-completion") => {
            let m = matches.subcommand_matches("shell-completion").unwrap();
            set_log_level(m, true);
            let mut file = std::fs::File::create(m.value_of("output-file").unwrap())
                .expect("failed to open file");
            let shell = m.value_of("shell").unwrap().parse::<Shell>().unwrap();
            info!("Generating completion script for shell {}", shell);
            app.gen_completions_to(
                "coverm", // We specify the bin name manually
                shell,    // Which shell to build completions for
                &mut file,
            ); // Where write the completions to
        }
        Some("cluster") => {
            galah::cluster_argument_parsing::run_cluster_subcommand(&matches);
        }
        _ => {
            app.print_help().unwrap();
            println!();
        }
    }
}

fn setup_mapping_index(
    reference_wise_params: &SingleReferenceMappingParameters,
    m: &clap::ArgMatches,
    mapping_program: MappingProgram,
) -> Option<Box<dyn coverm::mapping_index_maintenance::MappingIndex>> {
    match mapping_program {
        MappingProgram::BWA_MEM => Some(coverm::mapping_index_maintenance::generate_bwa_index(
            reference_wise_params.reference,
            None,
        )),
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_PB
        | MappingProgram::MINIMAP2_NO_PRESET => {
            if m.is_present("minimap2-reference-is-index") || reference_wise_params.len() == 1 {
                info!("Not pre-generating minimap2 index");
                if m.is_present("minimap2-reference-is-index") {
                    warn!(
                        "Minimap2 uses mapping parameters defined when the index was created, \
                    not parameters defined when mapping. Proceeding on the assumption that you \
                    passed the correct parameters when creating the minimap2 index."
                    );
                }
                None
            } else {
                Some(coverm::mapping_index_maintenance::generate_minimap2_index(
                    reference_wise_params.reference,
                    Some(m.value_of("threads").unwrap().parse::<usize>().unwrap()),
                    Some(m.value_of("minimap2-params").unwrap_or("")),
                    mapping_program,
                ))
            }
        }
    }
}

fn dereplicate(m: &clap::ArgMatches, genome_fasta_files: &Vec<String>) -> Vec<String> {
    info!(
        "Found {} genomes specified before dereplication",
        genome_fasta_files.len()
    );
    // Generate clusterer and check for dependencies
    let clusterer = galah::cluster_argument_parsing::generate_galah_clusterer(
        genome_fasta_files,
        &m,
        &galah_command_line_definition(),
    )
    .expect("Failed to parse galah clustering arguments correctly");

    info!("Dereplicating genome at {}% ANI ..", clusterer.ani * 100.);
    let cluster_indices = clusterer.cluster();
    info!(
        "Finished dereplication, finding {} representative genomes.",
        cluster_indices.len()
    );
    debug!("Found cluster indices: {:?}", cluster_indices);
    let reps = cluster_indices
        .iter()
        .map(|cluster| genome_fasta_files[cluster[0]].clone())
        .collect::<Vec<_>>();
    debug!("Found cluster representatives: {:?}", reps);

    if m.is_present("output-dereplication-clusters") {
        let path = m.value_of("output-dereplication-clusters").unwrap();
        info!("Writing dereplication cluster memberships to {}", path);
        let mut f =
            std::fs::File::create(path).expect("Error creating dereplication cluster output file");
        for cluster in cluster_indices.iter() {
            let rep = cluster[0];
            for member in cluster {
                writeln!(
                    f,
                    "{}\t{}",
                    genome_fasta_files[rep], genome_fasta_files[*member]
                )
                .expect("Failed to write a specific line to dereplication cluster file");
            }
        }
    }
    reps
}

fn parse_mapping_program(m: &clap::ArgMatches) -> MappingProgram {
    let mapping_program = match m.value_of("mapper") {
        Some("bwa-mem") => MappingProgram::BWA_MEM,
        Some("minimap2-sr") => MappingProgram::MINIMAP2_SR,
        Some("minimap2-ont") => MappingProgram::MINIMAP2_ONT,
        Some("minimap2-pb") => MappingProgram::MINIMAP2_PB,
        Some("minimap2-no-preset") => MappingProgram::MINIMAP2_NO_PRESET,
        None => DEFAULT_MAPPING_SOFTWARE_ENUM,
        _ => panic!(
            "Unexpected definition for --mapper: {:?}",
            m.value_of("mapper")
        ),
    };
    match mapping_program {
        MappingProgram::BWA_MEM => {
            external_command_checker::check_for_bwa();
        }
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_PB
        | MappingProgram::MINIMAP2_NO_PRESET => {
            external_command_checker::check_for_minimap2();
        }
    }
    return mapping_program;
}

struct EstimatorsAndTaker<'a> {
    estimators: Vec<CoverageEstimator>,
    taker: CoverageTakerType<'a>,
    columns_to_normalise: Vec<usize>,
    rpkm_column: Option<usize>,
    printer: CoveragePrinter,
}

fn extract_genomes_and_contigs_option(
    m: &clap::ArgMatches,
    genome_fasta_files: &Vec<&str>,
) -> Option<GenomesAndContigs> {
    match m.is_present("genome-definition") {
        true => Some(coverm::genome_parsing::read_genome_definition_file(
            m.value_of("genome-definition").unwrap(),
        )),
        false => Some(coverm::genome_parsing::read_genome_fasta_files(
            &genome_fasta_files,
        )),
    }
}

fn parse_all_genome_definitions(m: &clap::ArgMatches) -> Option<GenomesAndContigs> {
    if m.is_present("single-genome") || m.is_present("separator") {
        None
    } else if m.is_present("genome-definition") {
        Some(coverm::genome_parsing::read_genome_definition_file(
            m.value_of("genome-definition").unwrap(),
        ))
    } else {
        extract_genomes_and_contigs_option(
            m,
            &bird_tool_utils::clap_utils::parse_list_of_genome_fasta_files(&m, true)
                .expect("Failed to parse genome paths")
                .iter()
                .map(|s| s.as_str())
                .collect(),
        )
    }
}

fn parse_percentage(m: &clap::ArgMatches, parameter: &str) -> f32 {
    match m.is_present(parameter) {
        true => {
            let mut percentage = value_t!(m.value_of(parameter), f32).unwrap();
            if percentage >= 1.0 && percentage <= 100.0 {
                percentage = percentage / 100.0;
            } else if percentage < 0.0 || percentage > 100.0 {
                error!("Invalid alignment percentage: '{}'", percentage);
                process::exit(1);
            }
            info!("Using {} {}%", parameter, percentage * 100.0);
            percentage
        }
        false => 0.0,
    }
}

impl<'a> EstimatorsAndTaker<'a> {
    pub fn generate_from_clap(
        m: &clap::ArgMatches,
        stream: &'a mut std::io::Stdout,
    ) -> EstimatorsAndTaker<'a> {
        let mut estimators = vec![];
        let min_fraction_covered = parse_percentage(&m, "min-covered-fraction");
        let contig_end_exclusion = value_t!(m.value_of("contig-end-exclusion"), u64).unwrap();

        let methods: Vec<&str> = m.values_of("methods").unwrap().collect();
        let mut columns_to_normalise: Vec<usize> = vec![];

        let taker;
        let output_format = m.value_of("output-format").unwrap();
        let printer;
        let mut rpkm_column = None;

        if doing_metabat(&m) {
            estimators.push(CoverageEstimator::new_estimator_length());
            estimators.push(CoverageEstimator::new_estimator_mean(
                min_fraction_covered,
                contig_end_exclusion,
                false,
            ));
            estimators.push(CoverageEstimator::new_estimator_variance(
                min_fraction_covered,
                contig_end_exclusion,
            ));

            debug!("Cached regular coverage taker for metabat mode being used");
            taker = CoverageTakerType::new_cached_single_float_coverage_taker(estimators.len());
            printer = CoveragePrinter::MetabatAdjustedCoveragePrinter;
        } else {
            for (i, method) in methods.iter().enumerate() {
                match method {
                    &"mean" => {
                        estimators.push(CoverageEstimator::new_estimator_mean(
                            min_fraction_covered,
                            contig_end_exclusion,
                            false,
                        )); // TODO: Parameterise exclude_mismatches
                    }
                    &"coverage_histogram" => {
                        estimators.push(CoverageEstimator::new_estimator_pileup_counts(
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));
                    }
                    &"trimmed_mean" => {
                        let min = value_t!(m.value_of("trim-min"), f32).unwrap();
                        let max = value_t!(m.value_of("trim-max"), f32).unwrap();
                        if min < 0.0 || min > 1.0 || max <= min || max > 1.0 {
                            error!(
                                "error: Trim bounds must be between 0 and 1, and \
                                 min must be less than max, found {} and {}",
                                min, max
                            );
                            process::exit(1);
                        }
                        estimators.push(CoverageEstimator::new_estimator_trimmed_mean(
                            min,
                            max,
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));
                    }
                    &"covered_fraction" => {
                        estimators.push(CoverageEstimator::new_estimator_covered_fraction(
                            min_fraction_covered,
                        ));
                    }
                    &"covered_bases" => {
                        estimators.push(CoverageEstimator::new_estimator_covered_bases(
                            min_fraction_covered,
                        ));
                    }
                    &"rpkm" => {
                        if rpkm_column.is_some() {
                            error!("The RPKM column cannot be specified more than once");
                            process::exit(1);
                        }
                        rpkm_column = Some(i);
                        estimators.push(CoverageEstimator::new_estimator_rpkm(min_fraction_covered))
                    }
                    &"variance" => {
                        estimators.push(CoverageEstimator::new_estimator_variance(
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));
                    }
                    &"length" => {
                        estimators.push(CoverageEstimator::new_estimator_length());
                    }
                    &"relative_abundance" => {
                        columns_to_normalise.push(i);
                        estimators.push(CoverageEstimator::new_estimator_mean(
                            min_fraction_covered,
                            contig_end_exclusion,
                            false,
                        ));
                        // TODO: Parameterise exclude_mismatches
                    }
                    &"count" => {
                        estimators.push(CoverageEstimator::new_estimator_read_count());
                    }
                    &"reads_per_base" => {
                        estimators.push(CoverageEstimator::new_estimator_reads_per_base());
                    }
                    _ => unreachable!(),
                };
            }

            if methods.contains(&"coverage_histogram") {
                if methods.len() > 1 {
                    error!("Cannot specify the coverage_histogram method with any other coverage methods");
                    process::exit(1);
                } else {
                    debug!("Coverage histogram type coverage taker being used");
                    taker = CoverageTakerType::new_pileup_coverage_coverage_printer(stream);
                    printer = CoveragePrinter::StreamedCoveragePrinter;
                }
            } else if columns_to_normalise.len() == 0
                && rpkm_column.is_none()
                && output_format == "sparse"
            {
                debug!("Streaming regular coverage output");
                taker =
                    CoverageTakerType::new_single_float_coverage_streaming_coverage_printer(stream);
                printer = CoveragePrinter::StreamedCoveragePrinter;
            } else {
                debug!(
                    "Cached regular coverage taker with columns to normlise: {:?} and rpkm_column: {:?}",
                    columns_to_normalise, rpkm_column
                );
                taker = CoverageTakerType::new_cached_single_float_coverage_taker(estimators.len());
                printer = match output_format {
                    "sparse" => CoveragePrinter::SparseCachedCoveragePrinter,
                    "dense" => CoveragePrinter::DenseCachedCoveragePrinter {
                        entry_type: None,
                        estimator_headers: None,
                    },
                    _ => unreachable!(),
                }
            }
        }

        // Check that min-covered-fraction is being used as expected
        if min_fraction_covered != 0.0 {
            let die = |estimator_name| {
                error!(
                    "The '{}' coverage estimator cannot be used when \
                     --min-covered-fraction is > 0 as it does not calculate \
                     the covered fraction. You may wish to set the \
                     --min-covered-fraction to 0 and/or run this estimator \
                     separately.",
                    estimator_name
                );
                process::exit(1)
            };
            for e in &estimators {
                match e {
                    CoverageEstimator::ReadCountCalculator { .. } => die("counts"),
                    CoverageEstimator::ReferenceLengthCalculator { .. } => die("length"),
                    CoverageEstimator::ReadsPerBaseCalculator { .. } => die("reads_per_base"),
                    _ => {}
                }
            }
        }

        return EstimatorsAndTaker {
            estimators: estimators,
            taker: taker,
            columns_to_normalise: columns_to_normalise,
            rpkm_column: rpkm_column,
            printer: printer,
        };
    }

    pub fn print_headers(
        mut self,
        entry_type: &str,
        print_stream: &mut dyn std::io::Write,
    ) -> Self {
        let mut headers: Vec<String> = vec![];
        for e in self.estimators.iter() {
            for h in e.column_headers() {
                headers.push(h.to_string())
            }
        }
        for i in self.columns_to_normalise.iter() {
            headers[*i] = "Relative Abundance (%)".to_string();
        }
        self.printer
            .print_headers(&entry_type, headers, print_stream);
        return self;
    }
}

fn parse_separator(m: &clap::ArgMatches) -> Option<u8> {
    let single_genome = m.is_present("single-genome");
    if single_genome {
        Some("0".as_bytes()[0])
    } else if m.is_present("separator") {
        let separator_str = m.value_of("separator").unwrap().as_bytes();
        if separator_str.len() != 1 {
            eprintln!(
                "error: Separator can only be a single character, found {} ({}).",
                separator_str.len(),
                str::from_utf8(separator_str).unwrap()
            );
            process::exit(1);
        }
        Some(separator_str[0])
    } else if m.is_present("bam-files") || m.is_present("reference") {
        // Argument parsing enforces that genomes have been specified as FASTA
        // files.
        None
    } else {
        // Separator is set by CoverM and written into the generated reference
        // fasta file.
        Some(CONCATENATED_FASTA_FILE_SEPARATOR.as_bytes()[0])
    }
}

fn run_genome<
    'a,
    R: coverm::bam_generator::NamedBamReader,
    T: coverm::bam_generator::NamedBamReaderGenerator<R>,
>(
    bam_generators: Vec<T>,
    m: &clap::ArgMatches,
    estimators_and_taker: &'a mut EstimatorsAndTaker<'a>,
    separator: Option<u8>,
    genomes_and_contigs_option: &Option<GenomesAndContigs>,
) {
    let print_zeros = !m.is_present("no-zeros");
    let proper_pairs_only = m.is_present("proper-pairs-only");
    let single_genome = m.is_present("single-genome");
    let threads = m.value_of("threads").unwrap().parse().unwrap();
    let reads_mapped = match separator.is_some() || single_genome {
        true => coverm::genome::mosdepth_genome_coverage(
            bam_generators,
            separator.unwrap(),
            &mut estimators_and_taker.taker,
            print_zeros,
            &mut estimators_and_taker.estimators,
            proper_pairs_only,
            single_genome,
            threads,
        ),

        false => match genomes_and_contigs_option {
            Some(gc) => coverm::genome::mosdepth_genome_coverage_with_contig_names(
                bam_generators,
                gc,
                &mut estimators_and_taker.taker,
                print_zeros,
                proper_pairs_only,
                &mut estimators_and_taker.estimators,
                threads,
            ),
            None => unreachable!(),
        },
    };

    debug!("Finalising printing ..");
    estimators_and_taker.printer.finalise_printing(
        &estimators_and_taker.taker,
        &mut std::io::stdout(),
        Some(&reads_mapped),
        &estimators_and_taker.columns_to_normalise,
        estimators_and_taker.rpkm_column,
    );
}

fn doing_metabat(m: &clap::ArgMatches) -> bool {
    match m.subcommand_name() {
        Some("contig") | None => {
            if !m.is_present("methods") {
                return false;
            }
            let methods: Vec<&str> = m.values_of("methods").unwrap().collect();
            if methods.contains(&"metabat") {
                if methods.len() > 1 {
                    error!("Cannot specify the metabat method with any other coverage methods");
                    process::exit(1);
                } else {
                    return true;
                }
            }
            return false;
        }
        _ => {
            debug!("Not running in contig mode so cannot be in metabat mode");
            return false;
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
            min_aligned_length_single: match m.is_present("min-read-aligned-length") {
                true => value_t!(m.value_of("min-read-aligned-length"), u32).unwrap(),
                false => 0,
            },
            min_percent_identity_single: parse_percentage(&m, "min-read-percent-identity"),
            min_aligned_percent_single: parse_percentage(&m, "min-read-aligned-percent"),
            min_aligned_length_pair: match m.is_present("min-read-aligned-length-pair") {
                true => value_t!(m.value_of("min-read-aligned-length-pair"), u32).unwrap(),
                false => 0,
            },
            min_percent_identity_pair: parse_percentage(&m, "min-read-percent-identity-pair"),
            min_aligned_percent_pair: parse_percentage(&m, "min-read-aligned-percent-pair"),
        };
        if doing_metabat(&m) {
            debug!(
                "Setting single read percent identity threshold at 0.97 for \
                 MetaBAT adjusted coverage."
            );
            // we use >= where metabat uses >. Gah.
            f.min_percent_identity_single = 0.97001;
            f.flag_filters.include_improper_pairs = true;
            f.flag_filters.include_supplementary = true;
            f.flag_filters.include_secondary = true;
        }
        debug!("Filter parameters set as {:?}", f);
        return f;
    }

    pub fn doing_filtering(&self) -> bool {
        return self.min_percent_identity_single > 0.0
            || self.min_percent_identity_pair > 0.0
            || self.min_aligned_percent_single > 0.0
            || self.min_aligned_percent_pair > 0.0
            || self.min_aligned_length_single > 0
            || self.min_aligned_length_pair > 0;
    }
}

fn get_sharded_bam_readers<'a, 'b, T>(
    m: &'a clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &'a Option<NamedTempFile>,
    genome_exclusion: &'b T,
) -> Vec<ShardedBamReaderGenerator<'b, T>>
where
    T: GenomeExclusion,
{
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");
    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
    let params = MappingParameters::generate_from_clap(&m, mapping_program, &reference_tempfile);
    let mut bam_readers = vec![];
    let mut concatenated_reference_name: Option<String> = None;
    let mut concatenated_read_names: Option<String> = None;

    for reference_wise_params in params {
        let index = setup_mapping_index(&reference_wise_params, &m, mapping_program);

        let reference = reference_wise_params.reference;
        let reference_name = std::path::Path::new(reference)
            .file_name()
            .expect("Unable to convert reference to file name")
            .to_str()
            .expect("Unable to covert file name into str")
            .to_string();
        concatenated_reference_name = match concatenated_reference_name {
            Some(prev) => Some(format!("{}|{}", prev, reference_name)),
            None => Some(reference_name),
        };
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.is_present("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.value_of("bam-file-cache-directory").unwrap(),
                        match reference_tempfile {
                            Some(_) => CONCATENATED_REFERENCE_CACHE_STEM,
                            None => reference,
                        },
                        naming_readset,
                    );
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        for p in reference_wise_params {
            bam_readers.push(
                coverm::shard_bam_reader::generate_named_sharded_bam_readers_from_reads(
                    mapping_program,
                    match index {
                        Some(ref index) => index.index_path(),
                        None => reference,
                    },
                    p.read1,
                    p.read2,
                    p.read_format.clone(),
                    p.threads,
                    bam_file_cache(p.read1).as_ref().map(String::as_ref),
                    discard_unmapped,
                    p.mapping_options,
                ),
            );
            let name = &std::path::Path::new(p.read1)
                .file_name()
                .expect("Unable to convert read1 name to file name")
                .to_str()
                .expect("Unable to covert file name into str")
                .to_string();
            concatenated_read_names = match concatenated_read_names {
                Some(prev) => Some(format!("{}|{}", prev, name)),
                None => Some(name.to_string()),
            };
        }

        debug!("Finished BAM setup");
    }
    let gen = ShardedBamReaderGenerator {
        stoit_name: format!(
            "{}/{}",
            concatenated_reference_name.unwrap(),
            concatenated_read_names.unwrap()
        ),
        read_sorted_bam_readers: bam_readers,
        sort_threads: sort_threads,
        genome_exclusion: genome_exclusion,
    };
    return vec![gen];
}

fn get_streamed_bam_readers<'a>(
    m: &'a clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &'a Option<NamedTempFile>,
) -> Vec<BamGeneratorSet<StreamingNamedBamReaderGenerator>> {
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params = MappingParameters::generate_from_clap(&m, mapping_program, &reference_tempfile);
    let mut generator_set = vec![];
    for reference_wise_params in params {
        let mut bam_readers = vec![];
        let index = setup_mapping_index(&reference_wise_params, &m, mapping_program);

        let reference = reference_wise_params.reference;
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.is_present("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.value_of("bam-file-cache-directory").unwrap(),
                        match reference_tempfile {
                            Some(_) => CONCATENATED_REFERENCE_CACHE_STEM,
                            None => reference,
                        },
                        naming_readset,
                    );
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        for p in reference_wise_params {
            bam_readers.push(
                coverm::bam_generator::generate_named_bam_readers_from_reads(
                    mapping_program,
                    match index {
                        Some(ref index) => index.index_path(),
                        None => reference,
                    },
                    p.read1,
                    p.read2,
                    p.read_format.clone(),
                    p.threads,
                    bam_file_cache(p.read1).as_ref().map(String::as_ref),
                    discard_unmapped,
                    p.mapping_options,
                    reference_tempfile.is_none(),
                ),
            );
        }

        debug!("Finished BAM setup");
        let to_return = BamGeneratorSet {
            generators: bam_readers,
            index: index,
        };
        generator_set.push(to_return);
    }
    return generator_set;
}

fn generate_cached_bam_file_name(directory: &str, reference: &str, read1_path: &str) -> String {
    debug!(
        "Constructing BAM file cache name in directory {}, reference {}, read1_path {}",
        directory, reference, read1_path
    );
    std::path::Path::new(directory)
        .to_str()
        .expect("Unable to covert bam-file-cache-directory name into str")
        .to_string()
        + "/"
        + &std::path::Path::new(reference)
            .file_name()
            .expect("Unable to convert reference to file name")
            .to_str()
            .expect("Unable to covert file name into str")
            .to_string()
        + "."
        + &std::path::Path::new(read1_path)
            .file_name()
            .expect("Unable to convert read1 name to file name")
            .to_str()
            .expect("Unable to covert file name into str")
            .to_string()
        + ".bam"
}

fn setup_bam_cache_directory(cache_directory: &str) {
    let path = std::path::Path::new(cache_directory);
    if path.is_dir() {
        if path
            .metadata()
            .expect("Unable to read metadata for cache directory")
            .permissions()
            .readonly()
        {
            error!(
                "Cache directory {} does not appear to be writeable, not continuing",
                cache_directory
            );
            process::exit(1);
        } else {
            info!(
                "Writing BAM files to already existing directory {}",
                cache_directory
            )
        }
    } else {
        match path.parent() {
            Some(parent) => {
                let parent2 = match parent == std::path::Path::new("") {
                    true => std::path::Path::new("."),
                    false => parent,
                };
                if parent2
                    .canonicalize()
                    .expect(&format!(
                        "Unable to canonicalize parent of cache directory {}",
                        cache_directory
                    ))
                    .is_dir()
                {
                    if parent2
                        .metadata()
                        .expect(&format!(
                            "Unable to get metadata for parent of cache directory {}",
                            cache_directory
                        ))
                        .permissions()
                        .readonly()
                    {
                        error!(
                            "The parent directory of the (currently non-existent) \
                             cache directory {} is not writeable, not continuing",
                            cache_directory
                        );
                        process::exit(1);
                    } else {
                        info!("Creating cache directory {}", cache_directory);
                        std::fs::create_dir(path).expect("Unable to create cache directory");
                    }
                } else {
                    error!(
                        "The parent directory of the cache directory {} does not \
                         yet exist, so not creating that cache directory, and not continuing.",
                        cache_directory
                    );
                    process::exit(1);
                }
            }
            None => {
                error!("Cannot create root directory {}", cache_directory);
                process::exit(1);
            }
        }
    }
    // Test writing a tempfile to the directory, to test it actually is
    // writeable.
    let tf_result = tempfile::tempfile_in(path);
    if tf_result.is_err() {
        error!(
            "Failed to create test file in bam cache directory: {}",
            tf_result.err().unwrap()
        );
        process::exit(1);
    }
}

fn get_streamed_filtered_bam_readers(
    m: &clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &Option<NamedTempFile>,
    filter_params: &FilterParameters,
) -> Vec<BamGeneratorSet<StreamingFilteredNamedBamReaderGenerator>> {
    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params = MappingParameters::generate_from_clap(&m, mapping_program, &reference_tempfile);
    let mut generator_set = vec![];
    for reference_wise_params in params {
        let mut bam_readers = vec![];
        let index = setup_mapping_index(&reference_wise_params, &m, mapping_program);

        let reference = reference_wise_params.reference;
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.is_present("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.value_of("bam-file-cache-directory").unwrap(),
                        match reference_tempfile {
                            Some(_) => CONCATENATED_REFERENCE_CACHE_STEM,
                            None => reference,
                        },
                        naming_readset,
                    );
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        for p in reference_wise_params {
            bam_readers.push(
                coverm::bam_generator::generate_filtered_named_bam_readers_from_reads(
                    mapping_program,
                    match index {
                        Some(ref index) => index.index_path(),
                        None => reference,
                    },
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
                    p.mapping_options,
                    discard_unmapped,
                    reference_tempfile.is_none(),
                ),
            );
        }

        debug!("Finished BAM setup");
        let to_return = BamGeneratorSet {
            generators: bam_readers,
            index: index,
        };
        generator_set.push(to_return);
    }
    return generator_set;
}

fn run_contig<
    'a,
    R: coverm::bam_generator::NamedBamReader,
    T: coverm::bam_generator::NamedBamReaderGenerator<R>,
>(
    estimators_and_taker: &'a mut EstimatorsAndTaker<'a>,
    bam_readers: Vec<T>,
    print_zeros: bool,
    flag_filters: FlagFilter,
    threads: usize,
) {
    let reads_mapped = coverm::contig::contig_coverage(
        bam_readers,
        &mut estimators_and_taker.taker,
        &mut estimators_and_taker.estimators,
        print_zeros,
        flag_filters,
        threads,
    );

    debug!("Finalising printing ..");

    estimators_and_taker.printer.finalise_printing(
        &estimators_and_taker.taker,
        &mut std::io::stdout(),
        Some(&reads_mapped),
        &estimators_and_taker.columns_to_normalise,
        estimators_and_taker.rpkm_column,
    );
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
            builder.parse_filters(&env::var("RUST_LOG").unwrap());
        }
        if builder.try_init().is_err() {
            panic!("Failed to set log level - has it been specified multiple times?")
        }
    }
    if is_last {
        info!("CoverM version {}", crate_version!());
    }
}
