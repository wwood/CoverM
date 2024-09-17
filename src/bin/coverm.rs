extern crate coverm;
use coverm::bam_generator::*;
use coverm::cli::*;
use coverm::coverage_printer::*;
use coverm::coverage_takers::*;
use coverm::external_command_checker;
use coverm::filter;
use coverm::genome_exclusion::*;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::mapping_index_maintenance::check_reference_existence;
use coverm::mapping_parameters::*;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::shard_bam_reader::*;
use coverm::FlagFilter;
use coverm::OutputWriter;
use coverm::CONCATENATED_FASTA_FILE_SEPARATOR;

extern crate galah;
use galah::ClusterDistanceFinder;

extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;

use std::collections::HashSet;
use std::env;
use std::process;
use std::str;

extern crate clap;
use clap::*;

extern crate clap_complete;
use clap_complete::{generate, Shell};

#[macro_use]
extern crate log;

extern crate tempfile;
use tempfile::NamedTempFile;

extern crate bird_tool_utils;
use bird_tool_utils::clap_utils::set_log_level as set_log_level_bird_tool_utils;

const CONCATENATED_REFERENCE_CACHE_STEM: &str = "coverm-genome";
const DEFAULT_MAPPING_SOFTWARE_ENUM: MappingProgram = MappingProgram::MINIMAP2_SR;

fn main() {
    let mut app = build_cli();
    let matches = app.clone().get_matches();
    set_log_level(&matches, false);
    let mut print_stream;

    match matches.subcommand_name() {
        Some("genome") => {
            let m = matches.subcommand_matches("genome").unwrap();
            bird_tool_utils::clap_utils::print_full_help_if_needed(m, genome_full_help());
            set_log_level(m, true);
            print_stream = OutputWriter::generate(m.get_one::<String>("output-file").map(|x| &**x));

            let genome_names_content: Vec<u8>;

            let mut estimators_and_taker =
                EstimatorsAndTaker::generate_from_clap(m, print_stream.clone());
            estimators_and_taker =
                estimators_and_taker.print_headers("Genome", print_stream.clone());
            let filter_params = FilterParameters::generate_from_clap(m);
            let separator = parse_separator(m);

            let genomes_and_contigs_option_predereplication = if m.get_flag("sharded")
                && !m.contains_id("separator")
                && !m.get_flag("dereplicate")
                && !m.get_flag("single-genome")
            {
                parse_all_genome_definitions(m)
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
                Separator,
                None,
                GenomesAndContigs,
            }
            let genome_exclusion_type = {
                if m.get_flag("sharded") {
                    if m.contains_id("exclude-genomes-from-deshard") {
                        let filename = m.get_one::<String>("exclude-genomes-from-deshard").unwrap();
                        genome_names_content = std::fs::read(filename).unwrap_or_else(|_| {
                            panic!(
                                "Failed to open file '{}' containing list of excluded genomes",
                                filename
                            )
                        });
                        let mut genome_names_hash: HashSet<&[u8]> = HashSet::new();
                        for n in genome_names_content.split(|s| *s == b"\n"[0]) {
                            if n != b"" {
                                genome_names_hash.insert(n);
                            }
                        }
                        if genome_names_hash.is_empty() {
                            warn!("No genomes read in that are to be excluded from desharding process");
                            genome_exclusion_filter_non_type = Some(NoExclusionGenomeFilter {});
                            GenomeExclusionTypes::None
                        } else {
                            info!("Read in {} distinct genomes to exclude from desharding process e.g. '{}'",
                                  genome_names_hash.len(),
                                  std::str::from_utf8(genome_names_hash.iter().next().unwrap())
                                  .unwrap());
                            if let Some(s) = separator {
                                genome_exclusion_filter_separator_type =
                                    Some(SeparatorGenomeExclusionFilter {
                                        split_char: s,
                                        excluded_genomes: genome_names_hash,
                                    });
                                GenomeExclusionTypes::Separator
                            } else {
                                match genomes_and_contigs_option_predereplication {
                                    Some(ref gc) => {
                                        genome_exclusion_genomes_and_contigs =
                                            Some(GenomesAndContigsExclusionFilter {
                                                genomes_and_contigs: gc,
                                                excluded_genomes: genome_names_hash,
                                            });
                                        GenomeExclusionTypes::GenomesAndContigs
                                    }
                                    None => unreachable!(),
                                }
                            }
                        }
                    } else {
                        debug!("Not excluding any genomes during the deshard process");
                        genome_exclusion_filter_non_type = Some(NoExclusionGenomeFilter {});
                        GenomeExclusionTypes::None
                    }
                } else {
                    genome_exclusion_filter_non_type = Some(NoExclusionGenomeFilter {});
                    GenomeExclusionTypes::None
                }
            };

            if m.contains_id("bam-files") {
                let bam_files: Vec<&str> = m
                    .get_many::<String>("bam-files")
                    .unwrap()
                    .map(|s| &**s)
                    .collect();

                // Associate genomes and contig names, if required
                let genomes_and_contigs_option = parse_all_genome_definitions(m);

                if filter_params.doing_filtering() {
                    run_genome(
                        coverm::bam_generator::generate_filtered_bam_readers_from_bam_files(
                            bam_files,
                            filter_params.flag_filters,
                            filter_params.min_aligned_length_single,
                            filter_params.min_percent_identity_single,
                            filter_params.min_aligned_percent_single,
                            filter_params.min_mapq,
                            filter_params.min_aligned_length_pair,
                            filter_params.min_percent_identity_pair,
                            filter_params.min_aligned_percent_pair,
                        ),
                        m,
                        &mut estimators_and_taker,
                        separator,
                        &genomes_and_contigs_option,
                        &mut print_stream,
                    );
                } else if m.get_flag("sharded") {
                    external_command_checker::check_for_samtools();
                    let sort_threads = *m.get_one::<u16>("threads").unwrap();
                    // Seems crazy, but I cannot work out how to make this more
                    // DRY, without making GenomeExclusion into an enum.
                    match genome_exclusion_type {
                        GenomeExclusionTypes::None => {
                            run_genome(
                                coverm::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                                    bam_files,
                                    sort_threads,
                                    &genome_exclusion_filter_non_type.unwrap()),
                                m,
                                &mut estimators_and_taker,
                                separator,
                                &genomes_and_contigs_option,
                                &mut print_stream,);
                        }
                        GenomeExclusionTypes::Separator => {
                            run_genome(
                                coverm::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                                    bam_files,
                                    sort_threads,
                                    &genome_exclusion_filter_separator_type.unwrap()),
                                m,
                                &mut estimators_and_taker,
                                separator,
                                &genomes_and_contigs_option,
                                &mut print_stream,);
                        }
                        GenomeExclusionTypes::GenomesAndContigs => {
                            run_genome(
                                coverm::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                                    bam_files,
                                    sort_threads,
                                    &genome_exclusion_genomes_and_contigs.unwrap()),
                                m,
                                &mut estimators_and_taker,
                                separator,
                                &genomes_and_contigs_option,
                                &mut print_stream,);
                        }
                    }
                } else {
                    run_genome(
                        coverm::bam_generator::generate_named_bam_readers_from_bam_files(bam_files),
                        m,
                        &mut estimators_and_taker,
                        separator,
                        &genomes_and_contigs_option,
                        &mut print_stream,
                    );
                }
            } else {
                let mapping_program = parse_mapping_program(m);
                external_command_checker::check_for_samtools();

                // If genomes defined by file, then potentially dereplicate. If
                // no reference is specified, make a new concatenated one.

                // Get the GenomesAndContigs from the dereplicated genomes or
                // list of genome fasta files.

                let genome_fasta_files_opt = {
                    match bird_tool_utils::clap_utils::parse_list_of_genome_fasta_files(m, false) {
                        Ok(paths) => {
                            if paths.is_empty() {
                                error!(
                                    "Genome paths were described, but ultimately none were found"
                                );
                                process::exit(1);
                            }
                            if m.contains_id("checkm-tab-table") || m.contains_id("genome-info") {
                                let genomes_after_filtering =
                                    galah::cluster_argument_parsing::filter_genomes_through_checkm(
                                        &paths,
                                        m,
                                        &coverm::cli::COVERM_CLUSTER_COMMAND_DEFINITION,
                                    )
                                    .expect("Error parsing CheckM-related options");
                                info!(
                                    "After filtering by CheckM, {} genomes remained",
                                    genomes_after_filtering.len()
                                );
                                if genomes_after_filtering.is_empty() {
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

                let (concatenated_genomes, genomes_and_contigs_option) =
                    match m.contains_id("reference") {
                        true => {
                            check_reference_existence(
                                m.get_one::<String>("reference").unwrap(),
                                &mapping_program,
                            );
                            match genome_fasta_files_opt {
                                Some(genome_paths) => (
                                    None,
                                    extract_genomes_and_contigs_option(
                                        m,
                                        &genome_paths.iter().map(|s| s.as_str()).collect(),
                                    ),
                                ),
                                None => match m.get_one::<String>("genome-definition") {
                                    Some(definition_path) => (
                                        None,
                                        Some(coverm::genome_parsing::read_genome_definition_file(
                                            definition_path,
                                        )),
                                    ),
                                    None => (None, None),
                                },
                            }
                        }
                        false => {
                            // Dereplicate if required
                            let dereplicated_genomes: Vec<String> = if m.get_flag("dereplicate") {
                                dereplicate(m, &genome_fasta_files_opt.unwrap())
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
                            None,
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
                        &mut print_stream,
                    );
                } else if m.get_flag("sharded") {
                    match genome_exclusion_type {
                        GenomeExclusionTypes::None => {
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
                                &mut print_stream,
                            );
                        }
                        GenomeExclusionTypes::Separator => {
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
                                &mut print_stream,
                            );
                        }
                        GenomeExclusionTypes::GenomesAndContigs => {
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
                                &mut print_stream,
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
                        &mut print_stream,
                    );
                };
            }
        }
        Some("filter") => {
            let m = matches.subcommand_matches("filter").unwrap();
            bird_tool_utils::clap_utils::print_full_help_if_needed(m, filter_full_help());
            set_log_level(m, true);

            let bam_files: Vec<&str> = m
                .get_many::<String>("bam-files")
                .unwrap()
                .map(|x| &**x)
                .collect();
            let output_bam_files: Vec<&str> = m
                .get_many::<String>("output-bam-files")
                .unwrap()
                .map(|x| &**x)
                .collect();
            if bam_files.len() != output_bam_files.len() {
                error!("The number of input BAM files must be the same as the number output");
                process::exit(1);
            }

            let filter_params = FilterParameters::generate_from_clap(m);

            let num_threads: u16 = *m.get_one::<u16>("threads").unwrap();

            for (bam, output) in bam_files.iter().zip(output_bam_files.iter()) {
                let reader = bam::Reader::from_path(bam)
                    .unwrap_or_else(|_| panic!("Unable to find BAM file {}", bam));
                let header = bam::header::Header::from_template(reader.header());
                let mut writer =
                    bam::Writer::from_path(output, &header, rust_htslib::bam::Format::Bam)
                        .unwrap_or_else(|_| panic!("Failed to write BAM file {}", output));
                writer
                    .set_threads(num_threads as usize)
                    .expect("Failed to set num threads in writer");
                let mut filtered = filter::ReferenceSortedBamFilter::new(
                    reader,
                    filter_params.flag_filters.clone(),
                    filter_params.min_aligned_length_single,
                    filter_params.min_percent_identity_single,
                    filter_params.min_aligned_percent_single,
                    filter_params.min_mapq,
                    filter_params.min_aligned_length_pair,
                    filter_params.min_percent_identity_pair,
                    filter_params.min_aligned_percent_pair,
                    !m.get_flag("inverse"),
                );

                let mut record = bam::record::Record::new();
                loop {
                    match filtered.read(&mut record) {
                        None => {
                            break;
                        }
                        Some(Ok(())) => {}
                        Some(e) => {
                            panic!("Failure to read filtered BAM record: {:?}", e)
                        }
                    }

                    debug!("Writing.. {:?}", record.qname());
                    writer.write(&record).expect("Failed to write BAM record");
                }
            }
        }
        Some("contig") => {
            let m = matches.subcommand_matches("contig").unwrap();
            bird_tool_utils::clap_utils::print_full_help_if_needed(m, contig_full_help());
            set_log_level(m, true);
            let print_zeros = !m.get_flag("no-zeros");

            // Add metabat filtering params since we are running contig
            let mut filter_params1 = FilterParameters::generate_from_clap(m);
            filter_params1.add_metabat_filtering_if_required(m);
            let filter_params = filter_params1;

            let threads = *m.get_one::<u16>("threads").unwrap();
            print_stream = OutputWriter::generate(m.get_one::<String>("output-file").map(|x| &**x));

            let mut estimators_and_taker =
                EstimatorsAndTaker::generate_from_clap(m, print_stream.clone());
            estimators_and_taker =
                estimators_and_taker.print_headers("Contig", print_stream.clone());

            if m.contains_id("bam-files") {
                let bam_files: Vec<&str> = m
                    .get_many::<String>("bam-files")
                    .unwrap()
                    .map(|x| &**x)
                    .collect();
                if filter_params.doing_filtering() {
                    let bam_readers =
                        coverm::bam_generator::generate_filtered_bam_readers_from_bam_files(
                            bam_files,
                            filter_params.flag_filters.clone(),
                            filter_params.min_aligned_length_single,
                            filter_params.min_percent_identity_single,
                            filter_params.min_aligned_percent_single,
                            filter_params.min_mapq,
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
                        &mut print_stream,
                    );
                } else if m.get_flag("sharded") {
                    external_command_checker::check_for_samtools();
                    let bam_readers =
                        coverm::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                            bam_files,
                            threads,
                            &NoExclusionGenomeFilter {},
                        );
                    run_contig(
                        &mut estimators_and_taker,
                        bam_readers,
                        print_zeros,
                        filter_params.flag_filters,
                        threads,
                        &mut print_stream,
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
                        &mut print_stream,
                    );
                }
            } else {
                let mapping_program = parse_mapping_program(m);
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
                        &mut print_stream,
                    );
                } else if m.get_flag("sharded") {
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
                        &mut print_stream,
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
                        filter_params.flag_filters,
                        threads,
                        &mut print_stream,
                    );
                }
            }
        }
        Some("make") => {
            let m = matches.subcommand_matches("make").unwrap();
            bird_tool_utils::clap_utils::print_full_help_if_needed(m, make_full_help());
            set_log_level(m, true);

            let mapping_program = parse_mapping_program(m);
            external_command_checker::check_for_samtools();

            let output_directory = m.get_one::<String>("output-directory").unwrap();
            setup_bam_cache_directory(output_directory);
            let params = MappingParameters::generate_from_clap(m, mapping_program, &None);
            let mut generator_sets = vec![];
            let discard_unmapped_reads = m.get_flag("discard-unmapped");

            for reference_wise_params in params {
                let mut bam_readers = vec![];
                let index = setup_mapping_index(&reference_wise_params, m, mapping_program);

                // Ensure there is no duplication in output files e.g.
                // https://github.com/wwood/CoverM/issues/128
                let mut unique_names = HashSet::new();
                for p in reference_wise_params {
                    let name =
                        generate_cached_bam_file_name(output_directory, p.reference, p.read1);
                    bam_readers.push(
                        coverm::bam_generator::generate_bam_maker_generator_from_reads(
                            mapping_program,
                            index.as_ref(),
                            p.read1,
                            p.read2,
                            p.read_format.clone(),
                            p.threads,
                            &name.clone(),
                            discard_unmapped_reads,
                            p.mapping_options,
                        ),
                    );
                    if !unique_names.insert(name.clone()) {
                        error!("Duplicate output file name: {}", name);
                        std::process::exit(1);
                    }
                }

                debug!("Finished BAM setup");
                let to_return = BamGeneratorSet {
                    generators: bam_readers,
                    index,
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
            let mut file = std::fs::File::create(m.get_one::<String>("output-file").unwrap())
                .expect("failed to open output file");

            if let Some(generator) = m.get_one::<Shell>("shell").copied() {
                let mut cmd = build_cli();
                info!("Generating completion script for shell {}", generator);
                let name = cmd.get_name().to_string();
                generate(generator, &mut cmd, name, &mut file);
            }
        }
        Some("cluster") => {
            galah::cluster_argument_parsing::run_cluster_subcommand(
                &matches,
                "CoverM",
                crate_version!(),
            );
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
) -> Box<dyn coverm::mapping_index_maintenance::MappingIndex> {
    match mapping_program {
        MappingProgram::BWA_MEM | MappingProgram::BWA_MEM2 => {
            coverm::mapping_index_maintenance::generate_bwa_index(
                reference_wise_params.reference,
                None,
                mapping_program,
            )
        }
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_HIFI
        | MappingProgram::MINIMAP2_PB
        | MappingProgram::MINIMAP2_NO_PRESET => {
            if m.get_flag("minimap2-reference-is-index") || reference_wise_params.len() == 1 {
                info!("Not pre-generating minimap2 index");
                if m.get_flag("minimap2-reference-is-index") {
                    warn!(
                        "Minimap2 uses mapping parameters defined when the index was created, \
                    not parameters defined when mapping. Proceeding on the assumption that you \
                    passed the correct parameters when creating the minimap2 index."
                    );
                }
                Box::new(coverm::mapping_index_maintenance::VanillaIndexStruct::new(
                    reference_wise_params.reference,
                ))
            } else {
                coverm::mapping_index_maintenance::generate_minimap2_index(
                    reference_wise_params.reference,
                    Some(*m.get_one::<u16>("threads").unwrap()),
                    Some(
                        m.get_one::<String>("minimap2-params")
                            .unwrap_or(&"".to_string()),
                    ),
                    mapping_program,
                )
            }
        }
        MappingProgram::STROBEALIGN => {
            // Indexing once for a batch of readsets is not yet supported for strobealign
            info!("Not pre-generating strobealign index");
            if m.get_flag("strobealign-use-index") {
                warn!(
                    "Strobealign uses mapping parameters defined when the index was created, \
                not parameters defined when mapping. Proceeding on the assumption that you \
                passed the correct parameters when creating the strobealign index."
                );

                Box::new(
                    coverm::mapping_index_maintenance::PregeneratedStrobealignIndexStruct::new(
                        reference_wise_params.reference,
                    ),
                )
            } else {
                Box::new(coverm::mapping_index_maintenance::VanillaIndexStruct::new(
                    reference_wise_params.reference,
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
        m,
        &coverm::cli::COVERM_CLUSTER_COMMAND_DEFINITION,
    )
    .expect("Failed to parse galah clustering arguments correctly");

    let cluster_outputs = galah::cluster_argument_parsing::setup_galah_outputs(
        m,
        &coverm::cli::COVERM_CLUSTER_COMMAND_DEFINITION,
    );

    info!(
        "Dereplicating genomes at {}% ANI ..",
        clusterer.clusterer.get_ani_threshold()
    );
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

    galah::cluster_argument_parsing::write_galah_outputs(
        cluster_outputs,
        &cluster_indices,
        &clusterer.genome_fasta_paths,
    );

    reps
}

fn parse_mapping_program(m: &clap::ArgMatches) -> MappingProgram {
    let mapping_program = match m.get_one::<String>("mapper").map(|x| &**x) {
        Some("bwa-mem") => MappingProgram::BWA_MEM,
        Some("bwa-mem2") => MappingProgram::BWA_MEM2,
        Some("minimap2-sr") => MappingProgram::MINIMAP2_SR,
        Some("minimap2-ont") => MappingProgram::MINIMAP2_ONT,
        Some("minimap2-pb") => MappingProgram::MINIMAP2_PB,
        Some("minimap2-hifi") => MappingProgram::MINIMAP2_HIFI,
        Some("minimap2-no-preset") => MappingProgram::MINIMAP2_NO_PRESET,
        Some("strobealign") => MappingProgram::STROBEALIGN,
        None => DEFAULT_MAPPING_SOFTWARE_ENUM,
        _ => panic!(
            "Unexpected definition for --mapper: {:?}",
            m.get_one::<String>("mapper")
        ),
    };
    match mapping_program {
        MappingProgram::BWA_MEM => {
            external_command_checker::check_for_bwa();
        }
        MappingProgram::BWA_MEM2 => {
            external_command_checker::check_for_bwa_mem2();
        }
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_HIFI
        | MappingProgram::MINIMAP2_PB
        | MappingProgram::MINIMAP2_NO_PRESET => {
            external_command_checker::check_for_minimap2();
        }
        MappingProgram::STROBEALIGN => {
            external_command_checker::check_for_strobealign();
        }
    }
    mapping_program
}

struct EstimatorsAndTaker {
    estimators: Vec<CoverageEstimator>,
    taker: CoverageTakerType,
    columns_to_normalise: Vec<usize>,
    rpkm_column: Option<usize>,
    tpm_column: Option<usize>,
    printer: CoveragePrinter,
}

fn extract_genomes_and_contigs_option(
    m: &clap::ArgMatches,
    genome_fasta_files: &Vec<&str>,
) -> Option<GenomesAndContigs> {
    match m.contains_id("genome-definition") {
        true => Some(coverm::genome_parsing::read_genome_definition_file(
            m.get_one::<String>("genome-definition").unwrap(),
        )),
        false => Some(coverm::genome_parsing::read_genome_fasta_files(
            genome_fasta_files,
            m.get_flag("use-full-contig-names"),
        )),
    }
}

fn parse_all_genome_definitions(m: &clap::ArgMatches) -> Option<GenomesAndContigs> {
    if m.get_flag("single-genome") || m.contains_id("separator") {
        None
    } else if m.contains_id("genome-definition") {
        Some(coverm::genome_parsing::read_genome_definition_file(
            m.get_one::<String>("genome-definition").unwrap(),
        ))
    } else {
        extract_genomes_and_contigs_option(
            m,
            &bird_tool_utils::clap_utils::parse_list_of_genome_fasta_files(m, true)
                .expect("Failed to parse genome paths")
                .iter()
                .map(|s| s.as_str())
                .collect(),
        )
    }
}

fn parse_percentage(m: &clap::ArgMatches, parameter: &str) -> f32 {
    if m.contains_id(parameter) {
        let mut percentage: f32 = *m.get_one::<f32>(parameter).unwrap_or(&0.0);
        if (1.0..=100.0).contains(&percentage) {
            percentage /= 100.0;
        } else if !(0.0..=100.0).contains(&percentage) {
            error!("Invalid alignment percentage: '{}'", percentage);
            process::exit(1);
        }
        if m.value_source(parameter) == Some(clap::parser::ValueSource::CommandLine) {
            info!("Using {} {}%", parameter, percentage * 100.0);
        }
        percentage
    } else {
        0.0
    }
}

impl EstimatorsAndTaker {
    pub fn generate_from_clap(m: &clap::ArgMatches, stream: OutputWriter) -> EstimatorsAndTaker {
        let mut estimators = vec![];
        let min_fraction_covered = parse_percentage(m, "min-covered-fraction");
        let contig_end_exclusion = *m.get_one::<u64>("contig-end-exclusion").unwrap();

        let methods: Vec<&str> = m
            .get_many::<String>("methods")
            .unwrap()
            .map(|x| &**x)
            .collect();
        let mut columns_to_normalise: Vec<usize> = vec![];

        let taker;
        let output_format = m.get_one::<String>("output-format").unwrap().as_str();
        let printer;
        let mut rpkm_column = None;
        let mut tpm_column = None;

        if doing_metabat(m) {
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
                match *method {
                    "mean" => {
                        estimators.push(CoverageEstimator::new_estimator_mean(
                            min_fraction_covered,
                            contig_end_exclusion,
                            false,
                        )); // TODO: Parameterise exclude_mismatches
                    }
                    "coverage_histogram" => {
                        estimators.push(CoverageEstimator::new_estimator_pileup_counts(
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));
                    }
                    "trimmed_mean" => {
                        let min = parse_percentage(m, "trim-min");
                        let max = parse_percentage(m, "trim-max");
                        estimators.push(CoverageEstimator::new_estimator_trimmed_mean(
                            min,
                            max,
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));
                    }
                    "covered_fraction" => {
                        estimators.push(CoverageEstimator::new_estimator_covered_fraction(
                            min_fraction_covered,
                        ));
                    }
                    "covered_bases" => {
                        estimators.push(CoverageEstimator::new_estimator_covered_bases(
                            min_fraction_covered,
                        ));
                    }
                    "rpkm" => {
                        if rpkm_column.is_some() {
                            error!("The RPKM column cannot be specified more than once");
                            process::exit(1);
                        }
                        rpkm_column = Some(i);
                        estimators.push(CoverageEstimator::new_estimator_rpkm(min_fraction_covered))
                    }
                    "tpm" => {
                        if tpm_column.is_some() {
                            error!("The TPM column cannot be specified more than once");
                            process::exit(1);
                        }
                        tpm_column = Some(i);
                        estimators.push(CoverageEstimator::new_estimator_tpm(min_fraction_covered))
                    }
                    "variance" => {
                        estimators.push(CoverageEstimator::new_estimator_variance(
                            min_fraction_covered,
                            contig_end_exclusion,
                        ));
                    }
                    "length" => {
                        estimators.push(CoverageEstimator::new_estimator_length());
                    }
                    "relative_abundance" => {
                        columns_to_normalise.push(i);
                        estimators.push(CoverageEstimator::new_estimator_mean(
                            min_fraction_covered,
                            contig_end_exclusion,
                            false,
                        ));
                        // TODO: Parameterise exclude_mismatches
                    }
                    "count" => {
                        estimators.push(CoverageEstimator::new_estimator_read_count());
                    }
                    "reads_per_base" => {
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
            } else if columns_to_normalise.is_empty()
                && rpkm_column.is_none()
                && tpm_column.is_none()
                && output_format == "sparse"
            {
                debug!("Streaming regular coverage output");
                taker =
                    CoverageTakerType::new_single_float_coverage_streaming_coverage_printer(stream);
                printer = CoveragePrinter::StreamedCoveragePrinter;
            } else {
                debug!(
                    "Cached regular coverage taker with columns to normlise: {:?} and rpkm_column: {:?} and tpm_column: {:?}",
                    columns_to_normalise, rpkm_column, tpm_column
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

        EstimatorsAndTaker {
            estimators,
            taker,
            columns_to_normalise,
            rpkm_column,
            tpm_column,
            printer,
        }
    }

    pub fn print_headers(mut self, entry_type: &str, print_stream: OutputWriter) -> Self {
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
            .print_headers(entry_type, headers, print_stream);
        self
    }
}

fn parse_separator(m: &clap::ArgMatches) -> Option<u8> {
    let single_genome = m.get_flag("single-genome");
    if single_genome {
        Some("0".as_bytes()[0])
    } else if m.contains_id("separator") {
        m.get_one::<char>("separator").map(|c| *c as u8)
    } else if m.contains_id("bam-files") || m.contains_id("reference") {
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
    R: coverm::bam_generator::NamedBamReader,
    T: coverm::bam_generator::NamedBamReaderGenerator<R>,
>(
    bam_generators: Vec<T>,
    m: &clap::ArgMatches,
    estimators_and_taker: &mut EstimatorsAndTaker,
    separator: Option<u8>,
    genomes_and_contigs_option: &Option<GenomesAndContigs>,
    print_stream: &mut OutputWriter,
) {
    let print_zeros = !m.get_flag("no-zeros");
    let flag_filter = FilterParameters::generate_from_clap(m).flag_filters;
    let single_genome = m.get_flag("single-genome");
    let threads = *m.get_one::<u16>("threads").unwrap();
    let reads_mapped = match separator.is_some() || single_genome {
        true => coverm::genome::mosdepth_genome_coverage(
            bam_generators,
            separator.unwrap(),
            &mut estimators_and_taker.taker,
            print_zeros,
            &mut estimators_and_taker.estimators,
            &flag_filter,
            single_genome,
            threads,
        ),

        false => match genomes_and_contigs_option {
            Some(gc) => coverm::genome::mosdepth_genome_coverage_with_contig_names(
                bam_generators,
                gc,
                &mut estimators_and_taker.taker,
                print_zeros,
                &flag_filter,
                &mut estimators_and_taker.estimators,
                threads,
            ),
            None => unreachable!(),
        },
    };

    debug!("Finalising printing ..");
    estimators_and_taker.printer.finalise_printing(
        &estimators_and_taker.taker,
        print_stream,
        Some(&reads_mapped),
        &estimators_and_taker.columns_to_normalise,
        estimators_and_taker.rpkm_column,
        estimators_and_taker.tpm_column,
    );
}

fn doing_metabat(m: &clap::ArgMatches) -> bool {
    let methods: Vec<&str> = m
        .get_many::<String>("methods")
        .unwrap()
        .map(|x| &**x)
        .collect();
    if methods.contains(&"metabat") {
        if methods.len() > 1 {
            error!("Cannot specify the metabat method with any other coverage methods");
            process::exit(1);
        } else {
            return true;
        }
    }
    false
}

#[derive(Debug)]
struct FilterParameters {
    flag_filters: FlagFilter,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    min_mapq: u8, // 255 indicates no filtering
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
}
impl FilterParameters {
    pub fn generate_from_clap(m: &clap::ArgMatches) -> FilterParameters {
        let f = FilterParameters {
            flag_filters: FlagFilter {
                include_improper_pairs: !m.get_flag("proper-pairs-only"),
                include_secondary: m.get_flag("include-secondary"),
                include_supplementary: !m.get_flag("exclude-supplementary"),
            },
            min_aligned_length_single: *m.get_one::<u32>("min-read-aligned-length").unwrap_or(&0),
            min_percent_identity_single: parse_percentage(m, "min-read-percent-identity"),
            min_aligned_percent_single: parse_percentage(m, "min-read-aligned-percent"),
            min_mapq: *m.get_one::<u8>("min-mapq").unwrap_or(&255),
            min_aligned_length_pair: *m
                .get_one::<u32>("min-read-aligned-length-pair")
                .unwrap_or(&0),
            min_percent_identity_pair: parse_percentage(m, "min-read-percent-identity-pair"),
            min_aligned_percent_pair: parse_percentage(m, "min-read-aligned-percent-pair"),
        };
        debug!("Filter parameters set as {:?}", f);
        f
    }

    pub fn add_metabat_filtering_if_required(&mut self, m: &clap::ArgMatches) {
        if doing_metabat(m) {
            info!(
                "Setting single read percent identity threshold at 0.97 for \
                 MetaBAT adjusted coverage, and not filtering out supplementary, \
                 secondary and improper pair alignments"
            );
            // we use >= where metabat uses >. Gah.
            self.min_percent_identity_single = 0.97001;
            self.flag_filters.include_improper_pairs = true;
            self.flag_filters.include_supplementary = true;
            self.flag_filters.include_secondary = true;
        }
    }

    pub fn doing_filtering(&self) -> bool {
        self.min_percent_identity_single > 0.0
            || self.min_percent_identity_pair > 0.0
            || self.min_aligned_percent_single > 0.0
            || self.min_mapq < 255
            || self.min_aligned_percent_pair > 0.0
            || self.min_aligned_length_single > 0
            || self.min_aligned_length_pair > 0
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
    if m.contains_id("bam-file-cache-directory") {
        setup_bam_cache_directory(m.get_one::<String>("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.get_flag("discard-unmapped");
    let sort_threads = *m.get_one::<u16>("threads").unwrap();
    let params = MappingParameters::generate_from_clap(m, mapping_program, reference_tempfile);
    let mut bam_readers = vec![];
    let mut concatenated_reference_name: Option<String> = None;
    let mut concatenated_read_names: Option<String> = None;

    for reference_wise_params in params {
        let index = setup_mapping_index(&reference_wise_params, m, mapping_program);

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
            match m.contains_id("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.get_one::<String>("bam-file-cache-directory").unwrap(),
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
                    index.as_ref(),
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
        sort_threads,
        genome_exclusion,
    };
    vec![gen]
}

fn get_streamed_bam_readers(
    m: &clap::ArgMatches,
    mapping_program: MappingProgram,
    reference_tempfile: &Option<NamedTempFile>,
) -> Vec<BamGeneratorSet<StreamingNamedBamReaderGenerator>> {
    // Check the output BAM directory actually exists and is writeable
    if m.contains_id("bam-file-cache-directory") {
        setup_bam_cache_directory(m.get_one::<String>("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.get_flag("discard-unmapped");

    let params = MappingParameters::generate_from_clap(m, mapping_program, reference_tempfile);
    let mut generator_set = vec![];
    for reference_wise_params in params {
        let mut bam_readers = vec![];
        let index = setup_mapping_index(&reference_wise_params, m, mapping_program);

        let reference = reference_wise_params.reference;
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.contains_id("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.get_one::<String>("bam-file-cache-directory").unwrap(),
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
                    index.as_ref(),
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
            index,
        };
        generator_set.push(to_return);
    }
    generator_set
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
        + std::path::Path::new(reference)
            .file_name()
            .expect("Unable to convert reference to file name")
            .to_str()
            .expect("Unable to covert file name into str")
        + "."
        + std::path::Path::new(read1_path)
            .file_name()
            .expect("Unable to convert read1 name to file name")
            .to_str()
            .expect("Unable to covert file name into str")
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
                    .unwrap_or_else(|_| {
                        panic!(
                            "Unable to canonicalize parent of cache directory {}",
                            cache_directory
                        )
                    })
                    .is_dir()
                {
                    if parent2
                        .metadata()
                        .unwrap_or_else(|_| {
                            panic!(
                                "Unable to get metadata for parent of cache directory {}",
                                cache_directory
                            )
                        })
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
    if m.contains_id("bam-file-cache-directory") {
        setup_bam_cache_directory(m.get_one::<String>("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.get_flag("discard-unmapped");

    let params = MappingParameters::generate_from_clap(m, mapping_program, reference_tempfile);
    let mut generator_set = vec![];
    for reference_wise_params in params {
        let mut bam_readers = vec![];
        let index = setup_mapping_index(&reference_wise_params, m, mapping_program);

        let reference = reference_wise_params.reference;
        let bam_file_cache = |naming_readset| -> Option<String> {
            let bam_file_cache_path;
            match m.contains_id("bam-file-cache-directory") {
                false => None,
                true => {
                    bam_file_cache_path = generate_cached_bam_file_name(
                        m.get_one::<String>("bam-file-cache-directory").unwrap(),
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
                    index.as_ref(),
                    p.read1,
                    p.read2,
                    p.read_format.clone(),
                    p.threads,
                    bam_file_cache(p.read1).as_ref().map(String::as_ref),
                    filter_params.flag_filters.clone(),
                    filter_params.min_aligned_length_single,
                    filter_params.min_percent_identity_single,
                    filter_params.min_aligned_percent_single,
                    filter_params.min_mapq,
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
            index,
        };
        generator_set.push(to_return);
    }
    generator_set
}

fn run_contig<
    R: coverm::bam_generator::NamedBamReader,
    T: coverm::bam_generator::NamedBamReaderGenerator<R>,
>(
    estimators_and_taker: &mut EstimatorsAndTaker,
    bam_readers: Vec<T>,
    print_zeros: bool,
    flag_filters: FlagFilter,
    threads: u16,
    print_stream: &mut OutputWriter,
) {
    let reads_mapped = coverm::contig::contig_coverage(
        bam_readers,
        &mut estimators_and_taker.taker,
        &mut estimators_and_taker.estimators,
        print_zeros,
        &flag_filters,
        threads,
    );

    debug!("Finalising printing ..");

    estimators_and_taker.printer.finalise_printing(
        &estimators_and_taker.taker,
        print_stream,
        Some(&reads_mapped),
        &estimators_and_taker.columns_to_normalise,
        estimators_and_taker.rpkm_column,
        estimators_and_taker.tpm_column,
    );
}

fn set_log_level(matches: &clap::ArgMatches, is_last: bool) {
    set_log_level_bird_tool_utils(matches, is_last, "CoverM", crate_version!());
}
