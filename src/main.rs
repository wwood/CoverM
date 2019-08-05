extern crate coverm;
use coverm::mosdepth_genome_coverage_estimators::*;
use coverm::bam_generator::*;
use coverm::filter;
use coverm::external_command_checker;
use coverm::coverage_takers::*;
use coverm::mapping_parameters::*;
use coverm::coverage_printer::*;
use coverm::shard_bam_reader::*;
use coverm::FlagFilter;
use coverm::CONCATENATED_FASTA_FILE_SEPARATOR;
use coverm::genomes_and_contigs::GenomesAndContigs;
use coverm::genome_exclusion::*;

extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;

use std::env;
use std::str;
use std::process;
use std::collections::HashSet;
use std::path::Path;

extern crate clap;
use clap::*;

#[macro_use]
extern crate log;
extern crate env_logger;
use log::LevelFilter;
use env_logger::Builder;

extern crate tempfile;
use tempfile::NamedTempFile;
#[macro_use]
extern crate lazy_static;

const CONCATENATED_REFERENCE_CACHE_STEM: &str = "coverm-genome";

fn filter_full_help() -> &'static str {
    "coverm filter: Remove alignments with insufficient identity.

Only primary, non-supplementary alignments are considered, and output files
are grouped by reference, but not sorted by position.

Files (both required):
   -b, --bam-files <PATH> ..             Path to reference-sorted BAM file(s)
   -o, --output-bam-files <PATH> ..      Path to corresponding output file(s)

Thresholds:
   --min-read-aligned-length <INT>            Exclude reads with smaller numbers of
                                         aligned bases [default: 0]
   --min-read-percent-identity <FLOAT>        Exclude reads by overall percent
                                         identity e.g. 0.95 for 95%. [default 0.0]
   --min-read-aligned-percent <FLOAT>         Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned. [default 0.0]
   --min-read-aligned-length-pair <INT>       Exclude pairs with smaller numbers of
                                         aligned bases.
                                         Implies --proper-pairs-only. [default: 0]
   --min-read-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95%.
                                         Implies --proper-pairs-only. [default 0.0]
   --min-read-aligned-percent-pair <FLOAT>    Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned.
                                         Implies --proper-pairs-only. [default 0.0]
   --proper-pairs-only                   Require reads to be mapped as proper pairs

Other:
   -t, --threads <INT>                   Number of threads for output compression
                                         [default 1]
   --inverse                             Only keep reads which are unmapped or
                                         align below thresholds. Note that output
                                         records may still be marked as mapped
                                         if they do not meet the thresholds.
                                         [default false]
   --verbose                             Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Example usage:

  coverm filter -b in.bam -o out.bam --min-read-aligned-length 75

Ben J. Woodcroft <benjwoodcroft near gmail.com>"
}

fn contig_full_help() -> &'static str {
    "coverm contig: Calculate read coverage per-contig

Define mapping(s) (required):
  Either define BAM:
   -b, --bam-files <PATH> ..             Path to BAM file(s). These must be
                                         reference sorted (e.g. with samtools sort)
                                         unless --sharded is specified, in which
                                         case they must be read name sorted (e.g.
                                         with samtools sort -n).

  Or do mapping:
   -r, --reference <PATH> ..             FASTA file of contigs or BWA index stem
                                         e.g. concatenated genomes or assembly.
                                         If multiple references FASTA files are
                                         provided and --sharded is specified,
                                         then reads will be mapped to references
                                         separately as sharded BAMs
   -t, --threads <INT>                   Number of threads for mapping / sorting
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

Sharding i.e. multiple reference sets (optional):
   --sharded (experimental)              If -b/--bam-files was used:
                                           Input BAM files are read-sorted alignments
                                           of a set of reads mapped to multiple
                                           reference contig sets. Choose the best
                                           hit for each read pair.

                                         Otherwise if mapping was carried out:
                                           Map reads to each reference, choosing the
                                           best hit for each pair.

Alignment filtering (optional):
   --min-read-aligned-length <INT>            Exclude reads with smaller numbers of
                                         aligned bases [default: 0]
   --min-read-percent-identity <FLOAT>        Exclude reads by overall percent
                                         identity e.g. 0.95 for 95%. [default 0.0]
   --min-read-aligned-percent <FLOAT>         Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned. [default 0.0]
   --min-read-aligned-length-pair <INT>       Exclude pairs with smaller numbers of
                                         aligned bases.
                                         Implies --proper-pairs-only. [default: 0]
   --min-read-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95%.
                                         Implies --proper-pairs-only. [default 0.0]
   --min-read-aligned-percent-pair <FLOAT>    Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned.
                                         Implies --proper-pairs-only. [default 0.0]
   --proper-pairs-only                   Require reads to be mapped as proper pairs

Other arguments (optional):
   -m, --methods <METHOD> [METHOD ..]    Method(s) for calculating coverage.
                                         One or more (space separated) of:
                                           mean (default)
                                           trimmed_mean
                                           coverage_histogram
                                           covered_fraction
                                           covered_bases
                                           variance
                                           length
                                           count
                                           metabat (\"MetaBAT adjusted coverage\")
                                           reads_per_base
                                         A more thorough description of the different
                                         methods is available at
                                         https://github.com/wwood/CoverM

   --output-format FORMAT                Shape of output: 'sparse' for long format,
                                         'dense' for species-by-site.
                                         [default: dense]
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
   --bam-file-cache-directory            Output BAM files generated during
                                         alignment to this directory
   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Ben J. Woodcroft <benjwoodcroft near gmail.com>"
}

fn genome_full_help() -> &'static str {
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
   --genome-definition <FILE>            File containing list of
                                         genome_name<tab>contig
                                         lines to define the genome of each contig
   --single-genome                       All contigs are from the same genome

Define mapping(s) (required):
  Either define BAM:
   -b, --bam-files <PATH> ..             Path to BAM file(s). These must be
                                         reference sorted (e.g. with samtools sort)
                                         unless --sharded is specified, in which
                                         case they must be read name sorted (e.g.
                                         with samtools sort -n).

  Or do mapping:
   -r, --reference <PATH> ..             FASTA file of contigs or BWA index stem
                                         e.g. concatenated genomes or assembly.
                                         If multiple references FASTA files are
                                         provided and --sharded is specified,
                                         then reads will be mapped to references
                                         separately as sharded BAMs
   -t, --threads <INT>                   Number of threads for mapping / sorting
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

Sharding i.e. multiple reference sets (optional):
   --sharded (experimental)              If -b/--bam-files was used:
                                           Input BAM files are read-sorted alignments
                                           of a set of reads mapped to multiple
                                           reference contig sets. Choose the best
                                           hit for each read pair.

                                         Otherwise if mapping was carried out:
                                           Map reads to each reference, choosing the
                                           best hit for each pair.
   --exclude-genomes-from-deshard <FILE> Ignore genomes whose name appears in this
                                         newline-separated file when combining shards.


Alignment filtering (optional):
   --min-read-aligned-length <INT>            Exclude reads with smaller numbers of
                                         aligned bases [default: 0]
   --min-read-percent-identity <FLOAT>        Exclude reads by overall percent
                                         identity e.g. 0.95 for 95%. [default 0.0]
   --min-read-aligned-percent <FLOAT>         Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned. [default 0.0]
   --min-read-aligned-length-pair <INT>       Exclude pairs with smaller numbers of
                                         aligned bases.
                                         Implies --proper-pairs-only. [default: 0]
   --min-read-percent-identity-pair <FLOAT>   Exclude pairs by overall percent
                                         identity e.g. 0.95 for 95%.
                                         Implies --proper-pairs-only. [default 0.0]
   --min-read-aligned-percent-pair <FLOAT>    Exclude reads by percent aligned
                                         bases e.g. 0.95 means 95% of the read's
                                         bases must be aligned.
                                         Implies --proper-pairs-only. [default 0.0]
   --proper-pairs-only                   Require reads to be mapped as proper pairs

Other arguments (optional):
   -m, --methods <METHOD> [METHOD ..]    Method(s) for calculating coverage.
                                         One or more (space separated) of:
                                              relative_abundance (default)
                                              mean
                                              trimmed_mean
                                              coverage_histogram
                                              covered_fraction
                                              covered_bases
                                              variance
                                              length
                                              count
                                              reads_per_base
                                         A more thorough description of the different
                                         methods is available at
                                         https://github.com/wwood/CoverM

   --output-format FORMAT                Shape of output: 'sparse' for long format,
                                         'dense' for species-by-site.
                                         [default: dense]
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
   --bam-file-cache-directory            Output BAM files generated during
                                         alignment to this directory
   --discard-unmapped                    Exclude unmapped reads from cached BAM files.
   -v, --verbose                         Print extra debugging information
   -q, --quiet                           Unless there is an error, do not print
                                         log messages

Ben J. Woodcroft <benjwoodcroft near gmail.com>
"}

fn main(){
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

            let mut genome_names_content: Vec<u8>;

            let separator = parse_separator(m);
            let single_genome = m.is_present("single-genome");
            let genomes_and_contigs_option = match separator.is_some() || single_genome {
                true => None,
                false => {
                    match m.is_present("genome-definition") {
                        true => {
                            Some(coverm::read_genome_definition_file(
                                m.value_of("genome-definition").unwrap()))
                        },
                        false => {
                            let genome_fasta_files: Vec<String> = parse_list_of_genome_fasta_files(m);
                            info!("Reading contig names for {} genomes ..", genome_fasta_files.len());
                            Some(coverm::read_genome_fasta_files(
                                &genome_fasta_files.iter().map(|s| s.as_str()).collect()))
                        }
                    }
                }
            };

            if m.values_of("methods").unwrap().find(|&m| m=="kmer").is_some() {
                if m.values_of("methods").unwrap().len() > 1 {
                    error!("The kmer method cannot be specified with other coverage calculation methods.");
                    process::exit(1);
                }

                match genomes_and_contigs_option {
                    None => {
                        // TODO: Implement.
                        error!("The kmer method must be used with multiple input \
                                genomes and not a separator (for now)");
                        process::exit(1);
                    },
                    Some(genomes_and_contigs) => {
                        // TODO: Much of the following stanza code copied from
                        // contig mode. When ready, refactor.
                        let read_files = m.values_of("single");
                        match read_files {
                            Some(singles) => {
                                if singles.len() > 1 {
                                    error!("Using >1 single read files to map is not yet supported");
                                    process::exit(1);
                                }

                                let num_threads = value_t!(m.value_of("threads"), usize).unwrap();
                                let reference = m.value_of("reference").unwrap();
                                let potential_index_file = format!("{}.covermdb", reference);
                                let index = match Path::new(&potential_index_file).exists() {
                                    true => {
                                        info!("Using pre-existing index {}", potential_index_file);
                                        coverm::kmer_coverage::restore_index::<coverm::pseudoaligner::config::KmerType>(
                                            &potential_index_file)
                                    },
                                    false => {
                                        info!("No pre-existing index file found, generating one ..");
                                        coverm::kmer_coverage::generate_debruijn_index(
                                            reference,
                                            num_threads)
                                    }
                                };

                                coverm::genome_pseudoaligner::calculate_genome_kmer_coverage(
                                    singles.collect::<Vec<_>>()[0],
                                    None,
                                    num_threads,
                                    !m.is_present("no-zeros"),
                                    &index,
                                    &genomes_and_contigs,
                                );
                            },
                            None => {
                                error!("Using -m kmer you must (for now) use --single to input reads");
                                process::exit(1);
                            }
                        }
                    }
                }

            } else {

                let mut estimators_and_taker = EstimatorsAndTaker::generate_from_clap(
                    m, print_stream);
                estimators_and_taker = estimators_and_taker.print_headers(
                    &"Genome", &mut std::io::stdout());
                let filter_params = FilterParameters::generate_from_clap(m);

                // This would be better as a separate functi1on to make this function
                // smaller, but I find this hard because functions cannot return a
                // trait.
                let mut genome_exclusion_filter_separator_type: Option<SeparatorGenomeExclusionFilter> = None;
                let mut genome_exclusion_filter_non_type: Option<NoExclusionGenomeFilter> = None;
                let mut genome_exclusion_genomes_and_contigs: Option<GenomesAndContigsExclusionFilter> = None;
                enum GenomeExclusionTypes {
                    SeparatorType, NoneType, GenomesAndContigsType
                }
                let genome_exclusion_type = {
                    if m.is_present("sharded") {
                        if m.is_present("exclude-genomes-from-deshard") {
                            let filename = m.value_of("exclude-genomes-from-deshard").unwrap();
                            genome_names_content = std::fs::read(filename)
                                .expect(&format!(
                                    "Failed to open file '{}' containing list of excluded genomes",
                                    filename));
                            let mut genome_names_hash: HashSet<&[u8]> = HashSet::new();
                            for n in genome_names_content.split(|s| *s == b"\n"[0]) {
                                if n != b"" {
                                    genome_names_hash.insert(n);
                                }
                            }
                            if genome_names_hash.is_empty() {
                                warn!("No genomes read in that are to be excluded from desharding process");
                                genome_exclusion_filter_non_type = Some(NoExclusionGenomeFilter{});
                                GenomeExclusionTypes::NoneType
                            } else {
                                info!("Read in {} distinct genomes to exclude from desharding process e.g. '{}'",
                                      genome_names_hash.len(),
                                      std::str::from_utf8(genome_names_hash.iter().next().unwrap())
                                      .unwrap());
                                if separator.is_some() {
                                    genome_exclusion_filter_separator_type = Some(
                                        SeparatorGenomeExclusionFilter {
                                            split_char: separator.unwrap(),
                                            excluded_genomes: genome_names_hash
                                        });
                                    GenomeExclusionTypes::SeparatorType
                                } else {
                                    match genomes_and_contigs_option {
                                        Some(ref gc) => {
                                            genome_exclusion_genomes_and_contigs = Some(
                                                GenomesAndContigsExclusionFilter {
                                                    genomes_and_contigs: gc,
                                                    excluded_genomes: genome_names_hash,
                                                });
                                            GenomeExclusionTypes::GenomesAndContigsType
                                        },
                                        None => unreachable!()
                                    }
                                }
                            }
                        } else {
                            debug!("Not excluding any genomes during the deshard process");
                            genome_exclusion_filter_non_type = Some(NoExclusionGenomeFilter{});
                            GenomeExclusionTypes::NoneType
                        }
                    } else {
                        genome_exclusion_filter_non_type = Some(NoExclusionGenomeFilter{});
                        GenomeExclusionTypes::NoneType
                    }
                };

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
                            &mut estimators_and_taker,
                            separator,
                            &genomes_and_contigs_option);

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
                            },
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
                            },
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
                            },
                        }
                    } else {
                        run_genome(
                            coverm::bam_generator::generate_named_bam_readers_from_bam_files(
                                bam_files),
                            m,
                            &mut estimators_and_taker,
                            separator,
                            &genomes_and_contigs_option);
                    }
                } else {
                    external_command_checker::check_for_bwa();
                    external_command_checker::check_for_samtools();

                    // Generate a temporary file of concatenated genomes if needed.
                    let mut concatenated_genomes: Option<NamedTempFile> = None;
                    if !m.is_present("reference") && !m.is_present("bam-files") {
                        let list_of_genome_fasta_files = parse_list_of_genome_fasta_files(m);
                        info!("Generating concatenated reference FASTA file of {} genomes ..",
                              list_of_genome_fasta_files.len());
                        concatenated_genomes = Some(
                            coverm::bwa_index_maintenance::generate_concatenated_fasta_file(
                                &list_of_genome_fasta_files));
                    }

                    if filter_params.doing_filtering() {
                        debug!("Mapping and filtering..");
                        let generator_sets = get_streamed_filtered_bam_readers(
                            m, &concatenated_genomes, &filter_params);
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
                            &mut estimators_and_taker,
                            separator,
                            &genomes_and_contigs_option);
                    } else if m.is_present("sharded") {
                        match genome_exclusion_type {
                            GenomeExclusionTypes::NoneType => {
                                run_genome(
                                    get_sharded_bam_readers(
                                        m,
                                        &concatenated_genomes,
                                        &genome_exclusion_filter_non_type.unwrap()),
                                    m,
                                    &mut estimators_and_taker,
                                    separator,
                                    &genomes_and_contigs_option);
                            },
                            GenomeExclusionTypes::SeparatorType => {
                                run_genome(
                                    get_sharded_bam_readers(
                                        m,
                                        &concatenated_genomes,
                                        &genome_exclusion_filter_separator_type.unwrap()),
                                    m,
                                    &mut estimators_and_taker,
                                    separator,
                                    &genomes_and_contigs_option);
                            },
                            GenomeExclusionTypes::GenomesAndContigsType => {
                                run_genome(
                                    get_sharded_bam_readers(
                                        m,
                                        &concatenated_genomes,
                                        &genome_exclusion_genomes_and_contigs.unwrap()),
                                    m,
                                    &mut estimators_and_taker,
                                    separator,
                                    &genomes_and_contigs_option);
                            },
                        }
                    } else {
                        let generator_sets = get_streamed_bam_readers(m, &concatenated_genomes);
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
                            &mut estimators_and_taker,
                            separator,
                            &genomes_and_contigs_option);
                    };
                }
            }
        },
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
                let reader = bam::Reader::from_path(bam).expect(
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
                    filter_params.min_aligned_percent_pair,
                    !m.is_present("inverse"));

                let mut record = bam::record::Record::new();
                while filtered.read(&mut record).is_ok() {
                    debug!("Writing.. {:?}", record.qname());
                    writer.write(&record).expect("Failed to write BAM record");
                }
            }
        },
        Some("contig") => {
            let m = matches.subcommand_matches("contig").unwrap();
            if m.is_present("full-help") {
                println!("{}", contig_full_help());
                process::exit(1);
            }
            set_log_level(m, true);
            let print_zeros = !m.is_present("no-zeros");
            let filter_params = FilterParameters::generate_from_clap(m);

            let methods: Vec<&str> = m.values_of("methods").unwrap().collect();
            if methods.contains(&"kmer") {
                if methods.len() > 1 {
                    error!("The kmer method cannot be specified with other coverage calculation methods.");
                    process::exit(1);
                }
                let mapping_parameters = MappingParameters::generate_from_clap(&m, &None);

                let num_threads = value_t!(m.value_of("threads"), usize).unwrap();
                let reference = m.value_of("reference").unwrap();
                let potential_index_file = format!("{}.covermdb", reference);
                let index = match Path::new(&potential_index_file).exists() {
                    true => {
                        info!("Using pre-existing index {}", potential_index_file);
                        coverm::kmer_coverage::restore_index::<coverm::pseudoaligner::config::KmerType>(
                            &potential_index_file)
                    },
                    false => {
                        info!("No pre-existing index file found, generating one ..");
                        coverm::kmer_coverage::generate_debruijn_index(
                            reference,
                            num_threads)
                    }
                };

                let mapping_pairs = mapping_parameters.read1.len() > 0;
                coverm::kmer_coverage::calculate_contig_kmer_coverage(
                    if mapping_pairs { mapping_parameters.read1[0] } else { mapping_parameters.unpaired[0] },
                    if mapping_pairs { Some(mapping_parameters.read2[0]) } else { None },
                    num_threads,
                    !m.is_present("no-zeros"),
                    index,
                );
            }

            else {
                let mut estimators_and_taker = EstimatorsAndTaker::generate_from_clap(
                    m, &mut print_stream);
                estimators_and_taker = estimators_and_taker.print_headers(
                    &"Contig", &mut std::io::stdout());

                if m.is_present("bam-files") {
                    let bam_files: Vec<&str> = m.values_of("bam-files").unwrap().collect();
                    if filter_params.doing_filtering() {
                        let bam_readers = coverm::bam_generator::generate_filtered_bam_readers_from_bam_files(
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
                    } else if m.is_present("sharded") {
                        external_command_checker::check_for_samtools();
                        let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
                        let bam_readers = coverm::shard_bam_reader::generate_sharded_bam_reader_from_bam_files(
                            bam_files, sort_threads, &NoExclusionGenomeFilter{});
                        run_contig(
                            &mut estimators_and_taker,
                            bam_readers,
                            print_zeros,
                            filter_params.flag_filters);
                    } else {
                        let bam_readers = coverm::bam_generator::generate_named_bam_readers_from_bam_files(
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
                        let generator_sets = get_streamed_filtered_bam_readers(
                            m, &None, &filter_params);
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
                    } else if m.is_present("sharded") {
                        let generator_sets = get_sharded_bam_readers(
                            m, &None, &NoExclusionGenomeFilter{});
                        run_contig(
                            &mut estimators_and_taker,
                            generator_sets,
                            print_zeros,
                            filter_params.flag_filters);
                    } else {
                        debug!("Not filtering..");
                        let generator_sets = get_streamed_bam_readers(m, &None);
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
            }
        },
        Some("make") => {
            let m = matches.subcommand_matches("make").unwrap();
            set_log_level(m, true);
            external_command_checker::check_for_bwa();
            external_command_checker::check_for_samtools();

            let output_directory = m.value_of("output-directory").unwrap();
            setup_bam_cache_directory(output_directory);
            let params = MappingParameters::generate_from_clap(&m, &None);
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
        Some("debruijn_index") => {
            let m = matches.subcommand_matches("debruijn_index").unwrap();
            set_log_level(m, true);

            let reference = m.value_of("reference").unwrap();
            let output = format!("{}.covermdb", reference);
            let num_threads = value_t!(m.value_of("threads"), usize).unwrap();

            info!("Generating index ..");
            let index = coverm::kmer_coverage::generate_debruijn_index::<coverm::pseudoaligner::config::KmerType>(
                reference, num_threads);

            info!("Saving index ..");
            coverm::kmer_coverage::save_index(
                index, &output);
            info!("Saving complete");
        }
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

fn parse_percentage(m: &clap::ArgMatches, parameter: &str) -> f32 {
    match m.is_present(parameter) {
        true => {
            let mut percentage = value_t!(
                m.value_of(parameter), f32).unwrap();
            if percentage >= 1.0 && percentage <= 100.0 {
                percentage = percentage / 100.0;
            } else if percentage < 0.0 || percentage > 100.0 {
                error!("Invalid alignment percentage: '{}'", percentage);
                process::exit(1);
            }
            info!("Using {} {}%", parameter, percentage*100.0);
            percentage
        },
        false => 0.0
    }
}

impl<'a> EstimatorsAndTaker<'a> {
    pub fn generate_from_clap(
        m: &clap::ArgMatches, stream: &'a mut std::io::Stdout)
        -> EstimatorsAndTaker<'a> {
        let mut estimators = vec!();
        let min_fraction_covered = parse_percentage(&m, "min-covered-fraction");
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
                            error!("error: Trim bounds must be between 0 and 1, and \
                                    min must be less than max, found {} and {}", min, max);
                            process::exit(1);
                        }
                        estimators.push(CoverageEstimator::new_estimator_trimmed_mean(
                            min, max, min_fraction_covered, contig_end_exclusion));
                    },
                    &"covered_fraction" => {
                        estimators.push(CoverageEstimator::new_estimator_covered_fraction(
                            min_fraction_covered, contig_end_exclusion));
                    },
                    &"covered_bases" => {
                        estimators.push(CoverageEstimator::new_estimator_covered_bases(
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
                    &"reads_per_base" => {
                        estimators.push(CoverageEstimator::new_estimator_reads_per_base());
                    }
                    _ => unreachable!()
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
                    _ => unreachable!()
                }
            }

        }

        // Check that min-covered-fraction is being used as expected
        if min_fraction_covered != 0.0 {
            let die = |estimator_name| {
                error!("The '{}' coverage estimator cannot be used when \
                        --min-covered-fraction is > 0 as it does not calculate \
                        the covered fraction. You may wish to set the \
                        --min-covered-fraction to 0 and/or run this estimator \
                        separately.", estimator_name);
                process::exit(1)
            };
            for e in &estimators {
                match e {
                    CoverageEstimator::ReadCountCalculator{..} => { die("counts") },
                    CoverageEstimator::ReferenceLengthCalculator{..} => { die("length") },
                    CoverageEstimator::ReadsPerBaseCalculator{..} => { die("reads_per_base") }
                    _ => {}
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
                str::from_utf8(separator_str).unwrap());
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


fn parse_list_of_genome_fasta_files(m: &clap::ArgMatches) -> Vec<String> {
    match m.is_present("genome-fasta-files") {
        true => {
            m.values_of("genome-fasta-files").unwrap().map(|s| s.to_string()).collect()
        },
        false => {
            if m.is_present("genome-fasta-directory") {
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
                if genome_fasta_files.len() == 0 {
                    error!("Found 0 genomes from the genome-fasta-directory, cannot continue.");
                    process::exit(1);
                }
                genome_fasta_files // return
            } else {
                error!("Either a separator (-s) or path(s) to genome FASTA files \
                        (with -d or -f) must be given");
                process::exit(1);
            }
        }
    }
}


fn run_genome<'a,
              R: coverm::bam_generator::NamedBamReader,
              T: coverm::bam_generator::NamedBamReaderGenerator<R>> (
    bam_generators: Vec<T>,
    m: &clap::ArgMatches,
    estimators_and_taker: &'a mut EstimatorsAndTaker<'a>,
    separator: Option<u8>,
    genomes_and_contigs_option: &Option<GenomesAndContigs>) {

    let print_zeros = !m.is_present("no-zeros");
    let proper_pairs_only = m.is_present("proper-pairs-only");
    let single_genome = m.is_present("single-genome");
    let reads_mapped = match separator.is_some() || single_genome {
        true => {
            coverm::genome::mosdepth_genome_coverage(
                bam_generators,
                separator.unwrap(),
	              &mut estimators_and_taker.taker,
                print_zeros,
                &mut estimators_and_taker.estimators,
                proper_pairs_only,
                single_genome)
        },

        false => {
            match genomes_and_contigs_option {
                Some(gc) =>
                    coverm::genome::mosdepth_genome_coverage_with_contig_names(
                        bam_generators,
                        gc,
                        &mut estimators_and_taker.taker,
                        print_zeros,
                        proper_pairs_only,
                        &mut estimators_and_taker.estimators),
                None => unreachable!()
            }
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
                    error!("Cannot specify the metabat method with any other coverage methods");
                    process::exit(1);
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
            min_aligned_length_single: match m.is_present("min-read-aligned-length") {
                true => value_t!(m.value_of("min-read-aligned-length"), u32).unwrap(),
                false => 0
            },
            min_percent_identity_single: parse_percentage(&m, "min-read-percent-identity"),
            min_aligned_percent_single: parse_percentage(&m, "min-read-aligned-percent"),
            min_aligned_length_pair: match m.is_present("min-read-aligned-length-pair") {
                true => value_t!(m.value_of("min-read-aligned-length-pair"), u32).unwrap(),
                false => 0
            },
            min_percent_identity_pair: parse_percentage(&m, "min-read-percent-identity-pair"),
            min_aligned_percent_pair: parse_percentage(&m, "min-read-aligned-percent-pair"),
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

fn get_sharded_bam_readers<'a, 'b, T>(
    m: &'a clap::ArgMatches,
    reference_tempfile: &'a Option<NamedTempFile>,
    genome_exclusion: &'b T)
    -> Vec<ShardedBamReaderGenerator<'b, T>>
where T: GenomeExclusion {

    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");
    let sort_threads = m.value_of("threads").unwrap().parse::<i32>().unwrap();
    let params = MappingParameters::generate_from_clap(&m, &reference_tempfile);
    let mut bam_readers = vec![];
    let mut concatenated_reference_name: Option<String> = None;
    let mut concatenated_read_names: Option<String> = None;

    for reference_wise_params in params {
        let index = coverm::bwa_index_maintenance::generate_bwa_index(
            reference_wise_params.reference);

        let reference = reference_wise_params.reference;
        let reference_name = std::path::Path::new(reference).file_name()
            .expect("Unable to convert reference to file name").to_str()
            .expect("Unable to covert file name into str").to_string();
        concatenated_reference_name = match concatenated_reference_name {
            Some(prev) => Some(format!("{}|{}",prev,reference_name)),
            None => Some(reference_name)
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
                            None => reference
                        },
                        naming_readset);
                    info!("Caching BAM file to {}", bam_file_cache_path);
                    Some(bam_file_cache_path)
                }
            }
        };

        for p in reference_wise_params {
            bam_readers.push(
                coverm::shard_bam_reader::generate_named_sharded_bam_readers_from_reads(
                    index.index_path(),
                    p.read1,
                    p.read2,
                    p.read_format.clone(),
                    p.threads,
                    bam_file_cache(p.read1).as_ref().map(String::as_ref),
                    discard_unmapped,
                    p.bwa_options));
            let name = &std::path::Path::new(p.read1).file_name()
                .expect("Unable to convert read1 name to file name").to_str()
                .expect("Unable to covert file name into str").to_string();
            concatenated_read_names = match concatenated_read_names {
                Some(prev) => Some(format!("{}|{}", prev, name)),
                None => Some(name.to_string())
            };
        }

        debug!("Finished BAM setup");
    };
    let gen = ShardedBamReaderGenerator {
        stoit_name: format!(
            "{}/{}",
            concatenated_reference_name.unwrap(),
            concatenated_read_names.unwrap()),
        read_sorted_bam_readers: bam_readers,
        sort_threads: sort_threads,
        genome_exclusion: genome_exclusion,
    };
    return vec!(gen);
}

fn get_streamed_bam_readers<'a>(
    m: &'a clap::ArgMatches,
    reference_tempfile: &'a Option<NamedTempFile>)
    -> Vec<BamGeneratorSet<StreamingNamedBamReaderGenerator>> {

    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params = MappingParameters::generate_from_clap(&m, &reference_tempfile);
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
                        match reference_tempfile {
                            Some(_) => CONCATENATED_REFERENCE_CACHE_STEM,
                            None => reference
                        },
                        naming_readset);
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
                    p.bwa_options,
                    reference_tempfile.is_none()));
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
    debug!("Constructing BAM file cache name in directory {}, reference {}, read1_path {}",
           directory, reference, read1_path);
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
                error!("Cache directory {} does not appear to be writeable, not continuing",
                       cache_directory);
                process::exit(1);
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
                            error!(
                                "The parent directory of the (currently non-existent) \
                                 cache directory {} is not writeable, not continuing",
                                cache_directory);
                            process::exit(1);
                        } else {
                            info!("Creating cache directory {}", cache_directory);
                            std::fs::create_dir(path).expect("Unable to create cache directory");
                        }
                } else {
                    error!("The parent directory of the cache directory {} does not \
                            yet exist, so not creating that cache directory, and not continuing.",
                           cache_directory);
                    process::exit(1);
                }
            },
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
        error!("Failed to create test file in bam cache directory: {}", tf_result.err().unwrap());
        process::exit(1);
    }
}

fn get_streamed_filtered_bam_readers(
    m: &clap::ArgMatches,
    reference_tempfile: &Option<NamedTempFile>,
    filter_params: &FilterParameters)
    -> Vec<BamGeneratorSet<StreamingFilteredNamedBamReaderGenerator>> {

    // Check the output BAM directory actually exists and is writeable
    if m.is_present("bam-file-cache-directory") {
        setup_bam_cache_directory(m.value_of("bam-file-cache-directory").unwrap());
    }
    let discard_unmapped = m.is_present("discard-unmapped");

    let params = MappingParameters::generate_from_clap(&m, &reference_tempfile);
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
                        match reference_tempfile {
                            Some(_) => CONCATENATED_REFERENCE_CACHE_STEM,
                            None => reference
                        },
                        naming_readset);
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
                    discard_unmapped,
                    reference_tempfile.is_none()));
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

fn build_cli() -> App<'static, 'static> {
    // specify static lazily because need to define it at runtime.
    lazy_static! {
        static ref CONTIG_HELP: String = format!(
            "
                            {}
              {}

{}

  coverm contig --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna

{}

  coverm contig --method metabat --bam-files my.bam
    --bam-file-cache-directory saved_bam_files

See coverm contig --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "coverm contig"),
            ansi_term::Colour::Green.paint(
                "Calculate coverage of individual contigs"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate mean coverage from reads and assembly:"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate MetaBAT adjusted coverage from a sorted BAM file, saving
the unfiltered BAM files in the saved_bam_files folder:")
        ).to_string();

        static ref GENOME_HELP: String = format!(
            "
                            {}
               {}

{}

  coverm genome --coupled read1.fastq.gz read2.fastq.gz
    --reference assembly.fna --separator '~'

{}

  coverm genome --bam-files my.bam --genome-fasta-directory genomes_directory/

See coverm genome --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "coverm genome"),
            ansi_term::Colour::Green.paint(
                "Calculate coverage of individual genomes"),
            ansi_term::Colour::Purple.paint(
                "Example: Map paired reads to a reference where the FASTA header separates
the genome name from the contig name with '~' e.g. >genome10~contig15"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate coverage of genomes defined as .fna files in
genomes_directory/ from a sorted BAM file:"),
        ).to_string();

        static ref FILTER_HELP: String = format!(
            "
                            {}
                     {}

{}

  coverm filter --bam-files input.bam --output-bam filtered.bam
    --min-read-aligned-length 50

{}

  coverm filter -b input.bam -o inverse_filtered.bam --inverse
    --min-read-percent-identity 0.95 --threads 16

See coverm filter --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint(
                "coverm filter"),
            ansi_term::Colour::Green.paint(
                "Filter BAM file alignments"),
            ansi_term::Colour::Purple.paint(
                "Example: Filter a BAM file by removing alignments shorter than 50bp:"),
            ansi_term::Colour::Purple.paint(
                "Example: Filter inverse: Keep alignments that have <95% alignment identity\n\
                 and those which do map at all. Note that the output BAM file will likely\n\
                 records that are still mapped, but align with < 95% identity. Use 16\n\
                 threads for output compression:"),
        ).to_string();
    }

    let make_help: &'static str =
        "coverm make: Generate BAM files through mapping.

Output (required):
   -o, --output-directory <DIR>          Where generated BAM files will go

Mapping parameters:
   -r, --reference <PATH>                FASTA file(s) of contig(s) or BWA index stem.
                                         If multiple reference FASTA files are provided,
                                         reads will be mapped to each reference separately
                                         and create sharded BAM files.
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

Main subcommands:
\tcontig\tCalculate coverage of contigs
\tgenome\tCalculate coverage of genomes

Less used utility subcommands:
\tmake\tGenerate BAM files through alignment
\tfilter\tRemove (or only keep) alignments with insufficient identity

Other options:
\t-V, --version\tPrint version information

Ben J. Woodcroft <benjwoodcroft near gmail.com>
")
        .global_setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(
            SubCommand::with_name("genome")
                .about("Calculate coverage of genomes")
                .help(GENOME_HELP.as_str())

                .arg(Arg::with_name("full-help")
                     .long("full-help"))

                .arg(Arg::with_name("bam-files")
                     .short("b")
                     .long("bam-files")
                     .multiple(true)
                     .takes_value(true))
                .arg(Arg::with_name("sharded")
                     .long("sharded")
                     .required(false))
                .arg(Arg::with_name("exclude-genomes-from-deshard")
                     .long("exclude-genomes-from-deshard")
                     .requires("sharded")
                     .takes_value(true))
                .arg(Arg::with_name("read1")
                     .short("-1")
                     .multiple(true)
                     .takes_value(true)
                     .requires("read2")
                     .required_unless_one(
                         &["bam-files","coupled","interleaved","single","full-help"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("read2")
                     .short("-2")
                     .multiple(true)
                     .takes_value(true)
                     .requires("read1")
                     .required_unless_one(
                         &["bam-files","coupled","interleaved","single","full-help"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("coupled")
                     .short("-c")
                     .long("coupled")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["bam-files","read1","interleaved","single","full-help"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("interleaved")
                     .long("interleaved")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["bam-files","read1","coupled","single","full-help"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("single")
                     .long("single")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["bam-files","read1","coupled","interleaved","full-help"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("reference")
                     .short("-r")
                     .long("reference")
                     .takes_value(true)
                     .multiple(true)
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
                     .requires("reference")) // TODO: Relax this for autoconcatenation
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
                         &["genome-fasta-files","genome-fasta-directory","single-genome",
                           "genome-definition","full-help"])
                     .takes_value(true))
                .arg(Arg::with_name("genome-fasta-files")
                     .short("f")
                     .long("genome-fasta-files")
                     .multiple(true)
                     .conflicts_with("separator")
                     .conflicts_with("genome-fasta-directory")
                     .conflicts_with("single-genome")
                     .required_unless_one(
                         &["separator","genome-fasta-directory","single-genome",
                           "genome-definition","full-help"])
                     .takes_value(true))
                .arg(Arg::with_name("genome-fasta-directory")
                     .short("d")
                     .long("genome-fasta-directory")
                     .conflicts_with("separator")
                     .conflicts_with("genome-fasta-files")
                     .conflicts_with("single-genome")
                     .required_unless_one(
                         &["genome-fasta-files","separator","single-genome",
                           "genome-definition","full-help"])
                     .takes_value(true))
                .arg(Arg::with_name("genome-fasta-extension")
                     .short("x")
                     .long("genome-fasta-extension")
                     // Unsure why, but uncommenting causes test failure - clap
                     // bug?
                     //.requires("genome-fasta-directory")
                     .default_value("fna")
                     .takes_value(true))
                .arg(Arg::with_name("genome-definition")
                     .long("genome-definition")
                     .conflicts_with("separator")
                     .conflicts_with("genome-fasta-files")
                     .conflicts_with("genome-fasta-directory")
                     .conflicts_with("single-genome")
                     .required_unless_one(
                         &["genome-fasta-files","separator","single-genome",
                           "genome-fasta-directory","full-help"])
                     .takes_value(true))
                .arg(Arg::with_name("single-genome")
                     .long("single-genome")
                     .conflicts_with("separator")
                     .conflicts_with("genome-fasta-files")
                     .conflicts_with("genome-fasta-directory")
                     .conflicts_with("genome-definition"))

                .arg(Arg::with_name("min-read-aligned-length")
                     .long("min-read-aligned-length")
                     .takes_value(true))
                .arg(Arg::with_name("min-read-percent-identity")
                     .long("min-read-percent-identity")
                     .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-percent")
                     .long("min-read-aligned-percent")
                     .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-length-pair")
                     .long("min-read-aligned-length-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))
                .arg(Arg::with_name("min-read-percent-identity-pair")
                     .long("min-read-percent-identity-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))
                .arg(Arg::with_name("min-read-aligned-percent-pair")
                     .long("min-read-aligned-percent-pair")
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
                         "covered_bases",
                         "variance",
                         "length",
                         "count",
                         "reads_per_base",
                         "kmer"])
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
                     .default_value("dense"))

                .arg(Arg::with_name("verbose")
                     .short("v")
                     .long("verbose"))
                .arg(Arg::with_name("quiet")
                     .short("q")
                     .long("quiet")))
        .subcommand(
            SubCommand::with_name("contig")
                .about("Calculate coverage of contigs")
                .help(CONTIG_HELP.as_str())

                .arg(Arg::with_name("full-help")
                     .long("full-help"))

                .arg(Arg::with_name("bam-files")
                     .short("b")
                     .long("bam-files")
                     .multiple(true)
                     .takes_value(true))
                .arg(Arg::with_name("sharded")
                    .long("sharded")
                    .required(false))
                .arg(Arg::with_name("read1")
                     .short("-1")
                     .multiple(true)
                     .takes_value(true)
                     .requires("read2")
                     .required_unless_one(
                         &["bam-files","coupled","interleaved","single","full-help"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("read2")
                     .short("-2")
                     .multiple(true)
                     .takes_value(true)
                     .requires("read1")
                     .required_unless_one(
                         &["bam-files","coupled","interleaved","single","full-help"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("coupled")
                     .short("-c")
                     .long("coupled")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["bam-files","read1","interleaved","single","full-help"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("interleaved")
                     .long("interleaved")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["bam-files","read1","coupled","single","full-help"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("single")
                     .long("single")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(
                         &["bam-files","read1","coupled","interleaved","full-help"])
                     .conflicts_with("bam-files"))
                .arg(Arg::with_name("reference")
                     .short("-r")
                     .long("reference")
                     .takes_value(true)
                     .multiple(true)
                     .required_unless_one(&["bam-files","full-help"])
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

                .arg(Arg::with_name("min-read-aligned-length")
                     .long("min-read-aligned-length")
                     .takes_value(true))
                .arg(Arg::with_name("min-read-percent-identity")
                     .long("min-read-percent-identity")
                     .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-percent")
                     .long("min-read-aligned-percent")
                     .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-length-pair")
                     .long("min-read-aligned-length-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))
                .arg(Arg::with_name("min-read-percent-identity-pair")
                     .long("min-read-percent-identity-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))
                .arg(Arg::with_name("min-read-aligned-percent-pair")
                     .long("min-read-aligned-percent-pair")
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
                         "covered_bases",
                         "variance",
                         "length",
                         "count",
                         "metabat",
                         "reads_per_base",
                         "kmer"])
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
                     .default_value("dense"))
                .arg(Arg::with_name("verbose")
                     .short("v")
                     .long("verbose"))
                .arg(Arg::with_name("quiet")
                     .short("q")
                     .long("quiet")))
        .subcommand(
            SubCommand::with_name("filter")
                .about("Remove alignments with insufficient identity")
                .help(FILTER_HELP.as_str())

                .arg(Arg::with_name("full-help")
                     .long("full-help"))

                .arg(Arg::with_name("bam-files")
                     .short("b")
                     .long("bam-files")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(&["full-help"]))
                .arg(Arg::with_name("output-bam-files")
                     .short("o")
                     .long("output-bam-files")
                     .multiple(true)
                     .takes_value(true)
                     .required_unless_one(&["full-help"]))
                .arg(Arg::with_name("inverse")
                     .long("inverse"))

                .arg(Arg::with_name("min-read-aligned-length")
                     .long("min-read-aligned-length")
                     .takes_value(true))
                .arg(Arg::with_name("min-read-percent-identity")
                     .long("min-read-percent-identity")
                     .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-percent")
                     .long("min-read-aligned-percent")
                     .takes_value(true))
                .arg(Arg::with_name("min-read-aligned-length-pair")
                     .long("min-read-aligned-length-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))
                .arg(Arg::with_name("min-read-percent-identity-pair")
                     .long("min-read-percent-identity-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))
                .arg(Arg::with_name("min-read-aligned-percent-pair")
                     .long("min-read-aligned-percent-pair")
                     .takes_value(true)
                     .requires("proper-pairs-only"))

                .arg(Arg::with_name("proper-pairs-only")
                     .long("proper-pairs-only"))
                .arg(Arg::with_name("threads")
                     .long("threads")
                     .short("t")
                     .default_value("1"))

                .arg(Arg::with_name("verbose")
                     // .short("v") // Do not use since could be confused with
                     // inverse (a la grep -v)
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
                     .multiple(true)
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
                     .requires("reference")))

        .subcommand(
            SubCommand::with_name("debruijn_index")
                .about("Generate a DeBruijn index file")
                .arg(Arg::with_name("reference")
                     .short("-r")
                     .long("reference")
                     .takes_value(true)
                     .required(true))
                .arg(Arg::with_name("threads")
                     .short("-t")
                     .long("threads")
                     .default_value("1")
                     .takes_value(true)));
}
