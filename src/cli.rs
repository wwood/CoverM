use bird_tool_utils::clap_utils::{add_clap_verbosity_flags, default_roff, monospace_roff};
use bird_tool_utils_man::prelude::{Author, Example, Flag, Manual, Opt, Section};
use clap::*;
use clap_complete::*;
use galah::cluster_argument_parsing::GalahClustererCommandDefinition;
use roff::bold as roff_bold;
use roff::Roff;

// See https://github.com/rust-cli/roff-rs/issues/19
fn bold(s: &str) -> String {
    Roff::new().text([roff_bold(s)]).to_roff()
}

const MAPPING_SOFTWARE_LIST: &[&str] = &[
    "bwa-mem",
    "bwa-mem2",
    "minimap2-sr",
    "minimap2-ont",
    "minimap2-pb",
    "minimap2-hifi",
    "minimap2-lr-hq",
    "minimap2-no-preset",
    "strobealign",
    "minibwa",
    "rammap-sr",
    "rammap-ont",
    "rammap-pb",
    "rammap-hifi",
    "rammap-lr-hq",
    "rammap-no-preset",
];
const DEFAULT_MAPPING_SOFTWARE: &str = "strobealign";
const DEFAULT_MAKEDB_MAPPING_SOFTWARE: &str = "minimap2-sr";

lazy_static! {
    pub static ref COVERM_CLUSTER_COMMAND_DEFINITION: GalahClustererCommandDefinition = {
        galah::cluster_argument_parsing::GalahClustererCommandDefinition {
            dereplication_ani_argument: "dereplication-ani".to_string(),
            dereplication_prethreshold_ani_argument: "dereplication-prethreshold-ani".to_string(),
            dereplication_quality_formula_argument: "dereplication-quality-formula".to_string(),
            dereplication_run_checkm2_argument: "run-checkm2".to_string(),
            dereplication_checkm2_db_path_argument: "checkm2-db-path".to_string(),
            dereplication_precluster_method_argument: "dereplication-precluster-method".to_string(),
            dereplication_cluster_method_argument: "dereplication-cluster-method".to_string(),
            dereplication_aligned_fraction_argument: "dereplication-aligned-fraction".to_string(),
            dereplication_small_genomes_argument: "dereplication-small-genomes".to_string(),
            dereplication_cluster_contigs_argument: "dereplication-cluster-contigs".to_string(),
            dereplication_small_contigs_argument: "dereplication-small-contigs".to_string(),
            dereplication_large_contigs_argument: "dereplication-large-contigs".to_string(),
            dereplication_fraglen_argument: "dereplication-fragment-length".to_string(),
            dereplication_low_memory_argument: "dereplication-low-memory".to_string(),
            dereplication_reference_genomes_argument: "dereplication-reference-genomes".to_string(),
            dereplication_reference_genomes_list_argument: "dereplication-reference-genomes-list"
                .to_string(),
            dereplication_output_cluster_definition_file: "dereplication-output-cluster-definition"
                .to_string(),
            dereplication_output_representative_fasta_directory:
                "dereplication-output-representative-fasta-directory".to_string(),
            dereplication_output_representative_fasta_directory_copy:
                "dereplication-output-representative-fasta-directory-copy".to_string(),
            dereplication_output_representative_list: "dereplication-output-representative-list"
                .to_string(),
        }
    };
}

fn add_mapping_options(manual: Manual) -> Manual {
    manual.custom(
        Section::new("Mapping algorithm options")
            .option(Opt::new("NAME").short("-p").long("--mapper").help(&format!(
                "Underlying mapping software used {}. One of: {}",
                default_roff("strobealign"),
                bird_tool_utils::clap_utils::table_roff(&[
                    &["name", "description"],
                    &[
                        &monospace_roff("strobealign"),
                        "strobealign using default parameters"
                    ],
                    &[
                        &monospace_roff("bwa-mem"),
                        "bwa mem using default parameters"
                    ],
                    &[
                        &monospace_roff("bwa-mem2"),
                        "bwa-mem2 using default parameters"
                    ],
                    &[
                        &monospace_roff("minimap2-sr"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x sr"))
                    ],
                    &[
                        &monospace_roff("minimap2-lr-hq"),
                        &format!("minimap2 with '{}' option (recommended for Oxford Nanopore R10+ reads)", &monospace_roff("-x lr:hq"))
                    ],
                    &[
                        &monospace_roff("minimap2-ont"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-ont"))
                    ],
                    &[
                        &monospace_roff("minimap2-pb"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-pb"))
                    ],
                    &[
                        &monospace_roff("minimap2-hifi"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x map-hifi"))
                    ],
                    &[
                        &monospace_roff("minimap2-no-preset"),
                        &format!("minimap2 with no '{}' option", &monospace_roff("-x"))
                    ],
                    &[
                        &monospace_roff("minibwa"),
                        "minibwa map using default parameters"
                    ],
                    &[
                        &monospace_roff("rammap-sr"),
                        &format!("rammap (a minimap2-compatible aligner) with '{}' option", &monospace_roff("-x sr"))
                    ],
                    &[
                        &monospace_roff("rammap-lr-hq"),
                        &format!("rammap with '{}' option", &monospace_roff("-x lr:hq"))
                    ],
                    &[
                        &monospace_roff("rammap-ont"),
                        &format!("rammap with '{}' option", &monospace_roff("-x map-ont"))
                    ],
                    &[
                        &monospace_roff("rammap-pb"),
                        &format!("rammap with '{}' option", &monospace_roff("-x map-pb"))
                    ],
                    &[
                        &monospace_roff("rammap-hifi"),
                        &format!("rammap with '{}' option", &monospace_roff("-x map-hifi"))
                    ],
                    &[
                        &monospace_roff("rammap-no-preset"),
                        &format!("rammap with no '{}' option", &monospace_roff("-x"))
                    ],
                ])
            )))
            .option(Opt::new("PARAMS").long("--minimap2-params").help(&format!(
                "Extra parameters to provide to minimap2, \
        both indexing command (if used) and for \
        mapping. Note that usage of this parameter \
        has security implications if untrusted input \
        is specified. '{}' is always specified to minimap2. \
        [default: none]",
                &monospace_roff("-a")
            )))
            .flag(Flag::new().long("--minimap2-reference-is-index").help(
                "Treat reference as a minimap2 database, not as a FASTA file. [default: not set]",
            ))
            .option(Opt::new("PARAMS").long("--bwa-params").help(
                "Extra parameters to provide to BWA or BWA-MEM2. Note \
        that usage of this parameter has security \
        implications if untrusted input is specified. \
        [default: none]",
            ))
            .option(Opt::new("PARAMS").long("--strobealign-params").help(
                "Extra parameters to provide to strobealign. Note \
        that usage of this parameter has security \
        implications if untrusted input is specified. \
        [default: none]",
            ))
            .option(Opt::new("PARAMS").long("--minibwa-params").help(
                "Extra parameters to provide to the 'minibwa map' \
        command. Note that usage of this parameter has security \
        implications if untrusted input is specified. \
        [default: none]",
            ))
            .option(Opt::new("PARAMS").long("--rammap-params").help(&format!(
                "Extra parameters to provide to rammap. Note \
        that usage of this parameter has security \
        implications if untrusted input is specified. '{}' \
        is always specified to rammap. [default: none]",
                &monospace_roff("-x sr -a")
            )))
            .flag(Flag::new().long("--strobealign-use-index").help(
                "Use a pregenerated index (one that has been created with 'strobealign --create-index'). The --reference option should be specified as the original FASTA file i.e. 'ref.fna' not 'ref.fna.r100.sti' [default: not set]",
            )),
    )
}

fn add_thresholding_options(manual: Manual) -> Manual {
    manual.custom(
        Section::new("Alignment thresholding")
            .option(
                Opt::new("INT")
                    .long("--min-read-aligned-length")
                    .help(&format!(
                        "Exclude reads with smaller numbers of \
        aligned bases. {}",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-percent-identity")
                    .help(&format!(
                        "Exclude reads by overall percent \
        identity e.g. 95 for 95%. {}",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-aligned-percent")
                    .help(&format!(
                        "Exclude reads by percent aligned \
        bases e.g. 95 means 95% of the read's \
        bases must be aligned. {}",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("INT")
                    .long("--min-read-aligned-length-pair")
                    .help(&format!(
                        "Exclude pairs with smaller numbers of \
        aligned bases. \
        Implies --proper-pairs-only. {}",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-percent-identity-pair")
                    .help(&format!(
                        "Exclude pairs by overall percent \
                identity e.g. 95 for 95%. \
                Implies --proper-pairs-only. {}",
                        default_roff("0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-aligned-percent-pair")
                    .help(&format!(
                        "Exclude reads by percent aligned \
                bases e.g. 95 means 95% of the read's \
                bases must be aligned. \
                Implies --proper-pairs-only. {}",
                        default_roff("0")
                    )),
            )
            .option(Opt::new("INT").long("--min-mapq").help(&format!(
                "Exclude reads with a mapping quality \
                below this value. An INT between 0 and 254. \
                When thresholding pairs, both reads \
                are excluded if either has a MAPQ below this value. \
                If 0, \
                reads with MAPQ of 255 (mapping quality unavailable) \
                are excluded. {}",
                default_roff("not set")
            )))
            .flag(
                Flag::new()
                    .long("--proper-pairs-only")
                    .help("Require reads to be mapped as proper pairs. [default: not set]"),
            )
            .flag(
                Flag::new()
                    .long("--exclude-supplementary")
                    .help("Exclude supplementary alignments. [default: not set]"),
            )
            .flag(
                Flag::new()
                    .long("--include-secondary")
                    .help("Include secondary alignments. [default: not set]"),
            ),
    )
}

fn read_mapping_params_section() -> Section {
    Section::new("Read mapping parameters")
        .option(
            Opt::new("PATH ..")
                .short("-1")
                .help("Forward FASTA/Q file(s) for mapping. These may be gzipped or not."),
        )
        .option(
            Opt::new("PATH ..")
                .short("-2")
                .help("Reverse FASTA/Q file(s) for mapping. These may be gzipped or not."),
        )
        .option(Opt::new("PATH ..").short("-c").long("--coupled").help(
            "One or more pairs of forward and reverse \
        possibly gzipped FASTA/Q files for mapping in order \
        <sample1_R1.fq.gz> <sample1_R2.fq.gz> \
        <sample2_R1.fq.gz> <sample2_R2.fq.gz> ..",
        ))
        .option(
            Opt::new("PATH ..")
                .long("--interleaved")
                .help("Interleaved FASTA/Q files(s) for mapping. These may be gzipped or not."),
        )
        .option(
            Opt::new("PATH ..")
                .long("--single")
                .help("Unpaired FASTA/Q files(s) for mapping. These may be gzipped or not."),
        )
}

fn add_help_options(manual: Manual) -> Manual {
    manual
        .flag(
            Flag::new()
                .short("-h")
                .long("--help")
                .help("Output a short usage message. [default: not set]"),
        )
        .flag(
            Flag::new()
                .long("--full-help")
                .help("Output a full help message and display in 'man'. [default: not set]"),
        )
        .flag(Flag::new().long("--full-help-roff").help(
            "Output a full help message in raw ROFF format for \
        conversion to other formats. [default: not set]",
        ))
}
fn add_help_options_to_section(section: Section) -> Section {
    section
        .flag(
            Flag::new()
                .short("-h")
                .long("--help")
                .help("Output a short usage message. [default: not set]"),
        )
        .flag(
            Flag::new()
                .long("--full-help")
                .help("Output a full help message and display in 'man'. [default: not set]"),
        )
        .flag(Flag::new().long("--full-help-roff").help(
            "Output a full help message in raw ROFF format for \
        conversion to other formats. [default: not set]",
        ))
}

fn sharding_section() -> Section {
    Section::new("Sharding").flag(Flag::new().long("--sharded").help(&format!(
        "If {} was used: \
        Input BAM files are read-sorted alignments \
        of a set of reads mapped to multiple \
        reference contig sets. Choose the best \
        hit for each read pair. Otherwise if mapping was carried out: \
        Map reads to each reference, choosing the \
        best hit for each pair. [default: not set]",
        monospace_roff("-b/--bam-files")
    )))
}

fn faq_section() -> Section {
    Section::new("Frequently asked questions (FAQ)").paragraph(&format!(
        "{} CoverM makes use of \
        the system temporary directory (often {}) to store intermediate files. This can cause \
        problems if the amount of storage available there is small or used by many programs. \
        To fix, set the {} environment variable e.g. to set it to use the current directory: {}\n\
        \n\
        {} Either is fine, CoverM determines which is being used by virtue of being less than \
        or greater than 1.",
        bold("Can the temporary directory used be changed?"),
        monospace_roff("/tmp"),
        monospace_roff("TMPDIR"),
        monospace_roff("TMPDIR=. coverm genome <etc>"),
        bold(
            "For thresholding arguments e.g. --dereplication-ani and --min-read-percent-identity, \
        should a percentage (e.g 97%) or fraction (e.g. 0.97) be specified?"
        )
    ))
}

fn add_verbosity_flags(manual: Manual) -> Manual {
    manual
        .flag(
            Flag::new()
                .short("-v")
                .long("--verbose")
                .help("Print extra debugging information. [default: not set]"),
        )
        .flag(Flag::new().short("-q").long("--quiet").help(
            "Unless there is an error, do not print \
    log messages. [default: not set]",
        ))
}
fn add_verbosity_flags_to_section(section: Section) -> Section {
    section
        .flag(
            Flag::new()
                .short("-v")
                .long("--verbose")
                .help("Print extra debugging information. [default: not set]"),
        )
        .flag(Flag::new().long("--quiet").help(
            "Unless there is an error, do not print \
    log messages. [default: not set]",
        ))
}

pub fn filter_full_help() -> Manual {
    let mut manual = Manual::new("coverm filter")
        .about(format!(
            "Threshold alignments with insufficient identity (version {})",
            crate_version!()
        ))
        .author(Author::new(crate::AUTHOR).email("benjwoodcroft near gmail.com"))
        .description(
            "Only primary, non-supplementary alignments are considered, and output files \
        are grouped by reference, but not sorted by position.",
        )
        .option(
            Opt::new("PATH ..")
                .short("-b")
                .long("--bam-files")
                .help("Path to reference-sorted BAM file(s). [required]"),
        )
        .option(
            Opt::new("PATH ..")
                .short("-o")
                .long("--output-bam-files")
                .help(" Path to corresponding output file(s). [required]"),
        );
    manual = add_thresholding_options(manual);
    manual = manual.option(Opt::new("INT").short("-t").long("--threads").help(&format!(
        "Number of threads for output compression. {}",
        default_roff("1")
    )));
    manual = manual.flag(Flag::new().long("--inverse").help(
        "Only keep reads which are unmapped or \
                align below thresholds. Note that output \
                records may still be marked as mapped \
                if they do not meet the thresholds. When used with \
                --proper-pairs-only, only proper pairs which fail alignment \
                thresholds are output i.e. it does not \"invert\" the \
                proper pairs flag. [default: not set]",
    ));
    manual = add_verbosity_flags(manual);
    manual = add_help_options(manual);

    manual = manual.example(
        Example::new()
            .text("Filter a BAM file by removing alignments shorter than 50bp")
            .command(
                "coverm filter --bam-files input.bam --output-bam-files filtered.bam \
                --min-read-aligned-length 50",
            ),
    );
    manual = manual.example(
        Example::new()
            .text(
                "Filter inverse: Keep alignments that have <95% alignment identity \
            and those which do map at all. Note that the output BAM file will likely \
            records that are still mapped, but align with < 95% identity. Use 16 \
            threads for output compression",
            )
            .command(
                "coverm filter -b input.bam -o inverse_filtered.bam --inverse \
                --min-read-percent-identity 95 --threads 16",
            ),
    );

    manual = manual.custom(faq_section());

    manual
}

pub fn make_full_help() -> Manual {
    let mut manual = Manual::new("coverm make")
        .about(format!(
            "Generate BAM files through mapping (version: {})",
            crate_version!()
        ))
        .custom_synopsis_expansion("<REFERENCE> <READ_DEFINITION> <OUTPUT> ..")
        .author(Author::new(crate::AUTHOR).email("benjwoodcroft near gmail.com"))
        .description(
            "coverm make generates BAM files by read mapping a set of reads against \
        a reference FASTA database.\n\n",
        );

    manual = manual.custom(read_mapping_params_section());

    manual = manual.custom(
        Section::new("Reference").option(Opt::new("PATH").short("-r").long("--reference").help(
            &format!(
                "FASTA file of contigs e.g. concatenated \
                        genomes or metagenome assembly, or minimap2 \
                        index \
                        (with {}), \
                        strobealign index (with {}), \
                        or BWA index stem (with {}). May be gzip-compressed. [required]",
                monospace_roff("--minimap2-reference-is-index"),
                monospace_roff("--strobealign-use-index"),
                monospace_roff("-p bwa-mem/bwa-mem2"),
            ),
        )),
    );

    manual = add_mapping_options(manual);

    manual = manual.custom(
            Section::new("Output")    
        .option(Opt::new("DIR").short("-o").long("--output-directory").help(
            "Where generated BAM files will go. The directory will be created if it does not exist. [required]",
        ))
        .flag(
            Flag::new()
                .long("--discard-unmapped")
                .help("Exclude unmapped reads from cached BAM files. [default: not set]"),
        ));

    manual = manual.example(
        Example::new()
            .text(
                "Map pair of read files to the combined_genomes.fna reference, \
                storing sorted BAM files in output_dir/",
            )
            .command("coverm make -r combined_genomes.fna -1 read1.fq -2 read2.fq -o output_dir"),
    );

    let mut general_section = Section::new("General options").option(
        Opt::new("INT").short("-t").long("--threads").help(&format!(
            "Number of threads for mapping and sorting. {}",
            default_roff("1")
        )),
    );
    general_section = add_help_options_to_section(general_section);
    general_section = add_verbosity_flags_to_section(general_section);
    manual = manual.custom(general_section);

    manual = manual.custom(faq_section());

    manual
}

pub fn makedb_full_help() -> Manual {
    let mut manual = Manual::new("coverm makedb")
        .about(format!(
            "Generate mapping database(s) from reference FASTA files (version: {})",
            crate_version!()
        ))
        .custom_synopsis_expansion("-r <REFERENCE> -p <MAPPER> .. -o <OUTPUT_DIRECTORY>")
        .author(Author::new(crate::AUTHOR).email("benjwoodcroft near gmail.com"))
        .description(
            "coverm makedb pre-generates one or more mapping databases (indexes) from \
        reference genome or contig FASTA files. The generated database can then be \
        supplied to 'coverm contig' or 'coverm genome' (or 'coverm make') to avoid \
        regenerating the index each time reads are mapped.\n\n\
        The reference can be specified directly as one or more pre-concatenated FASTA \
        files (--reference), or as a set of genome FASTA files (--genome-fasta-files, \
        --genome-fasta-directory or --genome-fasta-list), exactly as in 'coverm genome'. \
        When genomes are given, they are concatenated on the fly into a single reference \
        FASTA (contigs are renamed to 'genome<separator>contig') that is written to the \
        output directory as 'coverm_concatenated_genomes.fna' and used to build the \
        database(s). The resulting database can then be used with 'coverm genome -s ~' \
        (the concatenation separator) to recover per-genome coverage.\n\n\
        When genomes are given, they may optionally be dereplicated (--dereplicate) \
        and/or quality-filtered by CheckM before concatenation, exactly as in \
        'coverm genome'.\n\n\
        In all cases the reference passed back to 'coverm contig'/'coverm genome' is the \
        generated database *file* (as reported by makedb, e.g. \
        'db_dir/ref.fna.minimap2-sr.mmi'), not the -o output directory, and the matching \
        '-p/--mapper' must be given (the default mapper is strobealign, so omitting it \
        feeds the index to the wrong mapper). For minimap2 databases, pass the generated \
        '.mmi' file as the reference together with '-p minimap2-sr' (or the preset used) \
        and '--minimap2-reference-is-index'. For BWA databases, pass the generated prefix \
        as the reference together with the matching '-p bwa-mem' or '-p bwa-mem2'. For \
        rammap databases, pass the generated '.idx' file as the reference together with \
        '-p rammap-sr' (or the preset used); rammap auto-detects the index. For \
        strobealign databases, the reference FASTA is copied into the output directory \
        next to the index (strobealign reads the sequences from it at mapping time); pass \
        that copied FASTA as the reference together with '-p strobealign' and \
        '--strobealign-use-index'.\n\n\
        Multiple '-p/--mapper' values may be specified to create several databases in \
        one invocation, one per mapper. Likewise, multiple references may be given, in \
        which case a database is created for each combination of reference and mapper.\n\n",
        );

    manual = manual.custom(
        Section::new("Input")
            .option(Opt::new("PATH ..").short("-r").long("--reference").help(
                "FASTA file(s) of contigs e.g. concatenated genomes or metagenome assembly. \
                May be gzip-compressed. [required unless a genome input below is specified]",
            ))
            .option(
                Opt::new("PATH ..")
                    .short("-f")
                    .long("--genome-fasta-files")
                    .help(
                        "Path(s) to FASTA files of each genome e.g. 'pathA/genome1.fna \
                        pathB/genome2.fa'. These are concatenated into a single reference \
                        before the database is built.",
                    ),
            )
            .option(
                Opt::new("PATH")
                    .short("-d")
                    .long("--genome-fasta-directory")
                    .help(
                        "Directory containing FASTA files of each genome, which are \
                        concatenated into a single reference before the database is built.",
                    ),
            )
            .option(Opt::new("PATH").long("--genome-fasta-list").help(
                "File containing FASTA file paths, one per line, of each genome to \
                concatenate into a single reference before the database is built.",
            ))
            .option(
                Opt::new("EXT")
                    .short("-x")
                    .long("--genome-fasta-extension")
                    .help(
                        "File extension of genomes in the directory specified with \
                        -d/--genome-fasta-directory. [default: fna]",
                    ),
            ),
    );

    manual = manual.custom(
        Section::new("Database type")
            .option(
                Opt::new("NAME ..")
                    .short("-p")
                    .long("--mapper")
                    .help(&format!(
                        "Kind(s) of database to generate, one per mapping software. Specify more \
                than once (or as a space-separated list) to generate several databases. \
                {}. One of: {}",
                        default_roff("minimap2-sr"),
                        bird_tool_utils::clap_utils::table_roff(&[
                            &["name", "description"],
                            &[
                                &monospace_roff("minimap2-sr"),
                                &format!(
                                    "minimap2 index built with '{}'",
                                    &monospace_roff("-x sr")
                                )
                            ],
                            &[
                                &monospace_roff("minimap2-lr-hq"),
                                &format!(
                                    "minimap2 index built with '{}'",
                                    &monospace_roff("-x lr:hq")
                                )
                            ],
                            &[
                                &monospace_roff("minimap2-ont"),
                                &format!(
                                    "minimap2 index built with '{}'",
                                    &monospace_roff("-x map-ont")
                                )
                            ],
                            &[
                                &monospace_roff("minimap2-pb"),
                                &format!(
                                    "minimap2 index built with '{}'",
                                    &monospace_roff("-x map-pb")
                                )
                            ],
                            &[
                                &monospace_roff("minimap2-hifi"),
                                &format!(
                                    "minimap2 index built with '{}'",
                                    &monospace_roff("-x map-hifi")
                                )
                            ],
                            &[
                                &monospace_roff("minimap2-no-preset"),
                                &format!(
                                    "minimap2 index built with no '{}' option",
                                    &monospace_roff("-x")
                                )
                            ],
                            &[&monospace_roff("bwa-mem"), "BWA index (bwa index)"],
                            &[
                                &monospace_roff("bwa-mem2"),
                                "BWA-MEM2 index (bwa-mem2 index)"
                            ],
                            &[&monospace_roff("minibwa"), "minibwa index (minibwa index)",],
                            &[
                                &monospace_roff("rammap-sr"),
                                &format!(
                                    "rammap index ('{}'), built with the '{}' preset",
                                    &monospace_roff("rammap --dump-index"),
                                    &monospace_roff("-x sr")
                                ),
                            ],
                            &[
                                &monospace_roff("rammap-lr-hq"),
                                &format!(
                                    "rammap index built with the '{}' preset",
                                    &monospace_roff("-x lr:hq")
                                ),
                            ],
                            &[
                                &monospace_roff("rammap-ont"),
                                &format!(
                                    "rammap index built with the '{}' preset",
                                    &monospace_roff("-x map-ont")
                                ),
                            ],
                            &[
                                &monospace_roff("rammap-pb"),
                                &format!(
                                    "rammap index built with the '{}' preset",
                                    &monospace_roff("-x map-pb")
                                ),
                            ],
                            &[
                                &monospace_roff("rammap-hifi"),
                                &format!(
                                    "rammap index built with the '{}' preset",
                                    &monospace_roff("-x map-hifi")
                                ),
                            ],
                            &[
                                &monospace_roff("rammap-no-preset"),
                                &format!(
                                    "rammap index built with no '{}' option",
                                    &monospace_roff("-x")
                                ),
                            ],
                            &[
                                &monospace_roff("strobealign"),
                                &format!(
                            "strobealign index (strobealign --create-index). The reference \
                            FASTA is copied into the output directory alongside the index, \
                            since strobealign requires it at mapping time. Use it via '{}'",
                            &monospace_roff("--strobealign-use-index")
                        )
                            ],
                        ])
                    )),
            )
            .option(Opt::new("PARAMS").long("--minimap2-params").help(
                "Extra parameters to provide to the minimap2 indexing command. \
                Note that usage of this parameter has security implications if \
                untrusted input is specified. [default: none]",
            ))
            .option(Opt::new("PARAMS").long("--bwa-params").help(
                "Extra parameters to provide to the BWA or BWA-MEM2 indexing \
                command. Note that usage of this parameter has security \
                implications if untrusted input is specified. [default: none]",
            ))
            .option(Opt::new("PARAMS").long("--strobealign-params").help(
                "Extra parameters to provide to the 'strobealign --create-index' \
                command. Strobealign indexes are read-length specific: set the \
                canonical read length with e.g. '-r 150', or estimate it from an \
                example read dataset by passing a reads file (e.g. 'reads.fq'). \
                Note that usage of this parameter has security implications if \
                untrusted input is specified. [default: none]",
            )),
    );

    manual = manual.custom(Section::new("Output").option(
        Opt::new("DIR").short("-o").long("--output-directory").help(
            "Where the generated database(s) will be written. The directory will \
                    be created if it does not exist. [required]",
        ),
    ));

    manual = manual.example(
        Example::new()
            .text("Generate a short-read minimap2 database from a pre-concatenated FASTA")
            .command("coverm makedb -r combined_genomes.fna -p minimap2-sr -o db_dir"),
    );
    manual = manual.example(
        Example::new()
            .text(
                "Generate a short-read minimap2 database from a directory of genomes \
                (concatenated on the fly)",
            )
            .command("coverm makedb -d genomes_directory/ -x fna -p minimap2-sr -o db_dir"),
    );
    manual = manual.example(
        Example::new()
            .text(
                "Use the database built from genomes to calculate per-genome coverage \
                (the concatenation separator is '~')",
            )
            .command(
                "coverm genome -p minimap2-sr -r db_dir/coverm_concatenated_genomes.fna.minimap2-sr.mmi \
                --minimap2-reference-is-index -s '~' -1 read1.fq -2 read2.fq",
            ),
    );
    manual = manual.example(
        Example::new()
            .text("Use the generated minimap2 database when calculating contig coverage")
            .command(
                "coverm contig -p minimap2-sr -r db_dir/combined_genomes.fna.minimap2-sr.mmi \
                --minimap2-reference-is-index -1 read1.fq -2 read2.fq",
            ),
    );
    manual = manual.example(
        Example::new()
            .text("Generate several databases at once, one per mapper")
            .command("coverm makedb -r combined_genomes.fna -p minimap2-sr minimap2-ont bwa-mem -o db_dir"),
    );
    manual = manual.example(
        Example::new()
            .text(
                "Generate a strobealign database for 150bp reads (the reference FASTA is \
                copied into db_dir alongside the index)",
            )
            .command(
                "coverm makedb -r combined_genomes.fna -p strobealign \
                --strobealign-params '-r 150' -o db_dir",
            ),
    );
    manual = manual.example(
        Example::new()
            .text("Dereplicate a directory of genomes before building the database")
            .command(
                "coverm makedb -d genomes_directory/ -x fna --dereplicate \
                -p minimap2-sr -o db_dir",
            ),
    );

    let mut derep_section = Section::new("DEREPLICATION / GENOME CLUSTERING").flag(
        Flag::new().long("--dereplicate").help(
            "Dereplicate the input genomes via average nucleotide identity (ANI) before \
            concatenating them into the reference - choose a genome to represent all \
            within a small distance, using skani for clustering by default. \
            Dereplication occurs transparently through the Galah method \
            (https://github.com/wwood/galah) [default: not set]",
        ),
    );
    derep_section =
        galah::cluster_argument_parsing::add_dereplication_filtering_parameters_to_section(
            derep_section,
        );
    derep_section =
        galah::cluster_argument_parsing::add_dereplication_clustering_parameters_to_section(
            derep_section,
            &COVERM_CLUSTER_COMMAND_DEFINITION,
        );
    derep_section = galah::cluster_argument_parsing::add_dereplication_output_parameters_to_section(
        derep_section,
        &COVERM_CLUSTER_COMMAND_DEFINITION,
    );
    manual = manual.custom(derep_section);

    let mut general_section = Section::new("General options").option(
        Opt::new("INT").short("-t").long("--threads").help(&format!(
            "Number of threads used to generate the database(s). {}",
            default_roff("1")
        )),
    );
    general_section = add_help_options_to_section(general_section);
    general_section = add_verbosity_flags_to_section(general_section);
    manual = manual.custom(general_section);

    manual = manual.custom(faq_section());

    manual
}

pub fn contig_full_help() -> Manual {
    let mut manual = Manual::new("coverm contig")
        .about(format!("Calculate read coverage per-contig (version {})",crate_version!()))
        .custom_synopsis_expansion("<MAPPING_INPUT> ..")
        .author(Author::new(crate::AUTHOR).email("benjwoodcroft near gmail.com"))
        .description("coverm contig calculates the coverage of a set of reads on a set of contigs.\n\n\
        This process can be undertaken in several ways, for instance by specifying BAM files or raw reads as input, \
        using different mapping programs, thresholding read alignments, using different methods of calculating coverage \
        and printing the calculated coverage in various formats.\n\
        \n\
        The source code for CoverM is available at https://github.com/wwood/CoverM");

    manual = manual.custom(
        read_mapping_params_section().option(
            Opt::new("PATH")
                .short("-b")
                .long("--bam-files")
                .help(&format!(
                    "Path to BAM file(s). These must be \
                reference sorted (e.g. with samtools sort) \
                unless {} is specified, in which \
                case they must be read name sorted (e.g. \
                with {}). When specified, no read mapping algorithm is undertaken.",
                    monospace_roff("--sharded"),
                    monospace_roff("samtools sort -n"),
                )),
        ),
    );

    manual = manual.custom(
        Section::new("Reference").option(Opt::new("PATH").short("-r").long("--reference").help(
            &format!(
                "FASTA file of contigs e.g. concatenated \
                    genomes or metagenome assembly (gzip or uncompressed), or minimap2 \
                    index \
                    (with {}), \
                    strobealign index (with {}), \
                    or BWA index stem (with {}). \
                    If multiple references FASTA files are \
                    provided and {} is specified, \
                    then reads will be mapped to references \
                    separately as sharded BAMs. [required unless {} is specified]",
                monospace_roff("--minimap2-reference-is-index"),
                monospace_roff("--strobealign-use-index"),
                monospace_roff("-p bwa-mem/bwa-mem2"),
                monospace_roff("--sharded"),
                monospace_roff("-b/--bam-files")
            ),
        )),
    );

    manual = manual.custom(sharding_section());
    manual = add_mapping_options(manual);
    manual = add_thresholding_options(manual);

    manual = manual.custom(
        Section::new("Coverage calculation options")
            .option(Opt::new("METHOD").short("-m").long("--methods").help(
                &format!("Method(s) for calculating coverage {}. A more thorough description of the different methods is available at\n\
                https://github.com/wwood/CoverM#calculation-methods but briefly:\n\
                {}",
                default_roff("mean"),
                bird_tool_utils::clap_utils::table_roff(&[
                    &["method","description"],
                    &[&monospace_roff("mean"), "(default) Average number of aligned reads overlapping each position on the contig"],
                    &[&monospace_roff("trimmed_mean"), &format!("Average number of aligned reads overlapping each position after removing the most deeply and shallow-ly covered positions. See {}/{} to adjust.",
                        &monospace_roff("--trim-min"),
                        &monospace_roff("--trim-max"))],
                    &[&monospace_roff("strobealign-aemb"), 
                        &format!("Mean coverage as estimated by strobealign --aemb, which is faster than the {} method, giving similar but not identical values. See https://github.com/ksahlin/strobealign for details. Cannot currently be used with other methods.",
                        &monospace_roff("mean"))],
                    &[&monospace_roff("coverage_histogram"), "Histogram of coverage depths"],
                    &[&monospace_roff("covered_fraction"), "Proportion of bases covered by 1 or more reads"],
                    &[&monospace_roff("covered_bases"), "Number of bases covered by 1 or more reads"],
                    &[&monospace_roff("variance"), "Variance of coverage depths"],
                    &[&monospace_roff("length"), "Length of each contig in base pairs"],
                    &[&monospace_roff("count"), "Number of reads aligned to each contig. Note that supplementary alignments are not counted."],
                    &[&monospace_roff("metabat"), "(\"MetaBAT adjusted coverage\") Coverage as defined in Kang et al 2015 https://doi.org/10.7717/peerj.1165"],
                    &[&monospace_roff("reads_per_base"), "Number of reads aligned divided by the length of the contig"],
                    &[&monospace_roff("anir"), "Average BLAST-like identity of mapped reads"],
                    &[&monospace_roff("rpkm"), "Reads mapped per kilobase of contig, per million mapped reads"],
                    &[&monospace_roff("tpm"), "Transcripts Per Million as described in Li et al 2010 https://doi.org/10.1093/bioinformatics/btp692"],
                ]),
            )))
            .option(Opt::new("FRACTION").long("--min-covered-fraction").help(
                &format!("Contigs with less covered bases than this are \
                reported as having zero coverage. \
                {}", default_roff("0"))
            ))
            .option(Opt::new("INT").long("--contig-end-exclusion").help(
                &format!("Exclude bases at the ends of reference \
                sequences from calculation {}",
                default_roff("75"))
            ))
            .option(Opt::new("FRACTION").long("--trim-min").help(
                &format!("Remove this smallest fraction of positions \
                when calculating trimmed_mean {}",
                default_roff("5"))
            ))
            .option(Opt::new("FRACTION").long("--trim-max").help(
                &format!("Maximum fraction for trimmed_mean \
                calculations {}", default_roff("95"))
            )),
    );

    manual = manual.custom(
        Section::new("Per-gene coverage")
            .option(Opt::new("PATH").long("--gff").help(&format!(
                "GFF (or GTF) file defining genes/features. When specified, \
                coverage is reported once per feature (using the chosen \
                {}) rather than once per contig, with the gene identifier and \
                its contig in the first two output columns. Each non-comment \
                line is treated as a separate feature; the reported identifier \
                is taken from the {}, {}, {} or {} attribute. Reads are \
                assigned to a feature for read-count based methods (e.g. {}, \
                {}, {}) when their leftmost mapped position falls within the \
                feature's coordinates. Note that {} applies to the ends of \
                each feature, so {} may be appropriate for short features. \
                Cannot be used with the {} or {} methods. (--genes-file is an \
                alias). [default: not used]",
                monospace_roff("--methods"),
                monospace_roff("ID"),
                monospace_roff("locus_tag"),
                monospace_roff("gene_id"),
                monospace_roff("Name"),
                monospace_roff("count"),
                monospace_roff("rpkm"),
                monospace_roff("tpm"),
                monospace_roff("--contig-end-exclusion"),
                monospace_roff("--contig-end-exclusion 0"),
                monospace_roff("metabat"),
                monospace_roff("strobealign-aemb"),
            )))
            .option(Opt::new("TYPE").long("--gff-feature-type").help(
                "Only use features of this type (i.e. the third column of \
                the GFF) when --gff is specified. [default: use all features]",
            )),
    );

    manual = manual.custom(
        Section::new("Output")
            .option(Opt::new("FILE").short("-o").long("--output-file").help(
                "Output coverage values to this file, or '-' for STDOUT. \
                [default: output to STDOUT]",
            ))
            .option(Opt::new("FORMAT").long("--output-format").help(
                "Shape of output: 'sparse' for long format, \
    'dense' for species-by-site. \
    [default: dense]",
            ))
            .flag(Flag::new().long("--no-zeros").help(
                "Omit printing of genomes that have zero \
        coverage. [default: not set]",
            ))
            .option(
                Opt::new("DIRECTORY")
                    .long("--cache-unfiltered-bam-directory")
                    .help(
                        "Output BAM files generated during \
                alignment to this directory (--bam-file-cache-directory is an alias). The directory may or may not exist. Note that \
                BAM files in this directory contain all mappings, including those that later \
                are excluded by alignment thresholding (e.g. --min-read-percent-identity) or \
                genome-wise thresholding (e.g. --min-covered-fraction). \
                [default: not used]",
                    ),
            )
            .option(
                Opt::new("FILE")
                    .long("--cache-unfiltered-bam-files")
                    .help(
                        "Output BAM files generated during alignment to these files. The order of files should correspond to: single-ended reads (-s), -1/-2, --coupled, --interleaved. [default: not used]",
                    ),
            )
            .flag(
                Flag::new()
                    .long("--discard-unmapped")
                    .help("Exclude unmapped reads from cached BAM files. [default: not set]"),
            ),
    );

    manual = manual.example(
        Example::new()
            .text("Calculate mean coverage from reads and assembly")
            .command(
                "coverm contig --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna",
            ),
    );
    manual = manual.example(
        Example::new()
            .text(
                "Calculate MetaBAT adjusted coverage from a sorted BAM file, saving \
                the unfiltered BAM files in the saved_bam_files folder",
            )
            .command(
                "coverm contig --methods metabat --bam-files my.bam \
                --cache-unfiltered-bam-directory saved_bam_files",
            ),
    );

    let mut general_section = Section::new("General options").option(
        Opt::new("INT")
            .short("-t")
            .long("--threads")
            .help("Number of threads for mapping, sorting and reading. [default: 1]"),
    );
    general_section = add_help_options_to_section(general_section);
    general_section = add_verbosity_flags_to_section(general_section);
    manual = manual.custom(general_section);

    manual = manual.custom(faq_section());

    manual
}

pub fn genome_full_help() -> Manual {
    let mut manual = Manual::new("coverm genome")
        .about(format!("Calculate read coverage per-genome (version {})",crate_version!()))
        .custom_synopsis_expansion("<GENOME_DESCRIPTION> <MAPPING_INPUT> ..")
        .author(Author::new(crate::AUTHOR).email("benjwoodcroft near gmail.com"))

        .description("coverm genome calculates the coverage of a set of reads on a set of genomes.\n\n\
            This process can be undertaken in several ways, for instance by specifying BAM files or raw reads as input, \
            defining genomes in different input formats, dereplicating genomes before mapping, \
            using different mapping programs, thresholding read alignments, using different methods of calculating coverage \
            and printing the calculated coverage in various formats.\n\
            \n\
            The source code for CoverM is available at https://github.com/wwood/CoverM");

    manual = manual.custom(
        read_mapping_params_section().option(
            Opt::new("PATH")
                .short("-b")
                .long("--bam-files")
                .help(&format!(
                    "Path to BAM file(s). These must be \
                reference sorted (e.g. with samtools sort) \
                unless {} is specified, in which \
                case they must be read name sorted (e.g. \
                with {}). When specified, no read mapping algorithm is undertaken.",
                    monospace_roff("--sharded"),
                    monospace_roff("samtools sort -n")
                )),
        ),
    );

    manual = manual.custom(
        bird_tool_utils::clap_utils::add_genome_specification_to_section(
            Section::new("Genome definition"))
            .option(
                Opt::new("PATH").short("-r").long("--reference").help(
                    &format!("FASTA file of contigs e.g. concatenated \
                    genomes or metagenome assembly (gzip or uncompressed), or minimap2 \
                    index \
                    (with {}), \
                    strobealign index (with {}), \
                    or BWA index stem (with {}). \
                    If multiple reference FASTA files are \
                    provided and {} is specified, \
                    then reads will be mapped to references \
                    separately as sharded BAMs. {}: If genomic FASTA files are \
                    specified elsewhere (e.g. with {} or {}), then {} is not needed as a reference FASTA file can be derived \
                    by concatenating input genomes. In these situations, {} can \
                    be optionally specified if an alternate reference sequence set is desired.",
                    monospace_roff("--minimap2-reference-is-index"),
                    monospace_roff("--strobealign-use-index"),
                    monospace_roff("-p bwa-mem/bwa-mem2"),
                    monospace_roff("--sharded"),
                    bold("NOTE"),
                    monospace_roff("--genome-fasta-files"),
                    monospace_roff("--genome-fasta-directory"),
                    //monospace_roff("--bam-files"),
                    monospace_roff("--reference"),
                    monospace_roff("--reference"),
                )
            )
            )
            .option(
                Opt::new("CHARACTER")
                    .short("-s")
                    .long("--separator")
                    .help(
                        &format!("This character separates genome names from contig names in the reference file. Requires {}. \
                    [default: unspecified]", monospace_roff("--reference")))
            )
            .flag(
                Flag::new()
                    .long("--single-genome")
                    .help(&format!("All contigs are from the same genome. Requires {}. [default: not set]", monospace_roff("--reference")))
            )
            .option(
                Opt::new("FILE")
                    .long("--genome-definition")
                    .help(&format!("File containing list of \
                    genome_name<tab>contig lines to define the genome of each contig. Requires {}. [default: not set]", monospace_roff("--reference")))
            )
            .flag(
                Flag::new()
                    .long("--use-full-contig-names")
                    .help("Specify that the input BAM files have been generated with mapping software that includes the full name of each contig \
                    in the reference definition (i.e. characters after the space), so when reading in genomes, record contig names as such.")
            )
        );

    let mut derep_section = Section::new("DEREPLICATION / GENOME CLUSTERING").flag(
        Flag::new().long("--dereplicate").help(
            "Do genome dereplication via average nucleotide \
        identity (ANI) - choose a genome to represent \
        all within a small distance, using skani for \
        clustering by default. When this flag is used, \
        dereplication occurs transparently through the \
        Galah method (https://github.com/wwood/galah) [default: not set]",
        ),
    );
    derep_section =
        galah::cluster_argument_parsing::add_dereplication_filtering_parameters_to_section(
            derep_section,
        );
    derep_section =
        galah::cluster_argument_parsing::add_dereplication_clustering_parameters_to_section(
            derep_section,
            &COVERM_CLUSTER_COMMAND_DEFINITION,
        );
    derep_section = galah::cluster_argument_parsing::add_dereplication_output_parameters_to_section(
        derep_section,
        &COVERM_CLUSTER_COMMAND_DEFINITION,
    );
    // not yet implemented
    //galah::cluster_argument_parsing::add_dereplication_output_parameters_to_section(derep_section);
    manual = manual.custom(derep_section);

    manual = manual.custom(sharding_section().flag(
        Flag::new().long("--exclude-genomes-from-deshard").help(
            "Ignore genomes whose name appears in this newline-separated \
                file when combining shards. [default: not set]",
        ),
    ));

    manual = add_mapping_options(manual);

    manual = add_thresholding_options(manual);

    manual = manual.custom(
        Section::new("Coverage calculation options")
            .option(Opt::new("METHOD").short("-m").long("--methods").help(
                &format!("Method(s) for calculating coverage {}. A more thorough description of the different methods is available at\n\
                https://github.com/wwood/CoverM#calculation-methods but briefly:\n\
                {}",
                default_roff("relative_abundance"),
                bird_tool_utils::clap_utils::table_roff(&[
                    &["method","description"],
                    &[&monospace_roff("relative_abundance"), "(default) Percentage relative abundance of each genome, and the unmapped read percentage"],
                    &[&monospace_roff("mean"), "Average number of aligned reads overlapping each position on the genome"],
                    &[&monospace_roff("trimmed_mean"), &format!("Average number of aligned reads overlapping each position after removing the most deeply and shallow-ly covered positions. See {}/{} to adjust.",
                        &monospace_roff("--trim-min"),
                        &monospace_roff("--trim-max"))],    
                    &[&monospace_roff("coverage_histogram"), "Histogram of coverage depths"],
                    &[&monospace_roff("covered_fraction"), "Proportion of bases covered by 1 or more reads"],
                    &[&monospace_roff("covered_bases"), "Number of bases covered by 1 or more reads"],
                    &[&monospace_roff("variance"), "Variance of coverage depths"],
                    &[&monospace_roff("length"), "Length of each genome in base pairs"],
                    &[&monospace_roff("count"), "Number of reads aligned to each genome. Note that supplementary alignments are not counted."],
                    &[&monospace_roff("reads_per_base"), "Number of reads aligned divided by the length of the genome"],
                    &[&monospace_roff("anir"), "Average BLAST-like identity of mapped reads"],
                    &[&monospace_roff("rpkm"), "Reads mapped per kilobase of genome, per million mapped reads"],
                    &[&monospace_roff("tpm"), "Transcripts Per Million as described in Li et al 2010 https://doi.org/10.1093/bioinformatics/btp692"],
                ])
            )))
            .option(Opt::new("FRACTION").long("--min-covered-fraction").help(
                &format!("Genomes with less covered bases than this are \
                reported as having zero coverage. \
                {}", default_roff("10"))
            ))
            .option(Opt::new("INT").long("--contig-end-exclusion").help(
                &format!("Exclude bases at the ends of reference \
                sequences from calculation {}",
                default_roff("75"))
            ))
            .option(Opt::new("FRACTION").long("--trim-min").help(
                &format!("Remove this smallest fraction of positions \
                when calculating trimmed_mean {}",
                default_roff("5"))
            ))
            .option(Opt::new("FRACTION").long("--trim-max").help(
                &format!("Maximum fraction for trimmed_mean \
                calculations {}", default_roff("95"))
            )),
    );

    manual = manual.custom(
        Section::new("Per-gene coverage")
            .option(Opt::new("PATH").long("--gff").help(&format!(
                "GFF (or GTF) file defining genes/features. When specified, \
                coverage is reported once per feature (using the chosen \
                {}) rather than once per genome, with the gene identifier, its \
                contig and its genome in the first three output columns. The \
                genome each feature belongs to is determined from the genome \
                definition supplied on the command line (which is still \
                required), and features whose contig is not assigned to a \
                genome are omitted. The reported identifier is taken from the \
                {}, {}, {} or {} attribute. The contig names in the GFF must \
                match the reference sequence names in the BAM/reference; note \
                that when mapping to multiple {} the contigs are renamed to \
                {}, so the GFF should use matching names. Reads are assigned \
                to a feature for read-count based methods when their leftmost \
                mapped position falls within the feature's coordinates. Note \
                that {} applies to the ends of each feature, so {} may be \
                appropriate for short features. Cannot be used with {}. \
                (--genes-file is an alias). [default: not used]",
                monospace_roff("--methods"),
                monospace_roff("ID"),
                monospace_roff("locus_tag"),
                monospace_roff("gene_id"),
                monospace_roff("Name"),
                monospace_roff("--genome-fasta-files"),
                monospace_roff("genome~contig"),
                monospace_roff("--contig-end-exclusion"),
                monospace_roff("--contig-end-exclusion 0"),
                monospace_roff("--dereplicate"),
            )))
            .option(Opt::new("TYPE").long("--gff-feature-type").help(
                "Only use features of this type (i.e. the third column of \
                the GFF) when --gff is specified. [default: use all features]",
            )),
    );

    manual = manual.custom(
        Section::new("Output")
            .option(Opt::new("FILE").short("-o").long("--output-file").help(
                "Output coverage values to this file, or '-' for STDOUT. \
                [default: output to STDOUT]",
            ))
            .option(Opt::new("FORMAT").long("--output-format").help(&format!(
                "Shape of output: 'sparse' for long format, \
            'dense' for species-by-site. {}",
                default_roff("dense")
            )))
            .flag(Flag::new().long("--no-zeros").help(
                "Omit printing of genomes that have zero \
            coverage. [default: not set]",
            ))
            .option(
                Opt::new("DIRECTORY")
                    .long("--cache-unfiltered-bam-directory")
                    .help(
                        "Output BAM files generated during \
                alignment to this directory (--bam-file-cache-directory is an alias). The directory may or may not exist. Note that \
                BAM files in this directory contain all mappings, including those that later \
                are excluded by alignment thresholding (e.g. --min-read-percent-identity) or \
                genome-wise thresholding (e.g. --min-covered-fraction). \
                [default: not used]",
                    ),
            )
            .option(
                Opt::new("FILE")
                    .long("--cache-unfiltered-bam-files")
                    .help(
                        "Output BAM files generated during alignment to these files. The order of files should correspond to: single-ended reads (-s), -1/-2, --coupled, --interleaved. [default: not used]",
                    ),
            )
            .flag(
                Flag::new()
                    .long("--discard-unmapped")
                    .help("Exclude unmapped reads from cached BAM files. [default: not set]"),
            ),
    );

    manual = manual.example(
        Example::new()
            .text("Map paired reads to 2 genomes, and output relative abundances to output.tsv")
            .command(
                "coverm genome --coupled read1.fastq.gz read2.fastq.gz \
            --genome-fasta-files genome1.fna genome2.fna -o output.tsv",
            ),
    );
    manual = manual.example(
        Example::new()
            .text(
                "Calculate coverage of genomes defined as .fna files in \
                genomes_directory/ from a sorted BAM file",
            )
            .command(
                "coverm genome --bam-files my.bam --genome-fasta-directory genomes_directory/",
            ),
    );
    manual = manual.example(
        Example::new()
            .text("Dereplicate genomes at 99% ANI before mapping unpaired reads")
            .command(
                "coverm genome --genome-fasta-directory genomes/ --dereplicate \
                --single single_reads.fq.gz",
            ),
    );

    let mut general_section = Section::new("General options").option(
        Opt::new("INT").short("-t").long("--threads").help(&format!(
            "Number of threads for mapping, sorting and reading. {}",
            default_roff("1")
        )),
    );
    general_section = add_help_options_to_section(general_section);
    general_section = add_verbosity_flags_to_section(general_section);
    manual = manual.custom(general_section);

    manual = manual.custom(faq_section());

    manual
}

/// Add the `--dereplicate` flag plus the Galah dereplication tuning arguments and
/// CheckM-based genome quality filtering arguments to a subcommand, mirroring the
/// options offered by `coverm genome`. These operate on the set of genomes given
/// via `--genome-fasta-*`, so they conflict with a pre-built `--reference`.
fn add_dereplication_arguments(command: Command) -> Command {
    command
        .arg(
            Arg::new("dereplicate")
                .long("dereplicate")
                .conflicts_with("reference")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("dereplication-ani")
                .long("dereplication-ani")
                .value_parser(clap::value_parser!(f32))
                .default_value(galah::DEFAULT_ANI),
        )
        .arg(
            Arg::new("dereplication-prethreshold-ani")
                .long("dereplication-prethreshold-ani")
                .value_parser(clap::value_parser!(f32))
                .default_value(galah::DEFAULT_PRETHRESHOLD_ANI),
        )
        .arg(
            Arg::new("dereplication-quality-formula")
                .long("dereplication-quality-formula")
                .value_parser([
                    "completeness-4contamination",
                    "completeness-5contamination",
                    "Parks2020_reduced",
                    "dRep",
                ])
                .default_value(galah::DEFAULT_QUALITY_FORMULA),
        )
        .arg(
            Arg::new("dereplication-precluster-method")
                .long("dereplication-precluster-method")
                .value_parser(galah::PRECLUSTER_METHODS)
                .default_value(galah::DEFAULT_PRECLUSTER_METHOD),
        )
        .arg(
            Arg::new("dereplication-cluster-method")
                .long("dereplication-cluster-method")
                .value_parser(galah::CLUSTER_METHODS)
                .default_value(galah::DEFAULT_CLUSTER_METHOD),
        )
        .arg(
            Arg::new("dereplication-aligned-fraction")
                .long("dereplication-aligned-fraction")
                .default_value(galah::DEFAULT_ALIGNED_FRACTION)
                .value_parser(clap::value_parser!(f32)),
        )
        .arg(
            Arg::new("dereplication-small-genomes")
                .long("dereplication-small-genomes")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("dereplication-cluster-contigs")
                .long("dereplication-cluster-contigs")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("dereplication-small-contigs")
                .long("dereplication-small-contigs")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("dereplication-large-contigs")
                .long("dereplication-large-contigs")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("dereplication-fragment-length")
                .long("dereplication-fragment-length")
                .default_value(galah::DEFAULT_FRAGMENT_LENGTH)
                .value_parser(clap::value_parser!(u32)),
        )
        .arg(
            Arg::new("dereplication-low-memory")
                .long("dereplication-low-memory")
                .action(clap::ArgAction::SetTrue)
                .conflicts_with("dereplication-reference-genomes")
                .conflicts_with("dereplication-reference-genomes-list"),
        )
        .arg(
            Arg::new("dereplication-reference-genomes")
                .long("dereplication-reference-genomes")
                .num_args(1..)
                .conflicts_with("dereplication-low-memory")
                .conflicts_with("dereplication-reference-genomes-list"),
        )
        .arg(
            Arg::new("dereplication-reference-genomes-list")
                .long("dereplication-reference-genomes-list")
                .conflicts_with("dereplication-low-memory")
                .conflicts_with("dereplication-reference-genomes"),
        )
        .arg(
            Arg::new("dereplication-output-cluster-definition")
                .long("dereplication-output-cluster-definition"),
        )
        .arg(
            Arg::new("dereplication-output-representative-fasta-directory")
                .long("dereplication-output-representative-fasta-directory"),
        )
        .arg(
            Arg::new("dereplication-output-representative-fasta-directory-copy")
                .long("dereplication-output-representative-fasta-directory-copy"),
        )
        .arg(
            Arg::new("dereplication-output-representative-list")
                .long("dereplication-output-representative-list"),
        )
        .arg(
            Arg::new("checkm-tab-table")
                .long("checkm-tab-table")
                .conflicts_with("reference"),
        )
        .arg(
            Arg::new("checkm2-quality-report")
                .long("checkm2-quality-report")
                .conflicts_with("reference"),
        )
        .arg(
            Arg::new("genome-info")
                .long("genome-info")
                .conflicts_with("reference")
                .conflicts_with("checkm-tab-table"),
        )
        .arg(Arg::new("min-completeness").long("min-completeness"))
        .arg(Arg::new("max-contamination").long("max-contamination"))
        .arg(
            Arg::new("run-checkm2")
                .long("run-checkm2")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(Arg::new("checkm2-db-path").long("checkm2-db-path"))
}

pub fn build_cli() -> Command {
    // specify static lazily because need to define it at runtime.
    lazy_static! {
        static ref CONTIG_HELP: String = format!(
            "
                            {}
              {}

{}

  coverm contig --coupled read1.fastq.gz read2.fastq.gz --reference assembly.fna

{}

  coverm contig --methods metabat --bam-files my.bam
    --cache-unfiltered-bam-directory saved_bam_files

See coverm contig --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint("coverm contig"),
            ansi_term::Colour::Green.paint("Calculate coverage of individual contigs"),
            ansi_term::Colour::Purple
                .paint("Example: Calculate mean coverage from reads and assembly:"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate MetaBAT adjusted coverage from a sorted BAM file, saving
the unfiltered BAM files in the saved_bam_files folder:"
            )
        );
        static ref GENOME_HELP: String = format!(
            "
                            {}
               {}

{}

  coverm genome --coupled read1.fastq.gz read2.fastq.gz
    --genome-fasta-files genome1.fna genome2.fna -o output.tsv 

{}

  coverm genome --bam-files my.bam --genome-fasta-directory genomes_directory/

{}

  coverm genome --genome-fasta-directory genomes/ --dereplicate
    --single single_reads.fq.gz

See coverm genome --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint("coverm genome"),
            ansi_term::Colour::Green.paint("Calculate coverage of individual genomes"),
            ansi_term::Colour::Purple.paint(
                "Example: Map paired reads to 2 genomes, and output relative abundances\n\
                to output.tsv:"
            ),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate coverage of genomes defined as .fna files in\n\
                genomes_directory/ from a sorted BAM file:"
            ),
            ansi_term::Colour::Purple
                .paint("Example: Dereplicate genomes at 99% ANI before mapping unpaired reads:"),
        );
        static ref FILTER_HELP: String = format!(
            "
                            {}
                     {}

{}

  coverm filter --bam-files input.bam --output-bam-files filtered.bam
    --min-read-aligned-length 50

{}

  coverm filter -b input.bam -o inverse_filtered.bam --inverse
    --min-read-percent-identity 95 --threads 16

See coverm filter --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint("coverm filter"),
            ansi_term::Colour::Green.paint("Filter BAM file alignments"),
            ansi_term::Colour::Purple
                .paint("Example: Filter a BAM file by removing alignments shorter than 50bp:"),
            ansi_term::Colour::Purple.paint(
                "Example: Filter inverse: Keep alignments that have <95% alignment identity\n\
                 and those which do map at all. Note that the output BAM file will likely\n\
                 records that are still mapped, but align with < 95% identity. Use 16\n\
                 threads for output compression:"
            ),
        );
        static ref MAKE_HELP: String = format!(
            "
                            {}
                     {}

{}

  coverm make -r combined_genomes.fna -1 read1.fq -2 read2.fq -o output_dir

See coverm make --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint("coverm make"),
            ansi_term::Colour::Green.paint("Generate BAM files through mapping"),
            ansi_term::Colour::Purple.paint(
                "Example: Map pair of read files to the combined_genomes.fna reference,\n\
                storing sorted BAM files in output_dir/"
            ),
        );
        static ref MAKEDB_HELP: String = format!(
            "
                            {}
                     {}

{}

  coverm makedb -d genomes_directory/ -x fna -p minimap2-sr -o db_dir

See coverm makedb --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint("coverm makedb"),
            ansi_term::Colour::Green
                .paint("Generate mapping database(s) from reference or genome FASTA files"),
            ansi_term::Colour::Purple.paint(
                "Example: Generate a short-read minimap2 database from a directory of\n\
                genomes, storing it in db_dir/"
            ),
        );
    }

    let mut app = Command::new("coverm")
        .version(crate_version!())
        .author(crate::AUTHOR_AND_EMAIL)
        .about("Mapping coverage analysis for metagenomics")
        .args(&[
            arg!(-v --verbose "Print extra debug logging information"),
            arg!(-q --quiet "Unless there is an error, do not print logging information"),
        ])
        .override_help(
            "
Mapping coverage analysis for metagenomics

Usage: coverm <subcommand> ...

Main subcommands:
\tcontig\tCalculate coverage of contigs
\tgenome\tCalculate coverage of genomes

Less used utility subcommands:
\tmake\tGenerate BAM files through alignment
\tmakedb\tGenerate mapping database(s) from reference FASTA files
\tfilter\tRemove (or only keep) alignments with insufficient identity
\tcluster\tDereplicate and cluster genomes
\tshell-completion
\t\tGenerate shell completion scripts

Other options:
\t-V, --version\tPrint version information

Ben J. Woodcroft <benjwoodcroft near gmail.com>
",
        )
        .arg_required_else_help(true)
        .subcommand(
            add_clap_verbosity_flags(Command::new("genome"))
                .about("Calculate coverage of genomes")
                .override_help(GENOME_HELP.as_str())
                .arg(
                    Arg::new("full-help")
                        .long("full-help")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("full-help-roff")
                        .long("full-help-roff")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("bam-files")
                        .short('b')
                        .long("bam-files")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .num_args(1..),
                )
                .arg(
                    Arg::new("sharded")
                        .long("sharded")
                        .required(false)
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("exclude-genomes-from-deshard")
                        .long("exclude-genomes-from-deshard")
                        .requires("sharded"),
                )
                .arg(
                    Arg::new("read1")
                        .short('1')
                        .long("read1")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .requires("read2")
                        .required_unless_present_any([
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("read2")
                        .short('2')
                        .long("read2")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .num_args(1..)
                        .requires("read1")
                        .required_unless_present_any([
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("coupled")
                        .short('c')
                        .long("coupled")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any([
                            "bam-files",
                            "read1",
                            "interleaved",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("interleaved")
                        .long("interleaved")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any([
                            "bam-files",
                            "read1",
                            "coupled",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("single")
                        .long("single")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any([
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("reference")
                        .short('r')
                        .long("reference")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("cache-unfiltered-bam-directory")
                        .long("cache-unfiltered-bam-directory")
                        .visible_alias("bam-file-cache-directory")
                        .conflicts_with("bam-files")
                        .conflicts_with("cache-unfiltered-bam-files"),
                )
                .arg(
                    Arg::new("cache-unfiltered-bam-files")
                        .long("cache-unfiltered-bam-files")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .conflicts_with("bam-files")
                        .conflicts_with("cache-unfiltered-bam-directory"),
                )
                .group(ArgGroup::new("cache-unfiltered-bam").args([
                    "cache-unfiltered-bam-directory",
                    "cache-unfiltered-bam-files",
                ]))
                .arg(
                    Arg::new("threads")
                        .short('t')
                        .long("threads")
                        .value_parser(clap::value_parser!(u16))
                        .default_value("1"),
                )
                .arg(
                    Arg::new("mapper")
                        .short('p')
                        .long("mapper")
                        .value_parser(MAPPING_SOFTWARE_LIST.iter().collect::<Vec<_>>())
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::new("minimap2-params")
                        .long("minimap2-params")
                        .alias("minimap2-parameters")
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::new("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index")
                        .requires("reference")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("bwa-params")
                        .long("bwa-params")
                        .alias("bwa-parameters")
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                .arg(
                    Arg::new("minibwa-params")
                        .long("minibwa-params")
                        .alias("minibwa-parameters")
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                .arg(
                    Arg::new("strobealign-params")
                        .long("strobealign-params")
                        .alias("strobealign-parameters")
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                .arg(
                    Arg::new("strobealign-use-index")
                        .long("strobealign-use-index")
                        .requires("reference")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("rammap-params")
                        .long("rammap-params")
                        .visible_alias("rammap-parameters")
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                // TODO: Relax this for autoconcatenation
                .arg(
                    Arg::new("discard-unmapped")
                        .long("discard-unmapped")
                        .requires("cache-unfiltered-bam")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("separator")
                        .short('s')
                        .long("separator")
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("genome-fasta-directory")
                        .conflicts_with("single-genome")
                        .required_unless_present_any([
                            "genome-fasta-files",
                            "genome-fasta-directory",
                            "genome-fasta-list",
                            "single-genome",
                            "genome-definition",
                            "full-help",
                            "full-help-roff",
                        ])
                        .value_parser(clap::value_parser!(char)),
                )
                .arg(
                    Arg::new("genome-fasta-files")
                        .short('f')
                        .long("genome-fasta-files")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .conflicts_with("separator")
                        .conflicts_with("genome-fasta-directory")
                        .conflicts_with("single-genome")
                        .required_unless_present_any([
                            "separator",
                            "genome-fasta-directory",
                            "genome-fasta-list",
                            "single-genome",
                            "genome-definition",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::new("genome-fasta-directory")
                        .short('d')
                        .long("genome-fasta-directory")
                        .conflicts_with("separator")
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("single-genome")
                        .required_unless_present_any([
                            "genome-fasta-files",
                            "genome-fasta-list",
                            "separator",
                            "single-genome",
                            "genome-definition",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::new("genome-fasta-list")
                        .long("genome-fasta-list")
                        .conflicts_with("separator")
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("genome-fasta-directory")
                        .conflicts_with("single-genome")
                        .required_unless_present_any([
                            "genome-fasta-files",
                            "genome-fasta-directory",
                            "separator",
                            "single-genome",
                            "genome-definition",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::new("genome-fasta-extension")
                        .short('x')
                        .long("genome-fasta-extension")
                        // Unsure why, but uncommenting causes test failure - clap
                        // bug?
                        //.requires("genome-fasta-directory")
                        .default_value("fna"),
                )
                .arg(
                    Arg::new("genome-definition")
                        .long("genome-definition")
                        .conflicts_with("separator")
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("genome-fasta-directory")
                        .conflicts_with("single-genome")
                        .required_unless_present_any([
                            "genome-fasta-files",
                            "genome-fasta-list",
                            "separator",
                            "single-genome",
                            "genome-fasta-directory",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::new("single-genome")
                        .long("single-genome")
                        .conflicts_with("separator")
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("genome-fasta-directory")
                        .conflicts_with("genome-definition")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("use-full-contig-names")
                        .long("use-full-contig-names")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("min-mapq")
                        .long("min-mapq")
                        .value_parser(clap::value_parser!(u8)),
                )
                .arg(
                    Arg::new("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .value_parser(clap::value_parser!(u32))
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::new("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .value_parser(clap::value_parser!(f32))
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::new("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .value_parser(clap::value_parser!(f32))
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::new("methods")
                        .short('m')
                        .long("method")
                        .long("methods")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .value_parser([
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
                            "anir",
                            "rpkm",
                            "tpm",
                        ])
                        .default_value("relative_abundance"),
                )
                .arg(
                    Arg::new("trim-min")
                        .long("trim-min")
                        .default_value("5")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("trim-max")
                        .long("trim-max")
                        .default_value("95")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("min-covered-fraction")
                        .long("min-covered-fraction")
                        .default_value("10")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("75")
                        .value_parser(clap::value_parser!(u64)),
                )
                .arg(
                    Arg::new("no-zeros")
                        .long("no-zeros")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("proper-pairs-only")
                        .long("proper-pairs-only")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("exclude-supplementary")
                        .long("exclude-supplementary")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("include-secondary")
                        .long("include-secondary")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(Arg::new("output-file").long("output-file").short('o'))
                .arg(
                    Arg::new("output-format")
                        .long("output-format")
                        .value_parser(["sparse", "dense"])
                        .default_value("dense"),
                )
                .arg(
                    Arg::new("gff")
                        .long("gff")
                        .visible_alias("genes-file")
                        .conflicts_with("dereplicate"),
                )
                .arg(
                    Arg::new("gff-feature-type")
                        .long("gff-feature-type")
                        .requires("gff"),
                )
                .arg(
                    Arg::new("dereplicate")
                        .long("dereplicate")
                        .conflicts_with("reference")
                        .conflicts_with("bam-files")
                        .conflicts_with("separator")
                        .conflicts_with("single-genome")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("dereplication-ani")
                        .long("dereplication-ani")
                        .value_parser(clap::value_parser!(f32))
                        .default_value(galah::DEFAULT_ANI),
                )
                .arg(
                    Arg::new("dereplication-prethreshold-ani")
                        .long("dereplication-prethreshold-ani")
                        .value_parser(clap::value_parser!(f32))
                        .default_value(galah::DEFAULT_PRETHRESHOLD_ANI),
                )
                .arg(
                    Arg::new("dereplication-quality-formula")
                        .long("dereplication-quality-formula")
                        .value_parser([
                            "completeness-4contamination",
                            "completeness-5contamination",
                            "Parks2020_reduced",
                            "dRep",
                        ])
                        .default_value(galah::DEFAULT_QUALITY_FORMULA),
                )
                .arg(
                    Arg::new("dereplication-precluster-method")
                        .long("dereplication-precluster-method")
                        .value_parser(galah::PRECLUSTER_METHODS)
                        .default_value(galah::DEFAULT_PRECLUSTER_METHOD),
                )
                .arg(
                    Arg::new("dereplication-cluster-method")
                        .long("dereplication-cluster-method")
                        .value_parser(galah::CLUSTER_METHODS)
                        .default_value(galah::DEFAULT_CLUSTER_METHOD),
                )
                .arg(
                    Arg::new("dereplication-aligned-fraction")
                        .long("dereplication-aligned-fraction")
                        .default_value(galah::DEFAULT_ALIGNED_FRACTION)
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("dereplication-small-genomes")
                        .long("dereplication-small-genomes")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("dereplication-cluster-contigs")
                        .long("dereplication-cluster-contigs")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("dereplication-small-contigs")
                        .long("dereplication-small-contigs")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("dereplication-large-contigs")
                        .long("dereplication-large-contigs")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("dereplication-fragment-length")
                        .long("dereplication-fragment-length")
                        .default_value(galah::DEFAULT_FRAGMENT_LENGTH)
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("dereplication-low-memory")
                        .long("dereplication-low-memory")
                        .action(clap::ArgAction::SetTrue)
                        .conflicts_with("dereplication-reference-genomes")
                        .conflicts_with("dereplication-reference-genomes-list"),
                )
                .arg(
                    Arg::new("dereplication-reference-genomes")
                        .long("dereplication-reference-genomes")
                        .num_args(1..)
                        .conflicts_with("dereplication-low-memory")
                        .conflicts_with("dereplication-reference-genomes-list"),
                )
                .arg(
                    Arg::new("dereplication-reference-genomes-list")
                        .long("dereplication-reference-genomes-list")
                        .conflicts_with("dereplication-low-memory")
                        .conflicts_with("dereplication-reference-genomes"),
                )
                .arg(
                    Arg::new("dereplication-output-cluster-definition")
                        .long("dereplication-output-cluster-definition"),
                )
                .arg(
                    Arg::new("dereplication-output-representative-fasta-directory")
                        .long("dereplication-output-representative-fasta-directory"),
                )
                .arg(
                    Arg::new("dereplication-output-representative-fasta-directory-copy")
                        .long("dereplication-output-representative-fasta-directory-copy"),
                )
                .arg(
                    Arg::new("dereplication-output-representative-list")
                        .long("dereplication-output-representative-list"),
                )
                .arg(
                    Arg::new("checkm-tab-table")
                        .long("checkm-tab-table")
                        .conflicts_with("reference")
                        .conflicts_with("bam-files")
                        .conflicts_with("separator")
                        .conflicts_with("single-genome"),
                )
                .arg(
                    Arg::new("checkm2-quality-report")
                        .long("checkm2-quality-report")
                        .conflicts_with("reference")
                        .conflicts_with("bam-files")
                        .conflicts_with("separator")
                        .conflicts_with("single-genome"),
                )
                .arg(
                    Arg::new("genome-info")
                        .long("genome-info")
                        .conflicts_with("reference")
                        .conflicts_with("bam-files")
                        .conflicts_with("separator")
                        .conflicts_with("single-genome")
                        .conflicts_with("checkm-tab-table"),
                )
                .arg(
                    Arg::new("min-completeness").long("min-completeness"), // Ideally we would use
                                                                           // requires_any(["checkm2-quality-report",
                                                                           // "checkm-tab-table"]) or similar, but it is not
                                                                           // immediately obvious how to implement this, so instead
                                                                           // we manually check at runtime.
                )
                .arg(Arg::new("max-contamination").long("max-contamination"))
                .arg(
                    Arg::new("run-checkm2")
                        .long("run-checkm2")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(Arg::new("checkm2-db-path").long("checkm2-db-path")),
        )
        .subcommand(
            add_clap_verbosity_flags(Command::new("contig"))
                .about("Calculate coverage of contigs")
                .override_help(CONTIG_HELP.as_str())
                .arg(
                    Arg::new("full-help")
                        .long("full-help")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("full-help-roff")
                        .long("full-help-roff")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("bam-files")
                        .short('b')
                        .long("bam-files")
                        .action(clap::ArgAction::Append)
                        .num_args(1..),
                )
                .arg(
                    Arg::new("sharded")
                        .long("sharded")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("read1")
                        .short('1')
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .requires("read2")
                        .required_unless_present_any([
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("read2")
                        .short('2')
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .requires("read1")
                        .required_unless_present_any([
                            "bam-files",
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("coupled")
                        .short('c')
                        .long("coupled")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any([
                            "bam-files",
                            "read1",
                            "interleaved",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("interleaved")
                        .long("interleaved")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any([
                            "bam-files",
                            "read1",
                            "coupled",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("single")
                        .long("single")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any([
                            "bam-files",
                            "read1",
                            "coupled",
                            "interleaved",
                            "full-help",
                            "full-help-roff",
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("reference")
                        .short('r')
                        .long("reference")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(["bam-files", "full-help", "full-help-roff"])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::new("cache-unfiltered-bam-directory")
                        .long("cache-unfiltered-bam-directory")
                        .visible_alias("bam-file-cache-directory")
                        .conflicts_with("bam-files")
                        .conflicts_with("cache-unfiltered-bam-files"),
                )
                .arg(
                    Arg::new("cache-unfiltered-bam-files")
                        .long("cache-unfiltered-bam-files")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .conflicts_with("bam-files")
                        .conflicts_with("cache-unfiltered-bam-directory"),
                )
                .group(ArgGroup::new("cache-unfiltered-bam").args([
                    "cache-unfiltered-bam-directory",
                    "cache-unfiltered-bam-files",
                ]))
                .arg(
                    Arg::new("threads")
                        .short('t')
                        .long("threads")
                        .default_value("1")
                        .value_parser(clap::value_parser!(u16)),
                )
                .arg(
                    Arg::new("mapper")
                        .short('p')
                        .long("mapper")
                        .value_parser(MAPPING_SOFTWARE_LIST.iter().collect::<Vec<_>>())
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::new("minimap2-params")
                        .long("minimap2-params")
                        .alias("minimap2-parameters")
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::new("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index")
                        .requires("reference")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("bwa-params")
                        .long("bwa-params")
                        .alias("bwa-parameters")
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                .arg(
                    Arg::new("minibwa-params")
                        .long("minibwa-params")
                        .alias("minibwa-parameters")
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                .arg(
                    Arg::new("strobealign-params")
                        .long("strobealign-params")
                        .alias("strobealign-parameters")
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                .arg(
                    Arg::new("strobealign-use-index")
                        .long("strobealign-use-index")
                        .requires("reference")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("rammap-params")
                        .long("rammap-params")
                        .visible_alias("rammap-parameters")
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                .arg(
                    Arg::new("discard-unmapped")
                        .long("discard-unmapped")
                        .requires("cache-unfiltered-bam")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("min-mapq")
                        .long("min-mapq")
                        .value_parser(clap::value_parser!(u8)),
                )
                .arg(
                    Arg::new("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .value_parser(clap::value_parser!(u32))
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::new("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .value_parser(clap::value_parser!(f32))
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::new("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .value_parser(clap::value_parser!(f32))
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::new("methods")
                        .short('m')
                        .long("method")
                        .long("methods")
                        .value_parser([
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
                            "anir",
                            "rpkm",
                            "tpm",
                            "strobealign-aemb",
                        ])
                        .default_value("mean")
                        .action(clap::ArgAction::Append)
                        .num_args(1..),
                )
                .arg(
                    Arg::new("min-covered-fraction")
                        .long("min-covered-fraction")
                        .default_value("0")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("75")
                        .value_parser(clap::value_parser!(u64)),
                )
                .arg(
                    Arg::new("trim-min")
                        .long("trim-min")
                        .default_value("5")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("trim-max")
                        .long("trim-max")
                        .default_value("95")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("no-zeros")
                        .long("no-zeros")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("proper-pairs-only")
                        .long("proper-pairs-only")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("exclude-supplementary")
                        .long("exclude-supplementary")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("include-secondary")
                        .long("include-secondary")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(Arg::new("output-file").long("output-file").short('o'))
                .arg(
                    Arg::new("output-format")
                        .long("output-format")
                        .value_parser(["sparse", "dense"])
                        .default_value("dense"),
                )
                .arg(Arg::new("gff").long("gff").visible_alias("genes-file"))
                .arg(
                    Arg::new("gff-feature-type")
                        .long("gff-feature-type")
                        .requires("gff"),
                ),
        )
        .subcommand(
            Command::new("filter") // Do not use add_clap_verbosity_flags since -v shouldn't be used here, specify manually below
                .about("Remove alignments with insufficient identity")
                .override_help(FILTER_HELP.as_str())
                .arg(
                    Arg::new("full-help")
                        .long("full-help")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("full-help-roff")
                        .long("full-help-roff")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("bam-files")
                        .short('b')
                        .long("bam-files")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(["full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::new("output-bam-files")
                        .short('o')
                        .long("output-bam-files")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(["full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::new("inverse")
                        .long("inverse")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("min-mapq")
                        .long("min-mapq")
                        .value_parser(clap::value_parser!(u8)),
                )
                .arg(
                    Arg::new("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .value_parser(clap::value_parser!(f32)),
                )
                .arg(
                    Arg::new("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .value_parser(clap::value_parser!(u32))
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::new("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .value_parser(clap::value_parser!(f32))
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::new("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .value_parser(clap::value_parser!(f32))
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::new("proper-pairs-only")
                        .long("proper-pairs-only")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("exclude-supplementary")
                        .long("exclude-supplementary")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("include-secondary")
                        .long("include-secondary")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("threads")
                        .long("threads")
                        .short('t')
                        .value_parser(clap::value_parser!(u16))
                        .default_value("1"),
                )
                .arg(
                    Arg::new("verbose")
                        // .short( 'v') // Do not use since could be confused with
                        // inverse (a la grep -v)
                        .long("verbose")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("quiet")
                        .short('q')
                        .long("quiet")
                        .action(clap::ArgAction::SetTrue),
                ),
        )
        .subcommand(
            add_clap_verbosity_flags(Command::new("make"))
                .about("Generate BAM files through mapping")
                .override_help(MAKE_HELP.as_str())
                .arg(
                    Arg::new("full-help")
                        .long("full-help")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("full-help-roff")
                        .long("full-help-roff")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("output-directory")
                        .short('o')
                        .long("output-directory")
                        .required_unless_present_any(["full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::new("read1")
                        .short('1')
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .requires("read2")
                        .required_unless_present_any([
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::new("read2")
                        .short('2')
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .requires("read1")
                        .required_unless_present_any([
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::new("coupled")
                        .short('c')
                        .long("coupled")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any([
                            "read1",
                            "interleaved",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::new("interleaved")
                        .long("interleaved")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any([
                            "read1",
                            "coupled",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::new("single")
                        .long("single")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any([
                            "read1",
                            "coupled",
                            "interleaved",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::new("reference")
                        .short('r')
                        .long("reference")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .required_unless_present_any(["full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::new("threads")
                        .short('t')
                        .long("threads")
                        .default_value("1")
                        .value_parser(clap::value_parser!(u16)),
                )
                .arg(
                    Arg::new("discard-unmapped")
                        .long("discard-unmapped")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("mapper")
                        .short('p')
                        .long("mapper")
                        .value_parser(MAPPING_SOFTWARE_LIST.iter().collect::<Vec<_>>())
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::new("minimap2-params")
                        .long("minimap2-params")
                        .alias("minimap2-parameters")
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::new("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index")
                        .requires("reference")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("bwa-params")
                        .long("bwa-params")
                        .alias("bwa-parameters")
                        .conflicts_with("minimap2-params")
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                .arg(
                    Arg::new("minibwa-params")
                        .long("minibwa-params")
                        .alias("minibwa-parameters")
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                .arg(
                    Arg::new("strobealign-params")
                        .long("strobealign-params")
                        .alias("strobealign-parameters")
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                .arg(
                    Arg::new("strobealign-use-index")
                        .long("strobealign-use-index")
                        .requires("reference")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("rammap-params")
                        .long("rammap-params")
                        .visible_alias("rammap-parameters")
                        .allow_hyphen_values(true)
                        .requires("reference"),
                ),
        )
        .subcommand(add_dereplication_arguments(
            add_clap_verbosity_flags(Command::new("makedb"))
                .about("Generate mapping database(s) from reference FASTA files")
                .override_help(MAKEDB_HELP.as_str())
                .arg(
                    Arg::new("full-help")
                        .long("full-help")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("full-help-roff")
                        .long("full-help-roff")
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("reference")
                        .short('r')
                        .long("reference")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("genome-fasta-directory")
                        .conflicts_with("genome-fasta-list")
                        .required_unless_present_any([
                            "genome-fasta-files",
                            "genome-fasta-directory",
                            "genome-fasta-list",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::new("genome-fasta-files")
                        .short('f')
                        .long("genome-fasta-files")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .conflicts_with("genome-fasta-directory")
                        .conflicts_with("genome-fasta-list")
                        .required_unless_present_any([
                            "reference",
                            "genome-fasta-directory",
                            "genome-fasta-list",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::new("genome-fasta-directory")
                        .short('d')
                        .long("genome-fasta-directory")
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("genome-fasta-list")
                        .required_unless_present_any([
                            "reference",
                            "genome-fasta-files",
                            "genome-fasta-list",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::new("genome-fasta-list")
                        .long("genome-fasta-list")
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("genome-fasta-directory")
                        .required_unless_present_any([
                            "reference",
                            "genome-fasta-files",
                            "genome-fasta-directory",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::new("genome-fasta-extension")
                        .short('x')
                        .long("genome-fasta-extension")
                        .default_value("fna"),
                )
                .arg(
                    Arg::new("output-directory")
                        .short('o')
                        .long("output-directory")
                        .required_unless_present_any(["full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::new("mapper")
                        .short('p')
                        .long("mapper")
                        .action(clap::ArgAction::Append)
                        .num_args(1..)
                        .value_parser(MAPPING_SOFTWARE_LIST.iter().collect::<Vec<_>>())
                        .default_value(DEFAULT_MAKEDB_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::new("threads")
                        .short('t')
                        .long("threads")
                        .default_value("1")
                        .value_parser(clap::value_parser!(u16)),
                )
                .arg(
                    Arg::new("minimap2-params")
                        .long("minimap2-params")
                        .alias("minimap2-parameters")
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::new("bwa-params")
                        .long("bwa-params")
                        .alias("bwa-parameters")
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::new("strobealign-params")
                        .long("strobealign-params")
                        .alias("strobealign-parameters")
                        .allow_hyphen_values(true),
                ),
        ))
        .subcommand(
            add_clap_verbosity_flags(Command::new("shell-completion"))
                .about("Generate a shell completion script for coverm")
                .arg(
                    Arg::new("output-file")
                        .short('o')
                        .long("output-file")
                        .required(true),
                )
                .arg(
                    Arg::new("shell")
                        .long("shell")
                        .required(true)
                        .value_parser(value_parser!(Shell)),
                ),
        );

    app = galah::cluster_argument_parsing::add_cluster_subcommand(app);
    app
}
