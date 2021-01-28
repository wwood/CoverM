use bird_tool_utils::clap_utils::{default_roff, monospace_roff};
use bird_tool_utils_man::prelude::{Author, Example, Flag, Manual, Opt, Section};
use clap::*;
use galah::cluster_argument_parsing::GalahClustererCommandDefinition;
use roff::bold;

const MAPPING_SOFTWARE_LIST: &[&str] = &[
    "bwa-mem",
    "minimap2-sr",
    "minimap2-ont",
    "minimap2-pb",
    "minimap2-no-preset",
];
const DEFAULT_MAPPING_SOFTWARE: &str = "minimap2-sr";

lazy_static! {
    pub static ref COVERM_CLUSTER_COMMAND_DEFINITION: GalahClustererCommandDefinition = {
        galah::cluster_argument_parsing::GalahClustererCommandDefinition {
            dereplication_ani_argument: "dereplication-ani".to_string(),
            dereplication_prethreshold_ani_argument: "dereplication-prethreshold-ani".to_string(),
            dereplication_quality_formula_argument: "dereplication-quality-formula".to_string(),
            dereplication_precluster_method_argument: "dereplication-precluster-method".to_string(),
            dereplication_aligned_fraction_argument: "dereplication-aligned-fraction".to_string(),
            dereplication_fraglen_argument: "dereplication-fragment-length".to_string(),
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
                default_roff("minimap2-sr"),
                bird_tool_utils::clap_utils::table_roff(&[
                    &["name", "description"],
                    &[
                        &monospace_roff("minimap2-sr"),
                        &format!("minimap2 with '{}' option", &monospace_roff("-x sr"))
                    ],
                    &[
                        &monospace_roff("bwa-mem"),
                        &format!("bwa mem using default parameters")
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
                        &monospace_roff("minimap2-no-preset"),
                        &format!("minimap2 with no '{}' option", &monospace_roff("-x"))
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
                "Extra parameters to provide to BWA. Note \
        that usage of this parameter has security \
        implications if untrusted input is specified. \
        [default: none]",
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
        identity e.g. 0.95 for 95%. {}",
                        default_roff("0.0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-aligned-percent")
                    .help(&format!(
                        "Exclude reads by percent aligned \
        bases e.g. 0.95 means 95% of the read's \
        bases must be aligned. {}",
                        default_roff("0.0")
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
                identity e.g. 0.95 for 95%. \
                Implies --proper-pairs-only. {}",
                        default_roff("0.0")
                    )),
            )
            .option(
                Opt::new("FLOAT")
                    .long("--min-read-aligned-percent-pair")
                    .help(&format!(
                        "Exclude reads by percent aligned \
                bases e.g. 0.95 means 95% of the read's \
                bases must be aligned. \
                Implies --proper-pairs-only. {}",
                        default_roff("0.0")
                    )),
            )
            .flag(
                Flag::new()
                    .long("--proper-pairs-only")
                    .help("Require reads to be mapped as proper pairs. [default: not set]"),
            )
            .flag(
                Flag::new()
                    .long("--exclude-supplementary")
                    .help("Exclude supplementary alignments. [default: not set]"),
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
        To fix, set the {} environment variable e.g. to set it to use the current directory: {}\n\n\

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
        .flag(Flag::new().short("-q").long("--quiet").help(
            "Unless there is an error, do not print \
    log messages. [default: not set]",
        ))
}

pub fn filter_full_help() -> Manual {
    let mut manual = Manual::new("coverm filter")
        .about(&format!(
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
                if they do not meet the thresholds. [default: not set]",
    ));
    manual = add_verbosity_flags(manual);
    manual = add_help_options(manual);

    manual = manual.example(
        Example::new()
            .text("Filter a BAM file by removing alignments shorter than 50bp")
            .command(
                "coverm filter --bam-files input.bam --output-bam filtered.bam \
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
                --min-read-percent-identity 0.95 --threads 16",
            ),
    );

    manual = manual.custom(faq_section());

    manual
}

pub fn make_full_help() -> Manual {
    let mut manual = Manual::new("coverm make")
        .about(&format!(
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
                        or BWA index stem (with {}). [required]",
                monospace_roff("--minimap2-reference-is-index"),
                monospace_roff("-p bwa-mem"),
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
                    genomes or metagenome assembly, or minimap2 \
                    index \
                    (with {}), \
                    or BWA index stem (with {}). \
                    If multiple references FASTA files are \
                    provided and {} is specified, \
                    then reads will be mapped to references \
                    separately as sharded BAMs. [required unless {} is specified]",
                monospace_roff("--minimap2-reference-is-index"),
                monospace_roff("-p bwa-mem"),
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

                    &[&monospace_roff("coverage_histogram"), "Histogram of coverage depths"],
                    &[&monospace_roff("covered_bases"), "Number of bases covered by 1 or more reads"],
                    &[&monospace_roff("variance"), "Variance of coverage depths"],
                    &[&monospace_roff("length"), "Length of each contig in base pairs"],
                    &[&monospace_roff("count"), "Number of reads aligned toq each contig. Note that a single read may be aligned to multiple contigs with supplementary alignments"],
                    &[&monospace_roff("metabat"), "(\"MetaBAT adjusted coverage\") Coverage as defined in Kang et al 2015 https://doi.org/10.7717/peerj.1165"],
                    &[&monospace_roff("reads_per_base"), "Number of reads aligned divided by the length of the contig"],
                    &[&monospace_roff("rpkm"), "Reads mapped per kilobase of contig, per million mapped reads"],
                    &[&monospace_roff("tpm"), "Transcripts Per Million as described in Li et al 2010 https://doi.org/10.1093/bioinformatics/btp692"],
                ]),
            )))
            .option(Opt::new("FRACTION").long("--min-covered-fraction").help(
                &format!("Genomes with less coverage than this \
                reported as having zero coverage. \
                {}", default_roff("0.10"))
            ))
            .option(Opt::new("INT").long("--contig-end-exclusion").help(
                &format!("Exclude bases at the ends of reference \
                sequences from calculation {}",
                default_roff("75"))
            ))
            .option(Opt::new("FRACTION").long("--trim-min").help(
                &format!("Remove this smallest fraction of positions \
                when calculating trimmed_mean {}",
                default_roff("0.05"))
            ))
            .option(Opt::new("FRACTION").long("--trim-max").help(
                &format!("Maximum fraction for trimmed_mean \
                calculations {}", default_roff("0.95"))
            )),
    );

    manual = manual.custom(
        Section::new("Output")
            .option(Opt::new("FILE").short("o").long("--output-file").help(
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
                    .long("--bam-file-cache-directory")
                    .help(
                        "Output BAM files generated during \
                alignment to this directory. The directory may or may not exist. \
                [default: not used]",
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
                "coverm contig --method metabat --bam-files my.bam \
                --bam-file-cache-directory saved_bam_files",
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

    return manual;
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
                    genomes or metagenome assembly, or minimap2 \
                    index \
                    (with {}), \
                    or BWA index stem (with {}). \
                    If multiple references FASTA files are \
                    provided and {} is specified, \
                    then reads will be mapped to references \
                    separately as sharded BAMs. {}: If genomic FASTA files are \
                    specified elsewhere (e.g. with {} or {}), then {} is not needed as a reference FASTA file can be derived \
                    by concatenating input genomes. However, while not necessary, {} can \
                    be specified if an alternate reference sequence set is desired.",
                    monospace_roff("--minimap2-reference-is-index"),
                    monospace_roff("-p bwa-mem"),
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
        );

    let mut derep_section = Section::new("DEREPLICATION / GENOME CLUSTERING").flag(
        Flag::new().long("--dereplicate").help(
            "Do genome dereplication via average nucleotide \
        identity (ANI) - choose a genome to represent \
        all within a small distance, using Dashing for \
        preclustering and FastANI for final ANI \
        calculation. When this flag is used, dereplication occurs \
        transparently through the Galah method (https://github.com/wwood/galah) [default: not set]",
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
                    &[&monospace_roff("covered_bases"), "Number of bases covered by 1 or more reads"],
                    &[&monospace_roff("variance"), "Variance of coverage depths"],
                    &[&monospace_roff("length"), "Length of each genome in base pairs"],
                    &[&monospace_roff("count"), "Number of reads aligned toq each genome. Note that a single read may be aligned to multiple genomes with supplementary alignments"],
                    &[&monospace_roff("reads_per_base"), "Number of reads aligned divided by the length of the genome"],
                    &[&monospace_roff("rpkm"), "Reads mapped per kilobase of genome, per million mapped reads"],
                    &[&monospace_roff("tpm"), "Transcripts Per Million as described in Li et al 2010 https://doi.org/10.1093/bioinformatics/btp692"],
                ])
            )))
            .option(Opt::new("FRACTION").long("--min-covered-fraction").help(
                &format!("Genomes with less coverage than this \
                reported as having zero coverage. \
                {}", default_roff("0.10"))
            ))
            .option(Opt::new("INT").long("--contig-end-exclusion").help(
                &format!("Exclude bases at the ends of reference \
                sequences from calculation {}",
                default_roff("75"))
            ))
            .option(Opt::new("FRACTION").long("--trim-min").help(
                &format!("Remove this smallest fraction of positions \
                when calculating trimmed_mean {}",
                default_roff("0.05"))
            ))
            .option(Opt::new("FRACTION").long("--trim-max").help(
                &format!("Maximum fraction for trimmed_mean \
                calculations {}", default_roff("0.95"))
            )),
    );

    manual = manual.custom(
        Section::new("Output")
            .option(Opt::new("FILE").short("o").long("--output-file").help(
                "Output coverage values to this file, or '-' for STDOUT. \
                [default: output to STDOUT]",
            ))
            .option(Opt::new("FORMAT").long("--output-format").help(
                &format!("Shape of output: 'sparse' for long format, \
            'dense' for species-by-site. {}", default_roff("dense"))
            ))
            .flag(Flag::new().long("--no-zeros").help(
                "Omit printing of genomes that have zero \
            coverage. [default: not set]",
            ))
            .option(
                Opt::new("DIRECTORY")
                    .long("--bam-file-cache-directory")
                    .help(
                        "Output BAM files generated during \
                    alignment to this directory. The directory may or may not exist. [default: not set]",
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

pub fn build_cli() -> App<'static, 'static> {
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
            ansi_term::Colour::Green.paint("coverm contig"),
            ansi_term::Colour::Green.paint("Calculate coverage of individual contigs"),
            ansi_term::Colour::Purple
                .paint("Example: Calculate mean coverage from reads and assembly:"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate MetaBAT adjusted coverage from a sorted BAM file, saving
the unfiltered BAM files in the saved_bam_files folder:"
            )
        )
        .to_string();
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
        )
        .to_string();
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
        )
        .to_string();
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
        )
        .to_string();
    }

    let mut app = App::new("coverm")
        .version(crate_version!())
        .author(crate::AUTHOR_AND_EMAIL)
        .about("Mapping coverage analysis for metagenomics")
        .args_from_usage(
            "-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'",
        )
        .help(
            "
Mapping coverage analysis for metagenomics

Usage: coverm <subcommand> ...

Main subcommands:
\tcontig\tCalculate coverage of contigs
\tgenome\tCalculate coverage of genomes

Less used utility subcommands:
\tmake\tGenerate BAM files through alignment
\tfilter\tRemove (or only keep) alignments with insufficient identity
\tcluster\tDereplicate and cluster genomes
\tshell-completion
\t\tGenerate shell completion scripts

Other options:
\t-V, --version\tPrint version information

Ben J. Woodcroft <benjwoodcroft near gmail.com>
",
        )
        .global_setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(
            SubCommand::with_name("genome")
                .about("Calculate coverage of genomes")
                .help(GENOME_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(Arg::with_name("full-help-roff").long("full-help-roff"))
                .arg(
                    Arg::with_name("bam-files")
                        .short("b")
                        .long("bam-files")
                        .multiple(true)
                        .takes_value(true),
                )
                .arg(Arg::with_name("sharded").long("sharded").required(false))
                .arg(
                    Arg::with_name("exclude-genomes-from-deshard")
                        .long("exclude-genomes-from-deshard")
                        .requires("sharded")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("read1")
                        .short("-1")
                        .long("read1")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read2")
                        .required_unless_one(&[
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
                    Arg::with_name("read2")
                        .short("-2")
                        .long("read2")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read1")
                        .required_unless_one(&[
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
                    Arg::with_name("coupled")
                        .short("-c")
                        .long("coupled")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
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
                    Arg::with_name("interleaved")
                        .long("interleaved")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
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
                    Arg::with_name("single")
                        .long("single")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
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
                    Arg::with_name("reference")
                        .short("-r")
                        .long("reference")
                        .takes_value(true)
                        .multiple(true)
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("bam-file-cache-directory")
                        .long("bam-file-cache-directory")
                        .takes_value(true)
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("threads")
                        .short("-t")
                        .long("threads")
                        .default_value("1")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("mapper")
                        .short("p")
                        .long("mapper")
                        .possible_values(MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("minimap2-params")
                        .long("minimap2-params")
                        .long("minimap2-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index")
                        .requires("reference"),
                )
                .arg(
                    Arg::with_name("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true)
                        .requires("reference"),
                ) // TODO: Relax this for autoconcatenation
                .arg(
                    Arg::with_name("discard-unmapped")
                        .long("discard-unmapped")
                        .requires("bam-file-cache-directory"),
                )
                .arg(
                    Arg::with_name("separator")
                        .short("s")
                        .long("separator")
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("genome-fasta-directory")
                        .conflicts_with("single-genome")
                        .required_unless_one(&[
                            "genome-fasta-files",
                            "genome-fasta-directory",
                            "genome-fasta-list",
                            "single-genome",
                            "genome-definition",
                            "full-help",
                            "full-help-roff",
                        ])
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("genome-fasta-files")
                        .short("f")
                        .long("genome-fasta-files")
                        .multiple(true)
                        .conflicts_with("separator")
                        .conflicts_with("genome-fasta-directory")
                        .conflicts_with("single-genome")
                        .required_unless_one(&[
                            "separator",
                            "genome-fasta-directory",
                            "genome-fasta-list",
                            "single-genome",
                            "genome-definition",
                            "full-help",
                            "full-help-roff",
                        ])
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("genome-fasta-directory")
                        .short("d")
                        .long("genome-fasta-directory")
                        .conflicts_with("separator")
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("single-genome")
                        .required_unless_one(&[
                            "genome-fasta-files",
                            "genome-fasta-list",
                            "separator",
                            "single-genome",
                            "genome-definition",
                            "full-help",
                            "full-help-roff",
                        ])
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("genome-fasta-list")
                        .long("genome-fasta-list")
                        .conflicts_with("separator")
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("genome-fasta-directory")
                        .conflicts_with("single-genome")
                        .required_unless_one(&[
                            "genome-fasta-files",
                            "genome-fasta-directory",
                            "separator",
                            "single-genome",
                            "genome-definition",
                            "full-help",
                            "full-help-roff",
                        ])
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("genome-fasta-extension")
                        .short("x")
                        .long("genome-fasta-extension")
                        // Unsure why, but uncommenting causes test failure - clap
                        // bug?
                        //.requires("genome-fasta-directory")
                        .default_value("fna")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("genome-definition")
                        .long("genome-definition")
                        .conflicts_with("separator")
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("genome-fasta-directory")
                        .conflicts_with("single-genome")
                        .required_unless_one(&[
                            "genome-fasta-files",
                            "genome-fasta-list",
                            "separator",
                            "single-genome",
                            "genome-fasta-directory",
                            "full-help",
                            "full-help-roff",
                        ])
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("single-genome")
                        .long("single-genome")
                        .conflicts_with("separator")
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("genome-fasta-directory")
                        .conflicts_with("genome-definition"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("methods")
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
                            "rpkm",
                            "tpm",
                        ])
                        .default_value("relative_abundance"),
                )
                .arg(
                    Arg::with_name("trim-min")
                        .long("trim-min")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("trim-max")
                        .long("trim-max")
                        .default_value("0.95"),
                )
                .arg(
                    Arg::with_name("min-covered-fraction")
                        .long("min-covered-fraction")
                        .default_value("0.10"),
                )
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("75"),
                )
                .arg(Arg::with_name("no-zeros").long("no-zeros"))
                .arg(Arg::with_name("proper-pairs-only").long("proper-pairs-only"))
                .arg(Arg::with_name("exclude-supplementary").long("exclude-supplementary"))
                .arg(
                    Arg::with_name("output-file")
                        .takes_value(true)
                        .long("output-file")
                        .short("o"),
                )
                .arg(
                    Arg::with_name("output-format")
                        .long("output-format")
                        .possible_values(&["sparse", "dense"])
                        .default_value("dense"),
                )
                .arg(
                    Arg::with_name("dereplicate")
                        .long("dereplicate")
                        .conflicts_with("reference")
                        .conflicts_with("bam-files")
                        .conflicts_with("separator")
                        .conflicts_with("single-genome"),
                )
                .arg(
                    Arg::with_name("dereplication-ani")
                        .long("dereplication-ani")
                        .takes_value(true)
                        .default_value(galah::DEFAULT_ANI),
                )
                .arg(
                    Arg::with_name("dereplication-prethreshold-ani")
                        .long("dereplication-prethreshold-ani")
                        .takes_value(true)
                        .default_value(galah::DEFAULT_PRETHRESHOLD_ANI),
                )
                .arg(
                    Arg::with_name("dereplication-quality-formula")
                        .long("dereplication-quality-formula")
                        .possible_values(&[
                            "completeness-4contamination",
                            "completeness-5contamination",
                            "Parks2020_reduced",
                            "dRep",
                        ])
                        .default_value(galah::DEFAULT_QUALITY_FORMULA)
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("dereplication-precluster-method")
                        .long("dereplication-precluster-method")
                        .possible_values(&["dashing", "finch"])
                        .default_value(galah::DEFAULT_PRECLUSTER_METHOD)
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("dereplication-aligned-fraction")
                        .long("dereplication-aligned-fraction")
                        .default_value(galah::DEFAULT_ALIGNED_FRACTION)
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("dereplication-fragment-length")
                        .long("dereplication-fragment-length")
                        .default_value(galah::DEFAULT_FRAGMENT_LENGTH)
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("dereplication-output-cluster-definition")
                        .long("dereplication-output-cluster-definition")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("dereplication-output-representative-fasta-directory")
                        .long("dereplication-output-representative-fasta-directory")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("dereplication-output-representative-fasta-directory-copy")
                        .long("dereplication-output-representative-fasta-directory-copy")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("dereplication-output-representative-list")
                        .long("dereplication-output-representative-list")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("checkm-tab-table")
                        .long("checkm-tab-table")
                        .conflicts_with("reference")
                        .conflicts_with("bam-file")
                        .conflicts_with("separator")
                        .conflicts_with("single-genome")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("genome-info")
                        .long("genome-info")
                        .conflicts_with("reference")
                        .conflicts_with("bam-file")
                        .conflicts_with("separator")
                        .conflicts_with("single-genome")
                        .conflicts_with("checkm-tab-table")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-completeness")
                        .long("min-completeness")
                        .requires("checkm-tab-table")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("max-contamination")
                        .long("max-contamination")
                        .requires("checkm-tab-table")
                        .takes_value(true),
                )
                .arg(Arg::with_name("verbose").short("v").long("verbose"))
                .arg(Arg::with_name("quiet").short("q").long("quiet")),
        )
        .subcommand(
            SubCommand::with_name("contig")
                .about("Calculate coverage of contigs")
                .help(CONTIG_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(Arg::with_name("full-help-roff").long("full-help-roff"))
                .arg(
                    Arg::with_name("bam-files")
                        .short("b")
                        .long("bam-files")
                        .multiple(true)
                        .takes_value(true),
                )
                .arg(Arg::with_name("sharded").long("sharded").required(false))
                .arg(
                    Arg::with_name("read1")
                        .short("-1")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read2")
                        .required_unless_one(&[
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
                    Arg::with_name("read2")
                        .short("-2")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read1")
                        .required_unless_one(&[
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
                    Arg::with_name("coupled")
                        .short("-c")
                        .long("coupled")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
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
                    Arg::with_name("interleaved")
                        .long("interleaved")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
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
                    Arg::with_name("single")
                        .long("single")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
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
                    Arg::with_name("reference")
                        .short("-r")
                        .long("reference")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless_one(&["bam-files", "full-help", "full-help-roff"])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("bam-file-cache-directory")
                        .long("bam-file-cache-directory")
                        .takes_value(true)
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("threads")
                        .short("-t")
                        .long("threads")
                        .default_value("1")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("mapper")
                        .short("p")
                        .long("mapper")
                        .possible_values(MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("minimap2-params")
                        .long("minimap2-params")
                        .long("minimap2-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index")
                        .requires("reference"),
                )
                .arg(
                    Arg::with_name("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                .arg(
                    Arg::with_name("discard-unmapped")
                        .long("discard-unmapped")
                        .requires("bam-file-cache-directory"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("methods")
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
                            "rpkm",
                            "tpm",
                        ])
                        .default_value("mean"),
                )
                .arg(
                    Arg::with_name("min-covered-fraction")
                        .long("min-covered-fraction")
                        .default_value("0.0"),
                )
                .arg(
                    Arg::with_name("contig-end-exclusion")
                        .long("contig-end-exclusion")
                        .default_value("75"),
                )
                .arg(
                    Arg::with_name("trim-min")
                        .long("trim-min")
                        .default_value("0.05"),
                )
                .arg(
                    Arg::with_name("trim-max")
                        .long("trim-max")
                        .default_value("0.95"),
                )
                .arg(Arg::with_name("no-zeros").long("no-zeros"))
                .arg(Arg::with_name("proper-pairs-only").long("proper-pairs-only"))
                .arg(Arg::with_name("exclude-supplementary").long("exclude-supplementary"))
                .arg(
                    Arg::with_name("output-file")
                        .takes_value(true)
                        .long("output-file")
                        .short("o"),
                )
                .arg(
                    Arg::with_name("output-format")
                        .long("output-format")
                        .possible_values(&["sparse", "dense"])
                        .default_value("dense"),
                )
                .arg(Arg::with_name("verbose").short("v").long("verbose"))
                .arg(Arg::with_name("quiet").short("q").long("quiet")),
        )
        .subcommand(
            SubCommand::with_name("filter")
                .about("Remove alignments with insufficient identity")
                .help(FILTER_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(Arg::with_name("full-help-roff").long("full-help-roff"))
                .arg(
                    Arg::with_name("bam-files")
                        .short("b")
                        .long("bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&["full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::with_name("output-bam-files")
                        .short("o")
                        .long("output-bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&["full-help", "full-help-roff"]),
                )
                .arg(Arg::with_name("inverse").long("inverse"))
                .arg(
                    Arg::with_name("min-read-aligned-length")
                        .long("min-read-aligned-length")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity")
                        .long("min-read-percent-identity")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent")
                        .long("min-read-aligned-percent")
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("min-read-aligned-length-pair")
                        .long("min-read-aligned-length-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-percent-identity-pair")
                        .long("min-read-percent-identity-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(
                    Arg::with_name("min-read-aligned-percent-pair")
                        .long("min-read-aligned-percent-pair")
                        .takes_value(true)
                        .requires("proper-pairs-only"),
                )
                .arg(Arg::with_name("proper-pairs-only").long("proper-pairs-only"))
                .arg(Arg::with_name("exclude-supplementary").long("exclude-supplementary"))
                .arg(
                    Arg::with_name("threads")
                        .long("threads")
                        .short("t")
                        .default_value("1"),
                )
                .arg(
                    Arg::with_name("verbose")
                        // .short("v") // Do not use since could be confused with
                        // inverse (a la grep -v)
                        .long("verbose"),
                )
                .arg(Arg::with_name("quiet").short("q").long("quiet")),
        )
        .subcommand(
            SubCommand::with_name("make")
                .about("Generate BAM files through mapping")
                .help(MAKE_HELP.as_str())
                .arg(Arg::with_name("full-help").long("full-help"))
                .arg(Arg::with_name("full-help-roff").long("full-help-roff"))
                .arg(
                    Arg::with_name("output-directory")
                        .short("-o")
                        .long("output-directory")
                        .takes_value(true)
                        .required_unless_one(&["full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::with_name("read1")
                        .short("-1")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read2")
                        .required_unless_one(&[
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::with_name("read2")
                        .short("-2")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read1")
                        .required_unless_one(&[
                            "coupled",
                            "interleaved",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::with_name("coupled")
                        .short("-c")
                        .long("coupled")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "read1",
                            "interleaved",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::with_name("interleaved")
                        .long("interleaved")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "read1",
                            "coupled",
                            "single",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::with_name("single")
                        .long("single")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&[
                            "read1",
                            "coupled",
                            "interleaved",
                            "full-help",
                            "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::with_name("reference")
                        .short("-r")
                        .long("reference")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&["full-help", "full-help-roff"]),
                )
                .arg(
                    Arg::with_name("threads")
                        .short("-t")
                        .long("threads")
                        .default_value("1")
                        .takes_value(true),
                )
                .arg(Arg::with_name("discard-unmapped").long("discard-unmapped"))
                .arg(
                    Arg::with_name("mapper")
                        .short("p")
                        .long("mapper")
                        .possible_values(MAPPING_SOFTWARE_LIST)
                        .default_value(DEFAULT_MAPPING_SOFTWARE),
                )
                .arg(
                    Arg::with_name("minimap2-params")
                        .long("minimap2-params")
                        .long("minimap2-parameters")
                        .takes_value(true)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("minimap2-reference-is-index")
                        .long("minimap2-reference-is-index")
                        .requires("reference"),
                )
                .arg(
                    Arg::with_name("bwa-params")
                        .long("bwa-params")
                        .long("bwa-parameters")
                        .conflicts_with("minimap2-params")
                        .takes_value(true)
                        .allow_hyphen_values(true)
                        .requires("reference"),
                )
                .arg(
                    Arg::with_name("verbose")
                        // .short("v") // Do not use since could be confused with
                        // inverse (a la grep -v)
                        .long("verbose"),
                )
                .arg(Arg::with_name("quiet").short("q").long("quiet")),
        )
        .subcommand(
            SubCommand::with_name("shell-completion")
                .about("Generate a shell completion script for coverm")
                .arg(
                    Arg::with_name("output-file")
                        .short("-o")
                        .long("output-file")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::with_name("shell")
                        .long("shell")
                        .takes_value(true)
                        .required(true)
                        .possible_values(&Shell::variants()),
                )
                .arg(
                    Arg::with_name("verbose")
                        // .short("v") // Do not use since could be confused with
                        // inverse (a la grep -v)
                        .long("verbose"),
                )
                .arg(Arg::with_name("quiet").short("q").long("quiet")),
        );

    app = galah::cluster_argument_parsing::add_cluster_subcommand(app);
    return app;
}
