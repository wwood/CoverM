use clap::*;
use man::prelude::{Author, Flag, Manual, Opt};

const MAPPING_SOFTWARE_LIST: &[&str] = &[
    "bwa-mem",
    "minimap2-sr",
    "minimap2-ont",
    "minimap2-pb",
    "minimap2-no-preset",
];
const DEFAULT_MAPPING_SOFTWARE: &str = "minimap2-sr";

pub fn filter_full_help() -> &'static str {
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
   --exclude-supplementary               Exclude supplementary alignments

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

fn add_mapping_options(manual: Manual) -> Manual {
    manual.option(Opt::new("NAME").short("-p").long("--mapper").help(
        "Underlying mapping software used \
                    (\"minimap2-sr\", \"bwa-mem\", \"minimap2-ont\", \
                    \"minimap2-pb\", or \"minimap2-no-preset\"). \
                    minimap2 -sr, -ont, -pb, -no-preset specify \
                    '-x' preset of minimap2 to be used \
                    (with map-ont, map-pb for -ont, -pb). \
                    [default: \"minimap2-sr\"]",
    ))
}

fn add_thresholding_options(manual: Manual) -> Manual {
    manual
        .option(Opt::new("INT").long("--min-read-aligned-length").help(
            "Exclude reads with smaller numbers of \
        aligned bases [default: 0]",
        ))
        .option(Opt::new("FLOAT").long("--min-read-percent-identity").help(
            "Exclude reads by overall percent \
        identity e.g. 0.95 for 95%. [default 0.0]",
        ))
        .option(Opt::new("FLOAT").long("--min-read-aligned-percent").help(
            "Exclude reads by percent aligned \
        bases e.g. 0.95 means 95% of the read's \
        bases must be aligned. [default 0.0]",
        ))
        .option(Opt::new("INT").long("--min-read-aligned-length-pair").help(
            "Exclude pairs with smaller numbers of \
        aligned bases. \
        Implies --proper-pairs-only. [default: 0]",
        ))
        .option(
            Opt::new("FLOAT")
                .long("--min-read-percent-identity-pair")
                .help(
                    "Exclude pairs by overall percent \
                identity e.g. 0.95 for 95%. \
                Implies --proper-pairs-only. [default 0.0]",
                ),
        )
        .option(
            Opt::new("FLOAT")
                .long("--min-read-aligned-percent-pair")
                .help(
                    "Exclude reads by percent aligned \
                bases e.g. 0.95 means 95% of the read's \
                bases must be aligned. \
                Implies --proper-pairs-only. [default 0.0]",
                ),
        )
        .flag(
            Flag::new()
                .long("--proper-pairs-only")
                .help("Require reads to be mapped as proper pairs"),
        )
        .flag(
            Flag::new()
                .long("--exclude-supplementary")
                .help("Exclude supplementary alignments"),
        )
}

fn add_read_params(manual: Manual) -> Manual {
    manual
        .option(
            Opt::new("PATH ..")
                .short("-1")
                .help("Forward FASTA/Q file(s) for mapping"),
        )
        .option(
            Opt::new("PATH ..")
                .short("-2")
                .help("Reverse FASTA/Q file(s) for mapping"),
        )
        .option(Opt::new("PATH ..").short("-c").long("--coupled").help(
            "One or more pairs of forward and reverse \
        FASTA/Q files for mapping in order \
        <sample1_R1.fq.gz> <sample1_R2.fq.gz> \
        <sample2_R1.fq.gz> <sample2_R2.fq.gz> ..",
        ))
        .option(
            Opt::new("PATH ..")
                .long("--interleaved")
                .help("Interleaved FASTA/Q files(s) for mapping."),
        )
        .option(
            Opt::new("PATH ..")
                .long("--single")
                .help("Unpaired FASTA/Q files(s) for mapping."),
        )
        .option(Opt::new("PARAMS").long("--minimap2-params").help(
            "Extra parameters to provide to minimap2, \
        both indexing command (if used) and for \
        mapping. Note that usage of this parameter \
        has security implications if untrusted input \
        is specified. '-a' is always specified. \
        [default \"\"]",
        ))
        .flag(
            Flag::new()
                .long("--minimap2-reference-is-index")
                .help("Treat reference as a minimap2 database, not as a FASTA file."),
        )
        .option(Opt::new("PARAMS").long("--bwa-params").help(
            "Extra parameters to provide to BWA. Note \
        that usage of this parameter has security \
        implications if untrusted input is specified. \
        [default \"\"]",
        ))
}

fn add_help_options(manual: Manual) -> Manual {
    manual
        .flag(
            Flag::new()
                .short("-h")
                .long("--help")
                .help("Output a short usage message."),
        )
        .flag(
            Flag::new()
                .long("--full-help")
                .help("Output a full help message and display in 'man'."),
        )
        .flag(Flag::new().long("--full-help-roff").help(
            "Output a full help message in raw ROFF format for \
        conversion to other formats.",
        ))
}

fn add_sharding_option(manual: Manual) -> Manual {
    manual.flag(Flag::new().long("--sharded").help(
        "If -b/--bam-files was used: \
        Input BAM files are read-sorted alignments \
        of a set of reads mapped to multiple \
        reference contig sets. Choose the best \
        hit for each read pair. Otherwise if mapping was carried out: \
        Map reads to each reference, choosing the \
        best hit for each pair.",
    ))
}

fn add_verbosity_flags(manual: Manual) -> Manual {
    manual
        .flag(
            Flag::new()
                .short("-v")
                .long("--verbose")
                .help("Print extra debugging information"),
        )
        .flag(Flag::new().short("-q").long("--quiet").help(
            "Unless there is an error, do not print \
    log messages",
        ))
}

pub fn make_full_help() -> Manual {
    let mut manual = Manual::new("coverm make")
        .about("Generate BAM files through mapping")
        .author(Author::new("Ben J Woodcroft").email("benjwoodcroft near gmail.com"))
        .option(Opt::new("DIR").short("-o").long("--output-directory").help(
            "Where generated BAM files will go. The directory will be created if it does not exist.",
        ));

    manual = manual
        .option(Opt::new("PATH").short("-r").long("--reference").help(
            "FASTA file of contigs e.g. concatenated \
                    genomes or metagenome assembly, or minimap2 \
                    index \
                    (with --minimap2-reference-is-index), \
                    or BWA index stem (with -p bwa-mem). \
                    If multiple references FASTA files are \
                    provided and --sharded is specified, \
                    then reads will be mapped to references \
                    separately as sharded BAMs.",
        ))
        .option(
            Opt::new("INT")
                .short("-t")
                .long("--threads")
                .help("Number of threads for mapping."),
        );

    manual = add_read_params(manual);
    manual = add_mapping_options(manual);
    manual = add_help_options(manual);

    manual = manual.flag(
        Flag::new()
            .long("--discard-unmapped")
            .help("Exclude unmapped reads from cached BAM files."),
    );

    manual = add_verbosity_flags(manual);
    manual
}

pub fn contig_full_help() -> Manual {
    let mut manual = Manual::new("coverm contig")
        .about("Calculate read coverage per-contig")
        .author(Author::new("Ben J Woodcroft").email("benjwoodcroft near gmail.com"))
        .option(Opt::new("PATH").short("-b").long("--bam-files").help(
            "Path to BAM file(s). These must be \
                reference sorted (e.g. with samtools sort) \
                unless --sharded is specified, in which \
                case they must be read name sorted (e.g. \
                with samtools sort -n).",
        ));
    manual = add_help_options(manual);
    manual = add_mapping_options(manual);
    manual = add_thresholding_options(manual);
    manual = add_sharding_option(manual);
    manual = manual
        .option(Opt::new("METHOD").short("-m").long("--methods").help(
            "Method(s) for calculating coverage. \
            One or more (space separated) of: \
            mean (default), \
            trimmed_mean, \
            coverage_histogram, \
            covered_fraction, \
            covered_bases, \
            variance, \
            length, \
            count, \
            metabat (\"MetaBAT adjusted coverage\"), \
            reads_per_base, \
            rpkm. \
        A more thorough description of the different \
        methods is available at \
        https://github.com/wwood/CoverM",
        ))
        .option(Opt::new("FORMAT").long("--output-format").help(
            "Shape of output: 'sparse' for long format, \
        'dense' for species-by-site. \
        [default: dense]",
        ))
        .option(Opt::new("FRACTION").long("--min-covered-fraction").help(
            "Genomes with less coverage than this \
        reported as having zero coverage. \
        [default: 0.10]",
        ))
        .option(Opt::new("INT").long("--contig-end-exclusion").help(
            "Exclude bases at the ends of reference \
        sequences from calculation [default: 75]",
        ))
        .option(Opt::new("FRACTION").long("--trim-min").help(
            "Remove this smallest fraction of positions \
        when calculating trimmed_mean \
        [default: 0.05]",
        ))
        .option(Opt::new("FRACTION").long("--trim-max").help(
            "Maximum fraction for trimmed_mean \
        calculations [default: 0.95]",
        ))
        .flag(Flag::new().long("--no-zeros").help(
            "Omit printing of genomes that have zero \
        coverage",
        ))
        .option(
            Opt::new("DIRECTORY")
                .long("--bam-file-cache-directory")
                .help(
                    "Output BAM files generated during \
                alignment to this directory. The directory may or may not exist",
                ),
        )
        .flag(
            Flag::new()
                .long("--discard-unmapped")
                .help("Exclude unmapped reads from cached BAM files."),
        );

    manual = add_verbosity_flags(manual);
    manual = add_help_options(manual);

    return manual;
}

pub fn genome_full_help() -> Manual {
    let mut manual = Manual::new("coverm genome")
            .about("Calculate read coverage per-genome")
            .author(Author::new("Ben J Woodcroft").email("benjwoodcroft near gmail.com"))
            .option(
                Opt::new("CHARACTER")
                    .short("-s")
                    .long("--separator")
                    .help("This character separates genome names from contig names in the reference file")
            )
            .option(
                Opt::new("PATH ..")
                    .short("-f")
                    .long("--genome-fasta-files")
                    .help("Path(s) to FASTA files of each genome e.g. 'pathA/genome1.fna pathB/genome2.fa'")
            )
            .option(
                Opt::new("PATH")
                    .short("-d")
                    .long("--genome-fasta-directory")
                    .help("Directory containing FASTA files of each genome")
            )
            .option(
                Opt::new("EXT")
                    .short("-x")
                    .long("--genome-fasta-extension")
                    .help("File extension of genomes in the directory \
                        specified with -d/--genome-fasta-directory [default \"fna\"]")
            )
            .option(
                Opt::new("PATH")
                    .short("-d")
                    .long("--genome-fasta-list")
                    .help("File containing FASTA file paths, one per line")
            )
            .option(
                Opt::new("FILE")
                    .long("--genome-definition")
                    .help("File containing list of \
                    genome_name<tab>contig lines to define the genome of each contig")
            )
            .flag(
                Flag::new()
                    .long("--single-genome")
                    .help("All contigs are from the same genome")
            )

            .option(
                Opt::new("PATH")
                    .short("-b")
                    .long("--bam-files")
                    .help("Path to BAM file(s). These must be \
                        reference sorted (e.g. with samtools sort) \
                        unless --sharded is specified, in which \
                        case they must be read name sorted (e.g. \
                        with samtools sort -n).")
            );

    manual = add_mapping_options(manual);
    manual = add_help_options(manual);

    manual = manual
        .option(Opt::new("PATH").short("-r").long("--reference").help(
            "FASTA file of contigs e.g. concatenated \
                    genomes or metagenome assembly, or minimap2 \
                    index \
                    (with --minimap2-reference-is-index), \
                    or BWA index stem (with -p bwa-mem). \
                    If multiple references FASTA files are \
                    provided and --sharded is specified, \
                    then reads will be mapped to references \
                    separately as sharded BAMs.",
        ))
        .option(
            Opt::new("INT")
                .short("-t")
                .long("--threads")
                .help("Number of threads for mapping, sorting and reading."),
        );

    manual = add_read_params(manual);
    manual = add_sharding_option(manual);

    manual = manual.flag(Flag::new().long("--exclude-genomes-from-deshard").help(
        "Ignore genomes whose name appears in this newline-separated \
                file when combining shards.",
    ));
    manual = add_thresholding_options(manual);
    manual = manual
        .flag(Flag::new().long("--dereplicate").help(
            "Do genome dereplication via average nucleotide \
            identity (ANI) - choose a genome to represent \
            all within a small distance, using Dashing for \
            preclustering and FastANI for final ANI \
            calculation.",
        ))
        .option(
            Opt::new("FLOAT")
                .long("--dereplication-ani")
                .help("Overall ANI level to dereplicate at with FastANI."),
        )
        .option(Opt::new("PATH").long("--checkm-tab-table").help(
            "CheckM tab table for defining genome quality, \
            which is in turn used during clustering. \
            Genomes are scored as \
            completeness-4*contamination.",
        ))
        .option(Opt::new("PATH").long("--genome-info").help(
            "dRep style genome info tabl for defining \
            quality. Used like --checkm-tab-table.",
        ))
        .option(Opt::new("FLOAT").long("--min-completeness").help(
            "Ignore genomes with less completeness than \
            this percentage.",
        ))
        .option(Opt::new("FLOAT").long("--max-contamination").help(
            "Ignore genomes with more contamination than \
            this percentage.",
        ))
        .option(
            Opt::new("PATH")
                .long("--output-dereplication-clusters")
                .help(
                    "Output clustered genome information to this \
                    file as 'representative<TAB>member'",
                ),
        )
        .option(
            Opt::new("FLOAT")
                .long("--dereplication-prethreshold-ani")
                .help(
                    "Require at least this dashing-derived ANI \
                    for preclustering and to avoid FastANI on \
                    distant lineages within preclusters.",
                ),
        )
        .option(
            Opt::new("NAME")
                .long("--dereplication-quality-formula")
                .help(
                    "Scoring function for genome quality. See \
                    `coverm dereplicate -h`.",
                ),
        )
        .option(
            Opt::new("NAME")
                .long("--dereplication-precluster-method")
                .help(
                    "method of calculating rough ANI for \
                    dereplication. 'dashing' for HyperLogLog, \
                    'finch' for finch MinHash.",
                ),
        )
        .option(Opt::new("METHOD").short("-m").long("--methods").help(
            "Method(s) for calculating coverage. \
            One or more (space separated) of: \
                relative_abundance (default), \
                mean, \
                trimmed_mean, \
                coverage_histogram, \
                covered_fraction, \
                covered_bases, \
                variance, \
                length, \
                count, \
                reads_per_base, \
                rpkm. \
            A more thorough description of the different \
            methods is available at \
            https://github.com/wwood/CoverM",
        ))
        .option(Opt::new("FORMAT").long("--output-format").help(
            "Shape of output: 'sparse' for long format, \
            'dense' for species-by-site. \
            [default: dense]",
        ))
        .option(Opt::new("FRACTION").long("--min-covered-fraction").help(
            "Genomes with less coverage than this \
            reported as having zero coverage. \
            [default: 0.10]",
        ))
        .option(Opt::new("INT").long("--contig-end-exclusion").help(
            "Exclude bases at the ends of reference \
            sequences from calculation [default: 75]",
        ))
        .option(Opt::new("FRACTION").long("--trim-min").help(
            "Remove this smallest fraction of positions \
            when calculating trimmed_mean \
            [default: 0.05]",
        ))
        .option(Opt::new("FRACTION").long("--trim-max").help(
            "Maximum fraction for trimmed_mean \
            calculations [default: 0.95]",
        ))
        .flag(Flag::new().long("--no-zeros").help(
            "Omit printing of genomes that have zero \
            coverage",
        ))
        .option(
            Opt::new("DIRECTORY")
                .long("--bam-file-cache-directory")
                .help(
                    "Output BAM files generated during \
                    alignment to this directory. The directory may or may not exist",
                ),
        )
        .flag(
            Flag::new()
                .long("--discard-unmapped")
                .help("Exclude unmapped reads from cached BAM files."),
        );

    manual = add_verbosity_flags(manual);
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
    --genome-fasta-files genome1.fna genome2.fna

{}

  coverm genome --bam-files my.bam --genome-fasta-directory genomes_directory/

{}

  coverm genome --genome-fasta-directory genomes/ --dereplicate
    --single single_reads.fq.gz

See coverm genome --full-help for further options and further detail.
",
            ansi_term::Colour::Green.paint("coverm genome"),
            ansi_term::Colour::Green.paint("Calculate coverage of individual genomes"),
            ansi_term::Colour::Purple.paint("Example: Map paired reads to 2 genomes:"),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate coverage of genomes defined as .fna files in
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
        .author("Ben J. Woodcroft <benjwoodcroft near gmail.com>")
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
                    Arg::with_name("output-format")
                        .long("output-format")
                        .possible_values(&["sparse", "dense"])
                        .default_value("dense"),
                )

                .arg(Arg::with_name("dereplicate")
                    .long("dereplicate")
                    .conflicts_with("reference")
                    .conflicts_with("bam-files")
                    .conflicts_with("separator")
                    .conflicts_with("single-genome"))
                .arg(Arg::with_name("dereplication-ani")
                    .long("dereplication-ani")
                    .takes_value(true)
                    .default_value("99")
                )
                .arg(Arg::with_name("dereplication-prethreshold-ani")
                    .long("dereplication-prethreshold-ani")
                    .takes_value(true)
                    .default_value("95"))
                .arg(Arg::with_name("dereplication-quality-formula")
                    .long("dereplication-quality-formula")
                    .possible_values(&[
                        "completeness-4contamination",
                        "completeness-5contamination",
                        "Parks2020_reduced",
                        "dRep"])
                    .default_value("Parks2020_reduced")
                    .takes_value(true))
                .arg(Arg::with_name("dereplication-precluster-method")
                    .long("dereplication-precluster-method")
                    .help("method of calculating rough ANI. 'dashing' for HyperLogLog, 'finch' for finch MinHash")
                    .possible_values(&["dashing","finch"])
                    .default_value("dashing")
                    .takes_value(true))
                .arg(Arg::with_name("output-dereplication-clusters")
                    .long("output-dereplication-clusters")
                    .takes_value(true))
                .arg(Arg::with_name("checkm-tab-table")
                    .long("checkm-tab-table")
                    .conflicts_with("reference")
                    .conflicts_with("bam-file")
                    .conflicts_with("separator")
                    .conflicts_with("single-genome")
                    .takes_value(true))
                .arg(Arg::with_name("genome-info")
                    .long("genome-info")
                    .conflicts_with("reference")
                    .conflicts_with("bam-file")
                    .conflicts_with("separator")
                    .conflicts_with("single-genome")
                    .conflicts_with("checkm-tab-table")
                    .takes_value(true))
                .arg(Arg::with_name("min-completeness")
                    .long("min-completeness")
                    .requires("checkm-tab-table")
                    .takes_value(true))
                .arg(Arg::with_name("max-contamination")
                    .long("max-contamination")
                    .requires("checkm-tab-table")
                    .takes_value(true))


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
                .arg(
                    Arg::with_name("bam-files")
                        .short("b")
                        .long("bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&["full-help"]),
                )
                .arg(
                    Arg::with_name("output-bam-files")
                        .short("o")
                        .long("output-bam-files")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&["full-help"]),
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
                        .required_unless_one(&[
                            "full-help",
                            "full-help-roff",
                        ])
                )
                .arg(
                    Arg::with_name("read1")
                        .short("-1")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read2")
                        .required_unless_one(&["coupled", "interleaved", "single",
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
                        .required_unless_one(&["coupled", "interleaved", "single",
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
                        .required_unless_one(&["read1", "interleaved", "single",
                        "full-help",
                        "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::with_name("interleaved")
                        .long("interleaved")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&["read1", "coupled", "single",
                        "full-help",
                        "full-help-roff",
                        ]),
                )
                .arg(
                    Arg::with_name("single")
                        .long("single")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&["read1", "coupled", "interleaved",
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
                        .required_unless_one(&[
                            "full-help",
                            "full-help-roff",
                        ]),
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
                        .required(true)
                )
                .arg(
                    Arg::with_name("shell")
                        .long("shell")
                        .takes_value(true)
                        .required(true)
                        .possible_values(&Shell::variants())
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
