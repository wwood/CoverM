use clap::*;

const MAPPING_SOFTWARE_LIST: &[&str] = &["bwa-mem", "minimap2-sr", "minimap2-ont", "minimap2-pb","minimap2-no-preset"];
const DEFAULT_MAPPING_SOFTWARE: &str = "minimap2-sr";

const MAPPER_HELP: &'static str = 
"   -p, --mapper <NAME>                   Underlying mapping software used
                                         (\"minimap2-sr\", \"bwa-mem\", \"minimap2-ont\",
                                         \"minimap2-pb\", or \"minimap2-no-preset\").
                                         minimap2 -sr, -ont, -pb, -no-preset specify
                                         '-x' preset of minimap2 to be used
                                         (with map-ont, map-pb for -ont, -pb).
                                         [default: \"minimap2-sr\"]";

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

pub fn contig_full_help() -> &'static str {
    lazy_static! {
        static ref CONTIG_HELP: String = format!(
        "coverm contig: Calculate read coverage per-contig

Define mapping(s) (required):
  Either define BAM:
   -b, --bam-files <PATH> ..             Path to BAM file(s). These must be
                                         reference sorted (e.g. with samtools sort)
                                         unless --sharded is specified, in which
                                         case they must be read name sorted (e.g.
                                         with samtools sort -n).

  Or do mapping:
   -r, --reference <PATH> ..             FASTA file of contigs e.g. concatenated 
                                         genomes or metagenome assembly, or minimap2
                                         index
                                         (with --minimap2-reference-is-index),
                                         or BWA index stem (with -p bwa-mem).
                                         If multiple references FASTA files are
                                         provided and --sharded is specified,
                                         then reads will be mapped to references
                                         separately as sharded BAMs.
   -t, --threads <INT>                   Number of threads for mapping, sorting
                                         and reading.
   -1 <PATH> ..                          Forward FASTA/Q file(s) for mapping
   -2 <PATH> ..                          Reverse FASTA/Q file(s) for mapping
   -c, --coupled <PATH> <PATH> ..        One or more pairs of forward and reverse
                                         FASTA/Q files for mapping in order
                                         <sample1_R1.fq.gz> <sample1_R2.fq.gz>
                                         <sample2_R1.fq.gz> <sample2_R2.fq.gz> ..
   --interleaved <PATH> ..               Interleaved FASTA/Q files(s) for mapping.
   --single <PATH> ..                    Unpaired FASTA/Q files(s) for mapping.
{}
   --minimap2-params PARAMS              Extra parameters to provide to minimap2,
                                         both indexing command (if used) and for
                                         mapping. Note that usage of this parameter
                                         has security implications if untrusted input
                                         is specified. '-a' is always specified.
                                         [default \"\"]
   --minimap2-reference-is-index         Treat reference as a minimap2 database, not 
                                         as a FASTA file.
   --bwa-params PARAMS                   Extra parameters to provide to BWA. Note
                                         that usage of this parameter has security
                                         implications if untrusted input is specified.
                                         [default \"\"]

Sharding i.e. multiple reference sets (optional):
   --sharded                             If -b/--bam-files was used:
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
                                           rpkm
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

Ben J. Woodcroft <benjwoodcroft near gmail.com>
", MAPPER_HELP);
    }
    &CONTIG_HELP
}

pub fn genome_full_help() -> &'static str {
    lazy_static! {
        static ref GENOME_HELP: String = format!(
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
   -d, --genome-fasta-list <PATH>        File containing FASTA file paths, one per 
                                         line
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
{}
   -r, --reference <PATH> ..             FASTA file of contigs e.g. concatenated
                                         genomes or metagenome assembly, or minimap2
                                         index
                                         (with --minimap2-reference-is-index),
                                         or BWA index stem (with -p bwa-mem).
                                         If multiple references FASTA files are
                                         provided and --sharded is specified,
                                         then reads will be mapped to references
                                         separately as sharded BAMs.
   -t, --threads <INT>                   Number of threads for mapping, sorting
                                         and reading.
   -1 <PATH> ..                          Forward FASTA/Q file(s) for mapping
   -2 <PATH> ..                          Reverse FASTA/Q file(s) for mapping
   -c, --coupled <PATH> <PATH> ..        One or more pairs of forward and reverse
                                         FASTA/Q files for mapping in order
                                         <sample1_R1.fq.gz> <sample1_R2.fq.gz>
                                         <sample2_R1.fq.gz> <sample2_R2.fq.gz> ..
   --interleaved <PATH> ..               Interleaved FASTA/Q files(s) for mapping.
   --single <PATH> ..                    Unpaired FASTA/Q files(s) for mapping.
   --minimap2-params PARAMS              Extra parameters to provide to minimap2,
                                         both indexing command (if used) and for
                                         mapping. Note that usage of this parameter
                                         has security implications if untrusted input
                                         is specified. '-a' is always specified.
                                         [default \"\"]
   --minimap2-reference-is-index         Treat reference as a minimap2 database, not 
                                         as a FASTA file.
   --bwa-params PARAMS                   Extra parameters to provide to BWA. Note
                                         that usage of this parameter has security
                                         implications if untrusted input is specified.
                                         [default \"\"]

Sharding i.e. multiple reference sets (optional):
   --sharded                             If -b/--bam-files was used:
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

Dereplication (optional):
   --dereplicate                         Do genome dereplication via average nucleotide
                                         identity (ANI) - choose a genome to represent
                                         all within a small distance, using Dashing for
                                         preclustering and FastANI for final ANI 
                                         calculation.
   --dereplication-ani <FLOAT>           Overall ANI level to dereplicate at with
                                         FastANI.
   --checkm-tab-table                    CheckM tab table for defining genome quality,
                                         which is in turn used during clustering. 
                                         Genomes are scored as 
                                         completeness-4*contamination.
   --genome-info                         dRep style genome info tabl for defining
                                         quality. Used like --checkm-tab-table.
   --min-completeness <FLOAT>            Ignore genomes with less completeness than 
                                         this percentage.
   --max-contamination <FLOAT>           Ignore genomes with more contamination than
                                         this percentage.
   --output-dereplication-clusters <PATH>  Output clustered genome information to this
                                         file as 'representative<TAB>member'
   --dereplication-prethreshold-ani <FLOAT>  Require at least this dashing-derived ANI
                                         for preclustering and to avoid FastANI on
                                         distant lineages within preclusters.
   --dereplication-quality-formula <NAME>  Scoring function for genome quality. See
                                         `coverm dereplicate -h`.
   --dereplication-precluster-method <NAME>  method of calculating rough ANI for 
                                         dereplication. 'dashing' for HyperLogLog, 
                                         'finch' for finch MinHash.

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
                                              rpkm
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
", MAPPER_HELP);
    }
    &GENOME_HELP
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
    --reference assembly.fna --separator '~'

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
                "Example: Map paired reads to a reference where the FASTA header separates
the genome name from the contig name with '~' e.g. >genome10~contig15"
            ),
            ansi_term::Colour::Purple.paint(
                "Example: Calculate coverage of genomes defined as .fna files in
genomes_directory/ from a sorted BAM file:"
            ),
            ansi_term::Colour::Purple.paint(
                "Example: Dereplicate genomes at 99% ANI before mapping unpaired reads:"
            ),
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
    }

    lazy_static! {
        static ref MAKE_HELP: String = format!(
            "coverm make: Generate BAM files through mapping.

Output (required):
   -o, --output-directory <DIR>          Where generated BAM files will go

Mapping parameters:
   -r, --reference <PATH> ..             FASTA file of contigs e.g. concatenated 
                                         genomes or assembly, or minimap2 index
                                         (with --minimap2-reference-is-index),
                                         or BWA index stem (with -p bwa-mem).
                                         If multiple references FASTA files are
                                         provided and --sharded is specified,
                                         then reads will be mapped to references
                                         separately as sharded BAMs.
   -t, --threads <INT>                   Number of threads to use for mapping
   -1 <PATH> ..                          Forward FASTA/Q file(s) for mapping
   -2 <PATH> ..                          Reverse FASTA/Q file(s) for mapping
   -c, --coupled <PATH> <PATH> ..        One or more pairs of forward and reverse
                                         FASTA/Q files for mapping in order
                                         <sample1_R1.fq.gz> <sample1_R2.fq.gz>
                                         <sample2_R1.fq.gz> <sample2_R2.fq.gz> ..
   --interleaved <PATH> ..               Interleaved FASTA/Q files(s) for mapping.
   --single <PATH> ..                    Unpaired FASTA/Q files(s) for mapping.

{}
   --minimap2-params PARAMS              Extra parameters to provide to minimap2,
                                         both indexing command (if used) and for
                                         mapping. Note that usage of this parameter
                                         has security implications if untrusted input
                                         is specified. '-a' is always specified.
                                         [default \"\"]
   --minimap2-reference-is-index         Treat reference as a minimap2 database, not 
                                         as a FASTA file.
   --bwa-params PARAMS                   Extra parameters to provide to BWA. Note
                                         that usage of this parameter has security
                                         implications if untrusted input is specified.
                                         [default \"\"]
   --discard-unmapped                    Exclude unmapped reads from generated BAM files.

Example usage:

  coverm make -r combined_genomes.fna -1 read1.fq -2 read2.fq

Ben J. Woodcroft <benjwoodcroft near gmail.com>
", MAPPER_HELP);
    };

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
                        ])
                        .conflicts_with("bam-files"),
                )
                .arg(
                    Arg::with_name("reference")
                        .short("-r")
                        .long("reference")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless_one(&["bam-files", "full-help"])
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
                .arg(
                    Arg::with_name("output-directory")
                        .short("-o")
                        .long("output-directory")
                        .takes_value(true)
                        .required(true),
                )
                .arg(
                    Arg::with_name("read1")
                        .short("-1")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read2")
                        .required_unless_one(&["coupled", "interleaved", "single"]),
                )
                .arg(
                    Arg::with_name("read2")
                        .short("-2")
                        .multiple(true)
                        .takes_value(true)
                        .requires("read1")
                        .required_unless_one(&["coupled", "interleaved", "single"]),
                )
                .arg(
                    Arg::with_name("coupled")
                        .short("-c")
                        .long("coupled")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&["read1", "interleaved", "single"]),
                )
                .arg(
                    Arg::with_name("interleaved")
                        .long("interleaved")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&["read1", "coupled", "single"]),
                )
                .arg(
                    Arg::with_name("single")
                        .long("single")
                        .multiple(true)
                        .takes_value(true)
                        .required_unless_one(&["read1", "coupled", "interleaved"]),
                )
                .arg(
                    Arg::with_name("reference")
                        .short("-r")
                        .long("reference")
                        .multiple(true)
                        .takes_value(true)
                        .required(true),
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
