NAME
====

coverm genome - Calculate read coverage per-genome

SYNOPSIS
========

**coverm genome** \[FLAGS\] \[OPTIONS\]

FLAGS
=====

**\--single-genome**

:   All contigs are from the same genome

<!-- -->

**\--minimap2-reference-is-index **

:   Treat reference as a minimap2 database, not as a FASTA file.

<!-- -->

**\--sharded**

:   If -b/\--bam-files was used: Input BAM files are read-sorted
    alignments of a set of reads mapped to multiple reference contig
    sets. Choose the best hit for each read pair. Otherwise if mapping
    was carried out: Map reads to each reference, choosing the best hit
    for each pair.

<!-- -->

**\--exclude-genomes-from-deshard**

:   Ignore genomes whose name appears in this newline-separated file
    when combining shards.

<!-- -->

**\--proper-pairs-only**

:   Require reads to be mapped as proper pairs

<!-- -->

**\--exclude-supplementary**

:   Exclude supplementary alignments

<!-- -->

**\--dereplicate**

:   Do genome dereplication via average nucleotide identity (ANI) -
    choose a genome to represent all within a small distance, using
    Dashing for preclustering and FastANI for final ANI calculation.

<!-- -->

**\--no-zeros**

:   Omit printing of genomes that have zero coverage

<!-- -->

**\--discard-unmapped**

:   Exclude unmapped reads from cached BAM files.

<!-- -->

**-v**, **\--verbose**

:   Print extra debugging information

<!-- -->

**-q**, **\--quiet**

:   Unless there is an error, do not print log messages

OPTIONS
=======

**-s**, **\--separator**=*CHARACTER*

:   This character separates genome names from contig names in the
    reference file

<!-- -->

**-f**, **\--genome-fasta-files**=*PATH ..*

:   Path(s) to FASTA files of each genome e.g. \'pathA/genome1.fna
    pathB/genome2.fa\'

<!-- -->

**-d**, **\--genome-fasta-directory**=*PATH*

:   Directory containing FASTA files of each genome

<!-- -->

**-x**, **\--genome-fasta-extension**=*EXT*

:   File extension of genomes in the directory specified with
    -d/\--genome-fasta-directory \[default \"fna\"\]

<!-- -->

**-d**, **\--genome-fasta-list**=*PATH*

:   File containing FASTA file paths, one per line

<!-- -->

**\--genome-definition**=*FILE*

:   File containing list of genome\_name\<tab\>contig lines to define
    the genome of each contig

<!-- -->

**-b**, **\--bam-files**=*PATH*

:   Path to BAM file(s). These must be reference sorted (e.g. with
    samtools sort) unless \--sharded is specified, in which case they
    must be read name sorted (e.g. with samtools sort -n).

<!-- -->

**-p**, **\--mapper**=*NAME*

:   Underlying mapping software used (\"minimap2-sr\", \"bwa-mem\",
    \"minimap2-ont\", \"minimap2-pb\", or \"minimap2-no-preset\").
    minimap2 -sr, -ont, -pb, -no-preset specify \'-x\' preset of
    minimap2 to be used (with map-ont, map-pb for -ont, -pb). \[default:
    \"minimap2-sr\"\]

<!-- -->

**-r**, **\--reference**=*PATH*

:   FASTA file of contigs e.g. concatenated genomes or metagenome
    assembly, or minimap2 index (with \--minimap2-reference-is-index),
    or BWA index stem (with -p bwa-mem). If multiple references FASTA
    files are provided and \--sharded is specified, then reads will be
    mapped to references separately as sharded BAMs.

<!-- -->

**-t**, **\--threads**=*INT*

:   Number of threads for mapping, sorting and reading.

<!-- -->

**-1**=*PATH ..*

:   Forward FASTA/Q file(s) for mapping

<!-- -->

**-2**=*PATH ..*

:   Reverse FASTA/Q file(s) for mapping

<!-- -->

**-c**, **\--coupled**=*PATH ..*

:   One or more pairs of forward and reverse FASTA/Q files for mapping
    in order \<sample1\_R1.fq.gz\> \<sample1\_R2.fq.gz\>
    \<sample2\_R1.fq.gz\> \<sample2\_R2.fq.gz\> ..

<!-- -->

**\--interleaved**=*PATH ..*

:   Interleaved FASTA/Q files(s) for mapping.

<!-- -->

**\--single**=*PATH ..*

:   Unpaired FASTA/Q files(s) for mapping.

<!-- -->

**\--minimap2-params**=*PARAMS*

:   Extra parameters to provide to minimap2, both indexing command (if
    used) and for mapping. Note that usage of this parameter has
    security implications if untrusted input is specified. \'-a\' is
    always specified. \[default \"\"\]

<!-- -->

**\--bwa-params**=*PARAMS*

:   Extra parameters to provide to BWA. Note that usage of this
    parameter has security implications if untrusted input is specified.
    \[default \"\"\]

<!-- -->

**\--min-read-aligned-length**=*INT*

:   Exclude reads with smaller numbers of aligned bases \[default: 0\]

<!-- -->

**\--min-read-percent-identity**=*FLOAT*

:   Exclude reads by overall percent identity e.g. 0.95 for 95%.
    \[default 0.0\]

<!-- -->

**\--min-read-aligned-percent**=*FLOAT*

:   Exclude reads by percent aligned bases e.g. 0.95 means 95% of the
    read\'s bases must be aligned. \[default 0.0\]

<!-- -->

**\--min-read-aligned-length-pair**=*INT*

:   Exclude pairs with smaller numbers of aligned bases. Implies
    \--proper-pairs-only. \[default: 0\]

<!-- -->

**\--min-read-percent-identity-pair**=*FLOAT*

:   Exclude pairs by overall percent identity e.g. 0.95 for 95%. Implies
    \--proper-pairs-only. \[default 0.0\]

<!-- -->

**\--min-read-aligned-percent-pair**=*FLOAT*

:   Exclude reads by percent aligned bases e.g. 0.95 means 95% of the
    read\'s bases must be aligned. Implies \--proper-pairs-only.
    \[default 0.0\]

<!-- -->

**\--dereplication-ani**=*FLOAT*

:   Overall ANI level to dereplicate at with FastANI.

<!-- -->

**\--checkm-tab-table**=*PATH*

:   CheckM tab table for defining genome quality, which is in turn used
    during clustering. Genomes are scored as
    completeness-4\*contamination.

<!-- -->

**\--genome-info**=*PATH*

:   dRep style genome info tabl for defining quality. Used like
    \--checkm-tab-table.

<!-- -->

**\--min-completeness**=*FLOAT*

:   Ignore genomes with less completeness than this percentage.

<!-- -->

**\--max-contamination**=*FLOAT*

:   Ignore genomes with more contamination than this percentage.

<!-- -->

**\--output-dereplication-clusters**=*PATH*

:   Output clustered genome information to this file as
    \'representative\<TAB\>member\'

<!-- -->

**\--dereplication-prethreshold-ani**=*FLOAT*

:   Require at least this dashing-derived ANI for preclustering and to
    avoid FastANI on distant lineages within preclusters.

<!-- -->

**\--dereplication-quality-formula**=*NAME*

:   Scoring function for genome quality. See \`coverm dereplicate -h\`.

<!-- -->

**\--dereplication-precluster-method**=*NAME*

:   method of calculating rough ANI for dereplication. \'dashing\' for
    HyperLogLog, \'finch\' for finch MinHash.

<!-- -->

**-m**, **\--methods**=*METHOD*

:   Method(s) for calculating coverage. One or more (space separated)
    of: relative\_abundance (default), mean, trimmed\_mean,
    coverage\_histogram, covered\_fraction, covered\_bases, variance,
    length, count, reads\_per\_base, rpkm. A more thorough description
    of the different methods is available at
    https://github.com/wwood/CoverM

<!-- -->

**\--output-format**=*FORMAT*

:   Shape of output: \'sparse\' for long format, \'dense\' for
    species-by-site. \[default: dense\]

<!-- -->

**\--min-covered-fraction**=*FRACTION*

:   Genomes with less coverage than this reported as having zero
    coverage. \[default: 0.10\]

<!-- -->

**\--contig-end-exclusion**=*INT*

:   Exclude bases at the ends of reference sequences from calculation
    \[default: 75\]

<!-- -->

**\--trim-min**=*FRACTION*

:   Remove this smallest fraction of positions when calculating
    trimmed\_mean \[default: 0.05\]

<!-- -->

**\--trim-max**=*FRACTION*

:   Maximum fraction for trimmed\_mean calculations \[default: 0.95\]

<!-- -->

**\--bam-file-cache-directory**=*DIRECTORY*

:   Output BAM files generated during alignment to this directory. The
    directory may or may not exist

EXIT STATUS
===========

**0**

:   Successful program execution.

<!-- -->

**1**

:   Unsuccessful program execution.

<!-- -->

**101**

:   The program panicked.

AUTHOR
======

>     Ben J Woodcroft <benjwoodcroft near gmail.com>
