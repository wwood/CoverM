# CoverM

[![Travis](https://img.shields.io/travis/wwood/CoverM.svg?style=flat-square)](https://travis-ci.org/wwood/CoverM)

CoverM aims to be a configurable, easy to use and fast DNA read coverage and
relative abundance calculator focused on metagenomics applications.

CoverM calculates coverage of genomes/MAGs (`coverm genome`) or individual
contigs (`coverm contig`). Calculating coverage by read mapping, its input can
either be BAM files sorted by reference, or raw reads and reference FASTA
sequences.

## Installation

### Pre-compiled binary

To install CoverM, the most straightforward way is to download the statically
compiled binaries available on the [releases page](https://github.com/wwood/CoverM/releases).

### Compiling from source

CoverM can also be installed from source, using the cargo build system after
installing [Rust](https://www.rust-lang.org/).

```
cargo install coverm
```

### Dependencies
For the full suite of options, these additional programs must be installed:

* [samtools](https://github.com/samtools/samtools) >1.0 (tested with v1.9)
* [bwa](https://github.com/lh3/bwa)
* [tee](https://www.gnu.org/software/coreutils/), which is installed by default
  on most Linux operating systems.

## Usage
```
Main subcommands:
	contig	Calculate coverage of contigs
	genome	Calculate coverage of genomes

Less used utility subcommands:
	make	Generate BAM files through alignment
	filter	Remove (or only keep) alignments with insufficient identity

Other options:
	-V, --version	Print version information
```

For more detailed usage see `coverm <command> -h` or `coverm <command> --full-help`.

## Calculation methods

The `-m/--methods` flag specified the specific kind(s) of coverage that is 
to be calculated.

To illustrate, imagine a set of 3 pairs of reads, where only 1 aligns to a
single reference contig of length 1000bp:

```
read1_forward    ==================>
read1_reverse                         <==================
contig    ...-----------------------------------------------------....
                 |        |         |         |         |
position        200      210       220       230       240
```
The difference coverage measures would be:

| Method | Value | Formula | Explanation |
|--------------------|------------|-------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| mean | 0.02235294 | (10+9)/(1000-2*75) | The two reads have 10 and 9 bases aligned exactly, averaged over 1000-2*75 bp (length of contig minus 75bp from each end). |
| relative_abundance | 33.3% | 0.02235294/0.02235294*(2/6) | If the contig is considered a genome, then its mean coverage is 0.02235294. There is a total of 0.02235294 mean coverage across all genomes, and 2 out of 6 reads (1 out of 3 pairs) map. This coverage calculation is only available in 'genome' mode. |
| trimmed_mean | 0 | mean_coverage(mid-ranked-positions) | After removing the 5% of bases with highest coverage and 5% of bases with lowest coverage, all remaining positions have coverage 0. |
| covered_fraction | 0.02 | (10+10)/1000 | 20 bases are covered by any read, out of 1000bp. |
| covered_bases | 20 | 10+10 | 20 bases are covered. |
| variance | 0.01961962 | var({1;20},{0;980}) | Variance is calculated as the sample variance. |
| length | 1000 |  | The contig's length is 1000bp. |
| count | 2 |  | 2 reads are mapped. |
| reads_per_base | 0.002 | 2/1000 | 2 reads are mapped over 1000bp. |
| metabat | contigLen 1000, totalAvgDepth 0.02235294, bam depth 0.02235294, variance 0.01961962 | | Reproduction of the [MetaBAT](https://bitbucket.org/berkeleylab/metabat) 'jgi_summarize_bam_contig_depths' tool output, producing [identical output](https://bitbucket.org/berkeleylab/metabat/issues/48/jgi_summarize_bam_contig_depths-coverage). |
| coverage_histogram | 20 bases with coverage 1, 980 bases with coverage 0 | | The number of positions with each different coverage are tallied. |

Calculation of genome-wise coverage (`genome` mode) is similar to calculating
contig-wise (`contig` mode) coverage, except that the unit of reporting is
per-genome rather than per-contig. For calculation methods which exclude base
positions based on their coverage, all positions from all contigs are considered
together. For instance, if a 2000bp contig with all positions having 1X coverage
is in a genome with 2,000,000bp contig with no reads mapped, then the
trimmed_mean will be 0 as all positions in the 2000bp are in the top 5% of
positions sorted by coverage.

## License

CoverM is made available under GPL3+. See LICENSE.txt for details. Copyright Ben
Woodcroft.

Developed by Ben Woodcroft at the [Australian Centre for Ecogenomics](http://ecogenomic.org).
