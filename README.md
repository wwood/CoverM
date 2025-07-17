![CoverM logo](images/coverm.png)

- [CoverM](#coverm)
	- [Installation](#installation)
		- [Install through the bioconda package](#install-through-the-bioconda-package)
		- [Pre-compiled binary](#pre-compiled-binary)
		- [Compiling from source](#compiling-from-source)
		- [Development version](#development-version)
		- [Dependencies](#dependencies)
		- [Shell completion](#shell-completion)
	- [Usage](#usage)
	- [Demo](#demo)
	- [Calculation methods](#calculation-methods)
	- [Citation](#citation)
	- [License](#license)

# CoverM

[![Anaconda-Server Badge](https://anaconda.org/bioconda/coverm/badges/version.svg)](https://anaconda.org/bioconda/coverm)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/coverm/badges/downloads.svg)](https://anaconda.org/bioconda/coverm)

CoverM aims to be a configurable, easy to use and fast DNA read coverage and
relative abundance calculator focused on metagenomics applications.

CoverM calculates coverage of genomes/MAGs `coverm genome` ([help](https://wwood.github.io/CoverM/coverm-genome.html)) or individual
contigs `coverm contig` ([help](https://wwood.github.io/CoverM/coverm-contig.html)). Calculating coverage by read mapping, its input can
either be BAM files sorted by reference, or raw reads and reference genomes in various formats.

## Installation

### Install through the bioconda package

CoverM and its dependencies can be installed through the [bioconda](https://bioconda.github.io/user/install.html) conda channel. After initial setup of conda and the bioconda channel, it can be installed with

```
conda install coverm
```

### Pre-compiled binary

Statically compiled CoverM binaries available on the [releases page](https://github.com/wwood/CoverM/releases).
This installation method requires non-Rust dependencies to be installed separately - see the [dependencies section](#Dependencies).

### Compiling from source

CoverM can also be installed from source, using the cargo build system after
installing [Rust](https://www.rust-lang.org/).

```
cargo install coverm
```

### Development version
To run an unreleased version of CoverM, after installing
[Rust](https://www.rust-lang.org/) and any additional dependencies listed below:

```
git clone https://github.com/wwood/CoverM
cd CoverM
cargo run -- genome ...etc...
```

To run tests:

```
cargo build
cargo test
```

### Dependencies
For the full suite of options, additional programs must also be installed, when
installing from source or for development.

These can be installed using the conda YAML environment definition:
```
conda env create -n coverm -f coverm.yml
```

Or, these can be installed manually:

* [samtools](https://github.com/samtools/samtools) v1.9
* [tee](https://www.gnu.org/software/coreutils/), which is installed by default
  on most Linux operating systems.
* [man](http://man-db.nongnu.org/), which is installed by default on most Linux
  operating systems.

and one of these for mapping:
* [strobealign](https://github.com/ksahlin/StrobeAlign) v0.14.0
* [minimap2](https://github.com/lh3/minimap2) v2.21
* [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) v2.0

and one of these for genome dereplication:
* [skani](https://github.com/bluenote-1577/skani) v0.1.1
* [Dashing](https://github.com/dnbaker/dashing) v0.4.0
* [FastANI](https://github.com/ParBLiSS/FastANI) v1.3

### Shell completion
Completion scripts for various shells e.g. BASH can be generated. For example, to install the bash completion script system-wide (this requires root privileges):

```
coverm shell-completion --output-file coverm --shell bash
mv coverm /etc/bash_completion.d/
```

It can also be installed into a user's home directory (root privileges not required):

```
coverm shell-completion --shell bash --output-file /dev/stdout >>~/.bash_completion
```

In both cases, to take effect, the terminal will likely need to be restarted. To test, type `coverm gen` and it should complete after pressing the TAB key.

## Usage

CoverM operates in several modes. Detailed usage information including examples is given at the links below, or alternatively by using the `-h` or `--full-help` flags for each mode:
* [genome](https://wwood.github.io/CoverM/coverm-genome.html) - Calculate coverage of genomes
* [contig](https://wwood.github.io/CoverM/coverm-contig.html) - Calculate coverage of contigs

There are several utility modes as well:
* [make](https://wwood.github.io/CoverM/coverm-make.html) - Generate BAM files through alignment
* [filter](https://wwood.github.io/CoverM/coverm-filter.html) - Remove (or only keep) alignments with insufficient identity
* [cluster](https://wwood.github.io/CoverM/coverm-cluster.html) - Dereplicate and cluster genomes
* shell-completion - Generate shell completion scripts

## Demo

A common use case for CoverM is to calculate the coverage or relative abundance of a set of genomes in a metagenomic sample.
This can be useful for assessing the proportion of the diversity in the sample that is represented by the genomes, and for comparing the abundance of genomes across samples.
Here we will demonstrate how to calculate the coverage of 8 genomes in a sample of paired-end reads.
We will assume that you have already installed CoverM and its dependencies (see above if not).

First we will download the test dataset of 8 genomes and 1 sample of paired-end reads.
This sample and genomes come from permafrost samples in Abisko, Sweden (see [doi](https://doi.org/10.1101/2025.02.07.636677) for more details).
The sample has been sub-sampled to 100,000 reads to speed up this demo.

```bash
# Create a directory for the demo
mkdir ~/coverm_demo
cd ~/coverm_demo

# Download files
wget https://raw.githubusercontent.com/wwood/CoverM/refs/heads/main/demo/sample_1.1.fq.gz
wget https://raw.githubusercontent.com/wwood/CoverM/refs/heads/main/demo/sample_1.2.fq.gz
wget https://raw.githubusercontent.com/wwood/CoverM/refs/heads/main/demo/genome_1.fna
wget https://raw.githubusercontent.com/wwood/CoverM/refs/heads/main/demo/genome_2.fna
wget https://raw.githubusercontent.com/wwood/CoverM/refs/heads/main/demo/genome_3.fna
wget https://raw.githubusercontent.com/wwood/CoverM/refs/heads/main/demo/genome_4.fna
wget https://raw.githubusercontent.com/wwood/CoverM/refs/heads/main/demo/genome_5.fna
wget https://raw.githubusercontent.com/wwood/CoverM/refs/heads/main/demo/genome_6.fna
wget https://raw.githubusercontent.com/wwood/CoverM/refs/heads/main/demo/genome_7.fna
wget https://raw.githubusercontent.com/wwood/CoverM/refs/heads/main/demo/genome_8.fna
```

Next, we can actually run CoverM.
We use `--coupled` to specify the paired-end reads, and `--genome-fasta-files` to specify the genomes.
We also specify the number of threads to use with `-t`, and the metrics we want to calculate with `-m`.
We will calculate the mean coverage, relative abundance, and covered fraction for each genome in the sample.
The output will be written to `output_coverm.tsv`.

```bash
coverm genome \
  --coupled sample_1.1.fq.gz sample_1.2.fq.gz \
  --genome-fasta-files \
    genome_1.fna genome_2.fna genome_3.fna genome_4.fna \
    genome_5.fna genome_6.fna genome_7.fna genome_8.fna \
  -t 8 \
  -m mean relative_abundance covered_fraction \
  -o output_coverm.tsv
```

This should have created the file `output_coverm.tsv` and logged the following message:
`coverm::genome] In sample 'sample_1.1.fq.gz', found 48254 reads mapped out of 100000 total (48.25%)`.
This indicates that 48.25% of the reads from our sample mapped to the genomes. So our genomes represent about half of the diversity in the sample.

The file `output_coverm.tsv` should look like this:

| Genome   | sample_1.1.fq.gz Mean | sample_1.1.fq.gz Relative Abundance (%) | sample_1.1.fq.gz Covered Fraction |
|----------|-----------------------|-----------------------------------------|-----------------------------------|
| unmapped | NA                    | 51.746                                  | NA                                |
| genome_1 | 0.9410575             | 25.87694                                | 0.52770287                        |
| genome_2 | 0.40274984            | 11.074703                               | 0.27789244                        |
| genome_3 | 0.20988818            | 5.7714467                               | 0.15818907                        |
| genome_4 | 0.20114066            | 5.5309105                               | 0.1509256                         |
| genome_5 | 0                     | 0                                       | 0                                 |
| genome_6 | 0                     | 0                                       | 0                                 |
| genome_7 | 0                     | 0                                       | 0                                 |
| genome_8 | 0                     | 0                                       | 0                                 |

Looking in `output_coverm.tsv`, we find columns with the following headings:

- `Genome`: The name of the genome
- `sample_1.1.fq.gz Mean`: The mean read coverage from sample_1 across the given genome, i.e. the average height across the genome if reads aligned were stacked on top of each other.
- `sample_1.1.fq.gz Relative Abundance (%)`: The relative abundance of the genome within sample_1. This metric accounts for differing genome sizes by using the proportion of mean coverage rather than the proportion of reads.
- `sample_1.1.fq.gz Covered Fraction`: The proportion of the genome that is covered by at least one read.

Each row represents a genome, and the columns represent the metrics calculated for that genome for each provided sample.
For instance, the row for `genome_1` shows that the mean coverage of this genome is `0.941`, the relative abundance is `25.9`%, and the covered fraction is `0.528`.
Again, the row for `genome_5` shows that the mean coverage of this genome is `0.0`, the relative abundance is `0.0`%, and the covered fraction is `0.0`.
This indicates that `genome_1` is well represented in the sample, while `genome_5` is not present at all.
There are 3 other genomes with varying coverage, and 3 other genomes with 0 coverage.

You may have noticed that the covered fraction for most genomes is rather low. This is because the reads have been sub-sampled to 100,000 reads.
The full sample has 76,618,686 reads and produces covered fractions of 1 for all present genomes.
In contrast, the relative abundance estimations do not change much when calculated with the full set of reads.
The output from the full sample can be found [here](https://github.com/wwood/CoverM/blob/main/demo/output_coverm_full.tsv), or see below.

| Genome   | 713_CPN1-1-X0.1.fq.gz Mean | 713_CPN1-1-X0.1.fq.gz Relative Abundance (%) | 713_CPN1-1-X0.1.fq.gz Covered Fraction |
|----------|----------------------------|----------------------------------------------|----------------------------------------|
| unmapped | NA                         | 50.81075                                     | NA                                     |
| genome_1 | 732.317                    | 26.360964                                    | 0.999935                               |
| genome_2 | 304.20364                  | 10.950313                                    | 1                                      |
| genome_3 | 165.43732                  | 5.9551897                                    | 0.99993813                             |
| genome_4 | 159.31432                  | 5.734782                                     | 0.9999948                              |
| genome_5 | 0                          | 0                                            | 0                                      |
| genome_6 | 0                          | 0                                            | 0                                      |
| genome_7 | 0                          | 0                                            | 0                                      |
| genome_8 | 5.2227674                  | 0.18800214                                   | 0.29310852                             |

There is an additional row named `unmapped` which represents the metrics for the reads that did not map to any of the provided genomes.
This is only applicable to the relative abundance metric (among those we selected), and we can see that 51% of the reads were unmapped.

## Calculation methods

The `-m/--methods` flag specifies the specific kind(s) of coverage that are
to be calculated.

To illustrate, imagine a set of 3 pairs of reads, where only 1 aligns to a
single reference contig of length 1000bp:

```
read1_forward    ========>
read1_reverse                                  <====+====
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
| anir | 0.95 | (sum identity of reads)/number of reads | Average BLAST-like identity of mapped reads |
| metabat | contigLen 1000, totalAvgDepth 0.02235294, bam depth 0.02235294, variance 0.01961962 | | Reproduction of the [MetaBAT](https://bitbucket.org/berkeleylab/metabat) 'jgi_summarize_bam_contig_depths' tool output, producing [identical output](https://bitbucket.org/berkeleylab/metabat/issues/48/jgi_summarize_bam_contig_depths-coverage). |
| coverage_histogram | 20 bases with coverage 1, 980 bases with coverage 0 | | The number of positions with each different coverage are tallied. |
| rpkm | 1000000 | 2 * 10^9 / 1000 / 2 | Calculation here assumes no other reads map to other contigs. See https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/ for an explanation of RPKM and TPM. Note that this calculates the RPKM contigs (or genomes in `genome` mode), not individual genes.|
| tpm | 1000000 | rpkm/total_of_rpkm * 10^6 | Calculation here assumes no other reads map to other contigs. See RPKM above. Note that this calculates the TPM contigs (or genomes in `genome` mode), not individual genes.|

Calculation of genome-wise coverage (`genome` mode) is similar to calculating
contig-wise (`contig` mode) coverage, except that the unit of reporting is
per-genome rather than per-contig. For calculation methods which exclude base
positions based on their coverage, all positions from all contigs are considered
together. For instance, if a 2000bp contig with all positions having 1X coverage
is in a genome with 2,000,000bp contig with no reads mapped, then the
trimmed_mean will be 0 as all positions in the 2000bp are in the top 5% of
positions sorted by coverage.

## Citation

If you use CoverM in your research, please cite the following publication:

Aroney, S.T., Newell, R.J., Nissen, J.N., Camargo, A.P., Tyson, G.W. and Woodcroft, B.J., 2025. CoverM: Read alignment statistics for metagenomics, Bioinformatics 41:4, https://doi.org/10.1093/bioinformatics/btaf147.

## License

CoverM is made available under GPL3+. See LICENSE.txt for details. Copyright Ben
Woodcroft.

Developed by Ben Woodcroft at the Queensland University of Technology [Centre for Microbiome Research](https://research.qut.edu.au/cmr/).
