# CoverM

CoverM aims to be a configurable, easy to use and fast read coverage calculator focused on metagenomics applications. 

Calculating coverage by read mapping, its input can either be BAM files sorted by reference, or raw reads and reference FASTA sequences.

CoverM calculates the coverage either or genomes/MAGs (`coverm genome`) or individual contigs (`coverm contig`).

## Installation

To install CoverM, the most straightforward way is to download the statically compiled binaries available on the [releases page](https://github.com/wwood/CoverM/releases).

CoverM can also be installed from source, using the cargo build system after installing [Rust](https://www.rust-lang.org/). CoverM is not currently on crates.io.

## Usage
```
Main modes:
	contig	Calculate coverage of contigs
	genome	Calculate coverage of genomes

Utilities:
	make	Generate BAM files through alignment
	filter	Remove alignments with insufficient identity

Other options:
	-V, --version	Print version information
```

For detailed usage see `coverm <command> -h`.

## License

CoverM is made available under GPL3+. See LICENSE.txt for details. Copyright Ben
Woodcroft.

Developed by Ben Woodcroft at the [Australian Centre for Ecogenomics](http://ecogenomic.org).
