# CoverM

CoverM is a read coverage calculator focused on metagenomics applications. It
calculates coverage from sorted BAM files, comprised of reads mapped to either a
collection of genomes, or to each contig individually.

## Installation

To install CoverM, first install [Rust](https://www.rust-lang.org/). CoverM is currently tested on Rust
1.25.0, but likely works on other versions also.

To install CoverM itself:
```
cargo install coverm
```

## Usage
See `coverm -h`.

## License

CoverM is made available under GPL3+. See LICENSE.txt for details. Copyright Ben
Woodcroft.

Developed at the [Australian Centre for Ecogenomics](http://ecogenomic.org).
