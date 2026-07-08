# Changelog

## Unreleased

## Version 0.8.0

- Add strobealign-aemb, rammap and minibwa as mappers. Thanks @PoisonAlien for your suggestion
- Add ANIr coverage method. Thanks @bglindner for the suggestion
- Add mapq-based filtering through `--min-mapq`. Thanks @SDmetagenomics for the suggestion
- Add makedb mode to build persistent mapping indexes (databases) from reference genome or contig FASTA files.
- Add `--bam-file-cache-names` argument for explicit BAM cache filenames
- Various optimisations. Thanks @lavafroth
- Fixed bioconda link. Thanks @ThisIsNotANamepng
- Add demo. Thanks @R-Nurdiansyah for testing it out
- Add fasta.gz support for inputs. Thanks @luispedro, @jianshu93 for the suggestion
- Switch to crates version of which. Thanks @nickp60 for the suggestion
- Block usage of skani with ANI <85% ANI. Thanks @jianshu93 for the suggestion.

## Version 0.7.0

- Add CITATION.cff
- Add `--checkm2-quality-file`
- Use strobealign as default mapper, instead of minimap2
- Update Galah to v0.4.0
- Add minimap2-hifi
- Update to clap command line parser v4
- Support both BWA and BWA-MEM2. Thanks @jianshu93
- genome: Update `--reference` help. Thanks @gallardoalba
- genome: Add `--use-full-contig-names`. Thanks @shafferm
- genome_parsing: Allow compressed genomes
- contig: Always reset after each contig. Thanks @PandengWang
- contig: Fix documented default for `--min-covered-fraction`. Thanks Jiarui Sun
- make: Detect when output filenames clash. Thanks @akiledal
- filter: Detail `--proper-pairs`/`--invert` interactions. Thanks @Rridley7
- release: Drop symbols and use lto. Thanks @jakobnissen
- mapping_index: Allow missing bwa-mem ref when index exists. Thanks @fbeghini

## Version 0.6.1

- Croak when minimap2 finds unequal read pairs. Thanks Robert Hoelzle, Katherine Weigh.
- Standardise documentation to use percentage over fraction (both are still acceptable as parameters though). Thanks Steve Robbins.
- genome: Fix noisy logging.
- genome, contig: Better error msg for bad -r. Thanks @mcmahon-uw Katherine (Trina) McMahon.

## Version 0.6.0

- genome, contig: Add -o/--output-file option. Thanks @michoug
- Suggest a solution for dashing install problems. Thanks @gecko1990
- Add an FAQ section to the manual, showing how to use `$TMPDIR` to change the temporary directory used. Thanks @ShangjinTan
- genome: Fix `-x` so it handles a leading dot. Thanks @mcmahon-uw
- Include examples in the full help and HTML pages. Thanks @mcmahon-uw
- genome: Fix autoconcatenation when contig names clash.
- Add tpm calculation method. Thanks @apcamargo

## Version 0.5.0

- Overhaul the way `--full-help` is displayed by making it a man page, and publishing HTML versions. Add `man` as a dependency.
- dereplication: Update to galah 0.2.0 (this changes dereplication results), and add new ways to output dereplication results and add thresholding options
- Fix bug in `genome` when used with `--no-zeroes` that caused incorrect coverage estimates
- Include supplementary alignments by default (can be reverted with `--exclude-supplementary-alignments`)
- Croak when input BAM file is unsorted
- `genome`: Add `--genome-fasta-list` input option

## Version 0.4.0

- `genome` mode: Incorporate on the fly dereplication via Galah — use `--dereplicate`.
- `dereplicate`: New utility mode
- `contig`: metabat output: Print 4 decimal places. Thanks @apcamargo
- Documentation improvements

## Version 0.3.2

- Fix for samtools 1.10
- More informative error messages when wrapped processes fail
- Logging updates and other minor fixes

## Version 0.3.1

- Fix integer overflow bug for mean and relative_abundance calculators. Thanks @Thexiyang
- Fixes for minimap2 for very large databases
- shell-completion: New mode generates shell completions
- Other bugfixes and documentation updates

## Version 0.3.0

- Default mapper is now minimap2
- Support for Nanopore and PacBio read mapping
- New coverage method RPKM
- Genomes can now be defined with a TSV file with lines `genome_name<tab>contig`
- Various documentation improvements and bugfixes

## Version 0.2.0-alpha7

- Correct bug in counting of total reads when filtering. Thanks Megan Clay.
- Add ability to shard reference sequence database (WIP, experimental). Thanks Rhys Newell.
- Documentation and UI improvements. Thanks Megan Clay, Steven Robbins.

## Version 0.2.0-alpha6

- Fix bug in dense output mode, and make it the default output mode. Thanks Steven Robbins.
- Change help messages so default is to have short colourful output
- Speed up mapping by running samtools sort in a temporary directory
- filter: Introduce --inverse for decontamination uses
- `covered_bases` and `reads_per_base` are new coverage statistics (methods)
- When running in genome mode and mapping, genomes are concatenated automatically so -r isn't required
- In genome mode, when a genome has too small a coverage for `--min-covered-fraction`, reads mapped to that genome are not counted as mapped when calculating relative abundance (same for contig mode)
- Filtering options now have `-read` to reduce potential confusion
- Other small bug and documentation fixes

## Version 0.2.0-alpha4

- Improve help messages
- Fix bug in dense output mode for genomes
