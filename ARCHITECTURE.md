# CoverM Architecture

## Overview

CoverM calculates read-mapping coverage and relative abundance of genomes/MAGs
(`coverm genome`) or individual contigs (`coverm contig`) for metagenomics. Its
core approach is a streaming, single-pass scan over a reference-sorted BAM file:
per reference sequence it builds a "ups-and-downs" delta array (mosdepth-style
coverage deltas — `+1` where an alignment starts, `-1` where it ends) and feeds
that array into a set of pluggable coverage estimators. Input can be pre-made
BAM files or raw reads, in which case CoverM shells out to a mapper (BWA,
minimap2, strobealign) and consumes its SAM output through a FIFO without ever
writing an intermediate BAM to disk.

## Module Map

Start at `src/bin/coverm.rs`, then follow the trait triad (reader → estimator →
taker → printer).

- **`src/bin/coverm.rs`** — The CLI binary and top-level orchestrator. `main()`
  dispatches on subcommand (`genome`, `contig`, `filter`, `make`, `makedb`,
  `cluster`, `shell-completion`). Key helpers: `run_genome(...)` (the genome
  workflow, called once per genome-definition strategy), `EstimatorsAndTaker`
  (bundles the chosen estimators + coverage taker + printer, built by
  `generate_from_clap`), `setup_mapping_index`, `dereplicate` (delegates to
  galah). This is where clap args become concrete estimator/reader choices.

- **`src/cli.rs`** — Pure clap command definitions, help text, and man-page
  generation (via `bird_tool_utils_man`). Large but mechanical; edit here to
  add/rename flags. Depended on by the binary only.

- **`src/bam_generator.rs`** — The heart of read ingestion. Defines the
  `NamedBamReader` and `NamedBamReaderGenerator<T>` traits and their
  implementations: `BamFileNamedReader` (pre-made BAM), `StreamingNamedBamReader`
  (spawns a mapper, reads its SAM via FIFO), and the filtered variants. The
  `MappingProgram` enum enumerates supported mappers. `generate_*` free functions
  are the factory API the binary calls.

- **`src/mosdepth_genome_coverage_estimators.rs`** — The `CoverageEstimator`
  enum: every coverage metric (Mean, TrimmedMean, Variance, RPKM, TPM,
  CoverageFraction, PileupCounts, ReadCount, etc.) is one variant. Each knows how
  to `add_contig(ups_and_downs, ...)`, `calculate_coverage(...)`, `setup()`
  (reset), and emit its `column_headers()`. This is where a new coverage metric
  is added.

- **`src/contig.rs`** — `contig_coverage(...)`: the per-contig scan loop. Reads
  records, maintains the ups-and-downs array per reference, and on reference
  change flushes accumulated values through the estimators into the coverage
  taker.

- **`src/genome.rs`** — The genome-mode equivalents:
  `mosdepth_genome_coverage_with_contig_names` (contig→genome mapping supplied)
  and `mosdepth_genome_coverage` (separator-based). Also contains the
  reads→mapping genome path. Aggregates per-contig coverage up to genome level.

- **`src/coverage_takers.rs`** — `CoverageTaker` trait + `CoverageTakerType`
  enum. A "taker" receives coverage values as the scan produces them. Variants
  either stream straight to output
  (`SingleFloatCoverageStreamingCoveragePrinter`) or cache into memory
  (`CachedSingleFloatCoverageTaker`) for output modes needing all values at once
  (dense matrix, metabat, TPM normalisation).

- **`src/coverage_printer.rs`** — `CoveragePrinter` enum: turns a (usually
  cached) taker into final formatted output — sparse, dense matrix, or
  MetaBAT-adjusted. `finalise_printing` handles RPKM/TPM column normalisation.

- **`src/filter.rs`** — `ReferenceSortedBamFilter`: wraps a BAM reader and drops
  reads/pairs failing alignment-length / percent-identity / MAPQ thresholds,
  pairing up mates using a `first_set` BTreeMap keyed by read name. Backs both
  the `filter` subcommand and inline filtering during coverage.

- **`src/mapping_parameters.rs`** — `MappingParameters` / `ReadFormat`: parses
  read-file arguments (coupled `-1/-2`, interleaved, single) into an iterable of
  per-reference, per-read-set mapping jobs.

- **`src/mapping_index_maintenance.rs`** — `MappingIndex` trait with
  `VanillaIndexStruct` (pre-existing index) and `TemporaryIndexStruct` (builds a
  mapper index in a temp dir that lives as long as the struct). Used before
  streaming mapping.

- **`src/genomes_and_contigs.rs`** — `GenomesAndContigs`: the contig-name →
  genome-index lookup (a `Vec<String>` of genome names + a `HashMap`). Shared
  verbatim with the cockatoo project (see file header). Serde-serializable for
  `makedb`.

- **`src/genome_parsing.rs`** — Builds `GenomesAndContigs` from genome FASTA
  files (`read_genome_fasta_files`) or a genome-definition file
  (`read_genome_definition_file`).

- **`src/genome_exclusion.rs`** — `GenomeExclusion` trait
  (`SeparatorGenomeExclusionFilter`, `GenomesAndContigsExclusionFilter`,
  `NoExclusionGenomeFilter`): decides whether a contig should be skipped, chiefly
  in sharded reading.

- **`src/shard_bam_reader.rs`** — `ReadSortedShardedBamReader`: merges multiple
  BAMs (one per reference shard) by picking, per read, the best-scoring alignment
  (via the `AS` tag) across shards while respecting genome exclusion. Used for
  the "shard" mapping strategy against large references.

- **`src/genes.rs`** — Per-gene coverage: `Gene` / `GeneDefinitions` parse a GFF,
  then coverage is computed over gene intervals rather than whole contigs.

- **`src/strobealign_aemb.rs`** — Special path for strobealign's `--aemb` mode,
  which emits abundance directly; bypasses the ups-and-downs machinery and parses
  strobealign's CSV output.

- **`src/external_command_checker.rs`** — Verifies external tools (bwa, minimap2,
  samtools, strobealign) exist and meet version requirements before use.

## Key Types & Data Flow

The design is a pipeline of four trait/enum abstractions, decoupled so any reader
can feed any estimator into any taker into any printer:

1. **`NamedBamReader` / `NamedBamReaderGenerator<T>`** (`bam_generator.rs`) — the
   source. A generator is `start()`ed into a reader that yields
   `bam::record::Record`s plus a header. Implementations hide whether records
   come from a file, a live mapper subprocess, a filter, or a shard merger.

2. **`CoverageEstimator`** (`mosdepth_genome_coverage_estimators.rs`) — the
   metric. Accumulates state from ups-and-downs arrays across contigs and
   produces one or more `f32` coverage values. Multiple estimators run at once.

3. **`CoverageTaker` / `CoverageTakerType`** (`coverage_takers.rs`) — the sink.
   Receives `(stoit, entry, coverage-values)` tuples in scan order, either
   streaming them out or caching them.

4. **`CoveragePrinter`** (`coverage_printer.rs`) — final formatting of a cached
   taker, including cross-column normalisation (RPKM/TPM).

Supporting types: **`GenomesAndContigs`** maps contigs to genomes;
**`FlagFilter`** (in `lib.rs`) is the SAM-flag gate (secondary / supplementary /
improper-pair); **`ReadsMapped`** carries mapped/total read counts per sample for
normalisation.

**Data flow (reads → coverage):** raw reads + reference → `setup_mapping_index`
builds a `MappingIndex` → a `StreamingNamedBamReaderGenerator` spawns the mapper
writing SAM into a FIFO → records stream into `contig_coverage` /
`mosdepth_genome_coverage_*` → `FlagFilter` and optional
`ReferenceSortedBamFilter` gate each record → per-reference an ups-and-downs
`Vec<i32>` is built; on reference change it is handed to every
`CoverageEstimator` via `add_contig`, then `calculate_coverage` → results pushed
to the `CoverageTaker` → after the scan the `CoveragePrinter` finalises output.
The pre-made-BAM path is identical minus the mapping step.

## Design Decisions & Invariants

- **Ups-and-downs over full pileup.** Coverage is derived from a delta array
  (start `+1`, end `-1`, prefix-summed lazily inside estimators) rather than a
  per-base depth array, which keeps memory at one `i32` per reference base and
  the scan single-pass. (Inferred from the algorithm; the module is named after
  mosdepth, which uses the same trick.)

- **Reference-sorted BAM is a hard invariant.** Both the coverage scan and the
  filter assume records arrive grouped/ordered by reference `tid`; a contig is
  considered "done" the moment a record with a different `tid` appears. Feeding
  an unsorted or name-sorted BAM produces silently wrong results. The streaming
  mappers are configured to emit reference-sorted output for this reason.

- **Streaming via FIFO, no intermediate BAM.** `StreamingNamedBamReader` uses a
  named pipe (`nix` mkfifo) and holds child-process handles + temp log files so
  it can surface mapper stderr on failure. The `TempDir`/`NamedTempFile` fields
  are load-bearing: they must stay in scope for the pipe and index to remain
  valid — do not "simplify" them away.

- **Estimators are cloned per genome.** In genome mode the estimator vector is
  cloned once per genome so each genome accumulates independently within a single
  BAM pass (see `genome.rs`); estimators must therefore implement `Clone`
  faithfully and reset via `setup()`.

- **Mate pairing by name in a BTreeMap.** `ReferenceSortedBamFilter` holds unpaired
  first-mates in `first_set` until the second mate is seen; this relies on both
  mates falling within the same reference block of a sorted BAM.

- **`add_coverage_entry` on the `CoverageTaker` trait is an acknowledged hack**
  (see comment) used only by `PileupCountsGenomeCoverageEstimator`; it is the one
  place the otherwise-uniform taker interface leaks a metric-specific concern.

- **`GenomesAndContigs` is copy-shared with cockatoo** (file header) — keep the
  two in sync rather than diverging its API casually.

## Entry Points

- **Add a coverage metric:** add a variant to `CoverageEstimator`
  (`mosdepth_genome_coverage_estimators.rs`), implement its `add_contig` /
  `calculate_coverage` / `column_headers` / `setup`, then wire it into
  `EstimatorsAndTaker::generate_from_clap` and the CLI method list in `cli.rs`.
- **Add/rename a CLI flag:** `src/cli.rs` for the definition, `src/bin/coverm.rs`
  for reading it.
- **Support a new mapper:** extend `MappingProgram` (`bam_generator.rs`),
  `parse_mapping_program` / `check_mapping_program_dependencies`
  (`coverm.rs`), `run_index_command` (`mapping_index_maintenance.rs`), and
  `external_command_checker.rs`.
- **Change output format:** `CoveragePrinter` (`coverage_printer.rs`); add a
  caching mode in `coverage_takers.rs` if all values are needed at once.
- **Adjust read filtering:** `filter.rs` (thresholds) or `FlagFilter` in
  `lib.rs` (flag gating).
- **Coverage from pre-made BAMs vs reads:** both start in `run_genome` /
  the `contig` arm of `main`, diverging at which `generate_*` factory in
  `bam_generator.rs` is called.
- **Per-gene coverage:** `genes.rs`.

## Out of Scope

Update this file when modules are added, removed, or restructured — not for
routine function-level changes. It documents shape, not implementation detail.
