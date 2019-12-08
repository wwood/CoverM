extern crate assert_cli;
extern crate tempfile;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;
    extern crate tempfile;
    use std;
    use std::io::Write;

    #[test]
    fn test_filter_all_reads(){
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        Assert::main_binary()
            .with_args(&[
                "filter",
                "-b",
                "tests/data/2seqs.bad_read.1.bam",
                "-o",
                t]).succeeds().unwrap();
        Assert::command(&["samtools","view",t])
            .stdout().contains("1\t99\tseq1").unwrap();
    }

    #[test]
    fn test_filter_filter_out(){
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        Assert::main_binary()
            .with_args(&[
                "filter",
                "--min-read-percent-identity-pair",
                "0.99",
                "-b",
                "tests/data/2seqs.bad_read.1.bam",
                "-o",
                t,
                "--proper-pairs-only"]).succeeds().unwrap();
        Assert::command(&["samtools","view",t])
            .stdout().doesnt_contain("1\t99\tseq1").unwrap();
    }

    #[test]
    fn test_contig_tempdir_index_creation(){
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t_full = tf.path().to_str().unwrap();
        std::fs::copy("tests/data/7seqs.fna", t_full).unwrap();
        let t = tf.path().file_name().unwrap().to_str().unwrap();
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--output-format",
                "sparse",
                "--min-read-percent-identity-pair",
                "0.95",
                "--contig-end-exclusion",
                "0",
                "-r",
                t_full,
                "-1",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "-2",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--proper-pairs-only",
            ])
            .succeeds()
            .stdout().contains(format!(
                "{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11000\t0
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11010\t0
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome2~seq1\t1.2
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome3~random_sequence_length_11001\t0
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome4~random_sequence_length_11002\t0
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome5~seq2\t1.2
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome6~random_sequence_length_11003\t0",
                t, t, t, t, t, t, t).as_str()).unwrap();
    }

    #[test]
    #[ignore] // known failure, cannot currently take multiple references
    fn test_contig_multiple_references(){
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--output-format",
                "sparse",
                "--min-read-percent-identity-pair",
                "0.95",
                "--contig-end-exclusion",
                "0",
                "-r",
                "tests/data/7seqs.fna",
                "tests/data/2seqs.fasta",
                "-1",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "-2",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
            ])
            .succeeds()
            .stdout().contains("7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11000\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11010\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome2~seq1\t1.2
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome3~random_sequence_length_11001\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome4~random_sequence_length_11002\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome5~seq2\t1.2
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome6~random_sequence_length_11003\t0
2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz\tseq1\t1.2
2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz\tseq2\t1.2").unwrap();
    }

    #[test]
    #[ignore] // cannot currently take multiple references
    fn test_coupled_reads_input(){
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--min-read-percent-identity",
                "0.95",
                "--contig-end-exclusion",
                "0",
                "-r",
                "tests/data/7seqs.fna",
                "tests/data/2seqs.fasta",
                "-1",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "-2",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "-c",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
            ])
            .succeeds()
            .stdout().contains("7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11000\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11010\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome2~seq1\t1.2
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome3~random_sequence_length_11001\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome4~random_sequence_length_11002\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome5~seq2\t1.2
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome6~random_sequence_length_11003\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11000\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11010\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome2~seq1\t1.2
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome3~random_sequence_length_11001\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome4~random_sequence_length_11002\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome5~seq2\t1.2
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome6~random_sequence_length_11003\t0
2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz\tseq1\t1.2
2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz\tseq2\t1.2
2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz\tseq1\t1.2
2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz\tseq2\t1.2").unwrap();
    }

    #[test]
    fn test_unfiltered_interleaved_input(){
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--contig-end-exclusion",
                "0",
                "-r",
                "tests/data/2seqs.fasta",
                "--output-format",
                "sparse",
                "--interleaved",
                "tests/data/bad_reads.interleaved.fq",
            ])
            .succeeds()
            .stdout().contains(
                "2seqs.fasta/bad_reads.interleaved.fq\tseq1\t0.899\n\
                 2seqs.fasta/bad_reads.interleaved.fq\tseq2\t0"
            ).unwrap();
    }

    #[test]
    fn test_filtered_interleaved_input(){
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-r",
                "tests/data/2seqs.fasta",
                "--interleaved",
                "tests/data/bad_reads.all.interleaved.fa",
                "--min-read-aligned-length-pair",
                "300",
                "--contig-end-exclusion",
                "0",
                "--proper-pairs-only",
                "--output-format",
                "sparse",
            ])
            .succeeds()
            .stdout().contains("2seqs.fasta/bad_reads.all.interleaved.fa\tseq1\t1.2
2seqs.fasta/bad_reads.all.interleaved.fa\tseq2\t1.5").unwrap();
    }

    #[test]
    fn test_single_reads_input(){
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--output-format",
                "sparse",
                "-r",
                "tests/data/2seqs.fasta",
                "--single",
                "tests/data/bad_reads.interleaved.fq",
                "--contig-end-exclusion",
                "0",
            ])
            .succeeds()
            .stdout().contains(
                "2seqs.fasta/bad_reads.interleaved.fq\tseq1\t0.899\n\
                 2seqs.fasta/bad_reads.interleaved.fq\tseq2\t0"
            ).unwrap();
    }

    #[test]
    fn test_genome_coupled_read_input_argparsing(){
        // Had trouble (because of a bug in clap? with this previously)
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--output-format",
                "sparse",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--reference",
                "tests/data/7seqs.fna",
                "-s","~"]).succeeds().unwrap();
    }

    #[test]
    fn test_contig_coupled_read_input_argparsing(){
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--reference",
                "tests/data/7seqs.fna"]).succeeds().unwrap();
    }

    #[test]
    fn test_cache_bam_files(){
        let td = tempfile::TempDir::new().unwrap();
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--output-format",
                "sparse",
                "--reference",
                "tests/data/7seqs.fna",
                "--bam-file-cache-directory",
                td.path().to_str().unwrap()
            ]).succeeds().stdout().contains(
                "Sample	Contig	Mean
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome1~random_sequence_length_11000	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome1~random_sequence_length_11010	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome2~seq1	1.4117647
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome3~random_sequence_length_11001	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome4~random_sequence_length_11002	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome5~seq2	1.2435294
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome6~random_sequence_length_11003	0").unwrap();
        assert!(td.path()
                .join("7seqs.fna.reads_for_seq1_and_seq2.1.fq.gz.bam")
                .is_file());
    }

    #[test]
    fn test_non_existant_cache_bam_files(){
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--reference",
                "tests/data/7seqs.fna",
                "--bam-file-cache-directory",
                "/no/no/no/165"
            ]).fails().unwrap();
    }

    #[test]
    fn test_unwriteable_cache_bam_files(){
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--reference",
                "tests/data/7seqs.fna",
                "--bam-file-cache-directory",
                "/"
            ]).fails().unwrap();
    }

    #[test]
    fn test_make(){
        let td = tempfile::TempDir::new().unwrap();
        Assert::main_binary()
            .with_args(&[
                "make",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--reference",
                "tests/data/7seqs.fna",
                "--output-directory",
                td.path().to_str().unwrap()
            ]).succeeds().unwrap();
        assert!(td.path()
                .join("7seqs.fna.reads_for_seq1_and_seq2.1.fq.gz.bam")
                .is_file());
    }

    #[test]
    fn test_make_with_mkdir(){
        let td = tempfile::TempDir::new().unwrap();
        Assert::main_binary()
            .with_args(&[
                "make",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--reference",
                "tests/data/7seqs.fna",
                "--output-directory",
                format!("{}/unmade_directory",td.path().to_str().unwrap()).as_str()
            ]).succeeds().unwrap();
        assert!(td.path()
                .join("unmade_directory")
                .join("7seqs.fna.reads_for_seq1_and_seq2.1.fq.gz.bam")
                .is_file());
    }

    #[test]
    fn test_relative_abundance_all_mapped() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--output-format",
                "sparse",
                "--reference",
                "tests/data/7seqs.fna",
                "-s","~"]).succeeds().stdout().contains(
                "Sample	Genome	Relative Abundance (%)
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	unmapped	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome1	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome2	53.16792
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome3	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome4	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome5	46.832077
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome6	0"
            ).unwrap();
    }

    #[test]
    fn test_relative_abundance_and_mean() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-m",
                "relative_abundance",
                "mean",
                "-b",
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
                "--output-format",
                "sparse",
                "-s",
                "~"]).succeeds().stdout().contains(
                "Sample	Genome	Relative Abundance (%)	Mean
7seqs.reads_for_seq1_and_seq2	unmapped	0	NA
7seqs.reads_for_seq1_and_seq2	genome1	0	0
7seqs.reads_for_seq1_and_seq2	genome2	53.16792	1.4117647
7seqs.reads_for_seq1_and_seq2	genome3	0	0
7seqs.reads_for_seq1_and_seq2	genome4	0	0
7seqs.reads_for_seq1_and_seq2	genome5	46.832077	1.2435294
7seqs.reads_for_seq1_and_seq2	genome6	0	0").unwrap();
    }

    #[test]
    fn test_contig_dense_output_simple() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-b",
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
                "--output-format",
                "dense"]).succeeds().stdout().contains(
                "Contig	7seqs.reads_for_seq1_and_seq2 Mean
genome1~random_sequence_length_11000	0
genome1~random_sequence_length_11010	0
genome2~seq1	1.4117647
genome3~random_sequence_length_11001	0
genome4~random_sequence_length_11002	0
genome5~seq2	1.2435294
genome6~random_sequence_length_11003	0").unwrap();
    }

    #[test]
    fn test_genome_dense_output_simple() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-m",
                "relative_abundance",
                "-b",
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
                "-s",
                "~",
                "--output-format",
                "dense"]).succeeds().stdout().contains(
                "Genome	7seqs.reads_for_seq1_and_seq2 Relative Abundance (%)
unmapped	0
genome1	0
genome2	53.167923
genome3	0
genome4	0
genome5	46.832077
genome6	0").unwrap();
    }

    #[test]
    fn test_genome_unknown_reason() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-m",
                "relative_abundance",
                "mean",
                "variance",
                "-r",
                "tests/data/7seqs.fna",
                "--output-format",
                "sparse",
                "-c",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--single",
                "tests/data/reads_for_seq1_and_seq2.fna",
                "-s",
                "~"]).succeeds().stdout().contains(
                "Sample	Genome	Relative Abundance (%)	Mean	Variance
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	unmapped	0	NA	NA
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome1	0	0	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome2	53.16792	1.4117647	1.3049262
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome3	0	0	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome4	0	0	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome5	46.832077	1.2435294	0.6862065
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome6	0	0	0
7seqs.fna/reads_for_seq1_and_seq2.fna	unmapped	0	NA	NA
7seqs.fna/reads_for_seq1_and_seq2.fna	genome1	0	0	0
7seqs.fna/reads_for_seq1_and_seq2.fna	genome2	53.16792	1.4117647	1.3049262
7seqs.fna/reads_for_seq1_and_seq2.fna	genome3	0	0	0
7seqs.fna/reads_for_seq1_and_seq2.fna	genome4	0	0	0
7seqs.fna/reads_for_seq1_and_seq2.fna	genome5	46.832077	1.2435294	0.6862065
7seqs.fna/reads_for_seq1_and_seq2.fna	genome6	0	0	0
").unwrap();
    }

    #[test]
    fn test_contig_unknown_reason() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--output-format",
                "sparse",
                "-m",
                "mean",
                "variance",
                "-r",
                "tests/data/7seqs.fna",
                "-c",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--single",
                "tests/data/reads_for_seq1_and_seq2.fna"]).succeeds().stdout().contains(
                "Sample	Contig	Mean	Variance
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome1~random_sequence_length_11000	0	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome1~random_sequence_length_11010	0	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome2~seq1	1.4117647	1.3049262
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome3~random_sequence_length_11001	0	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome4~random_sequence_length_11002	0	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome5~seq2	1.2435294	0.6862065
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome6~random_sequence_length_11003	0	0
7seqs.fna/reads_for_seq1_and_seq2.fna	genome1~random_sequence_length_11000	0	0
7seqs.fna/reads_for_seq1_and_seq2.fna	genome1~random_sequence_length_11010	0	0
7seqs.fna/reads_for_seq1_and_seq2.fna	genome2~seq1	1.4117647	1.3049262
7seqs.fna/reads_for_seq1_and_seq2.fna	genome3~random_sequence_length_11001	0	0
7seqs.fna/reads_for_seq1_and_seq2.fna	genome4~random_sequence_length_11002	0	0
7seqs.fna/reads_for_seq1_and_seq2.fna	genome5~seq2	1.2435294	0.6862065
7seqs.fna/reads_for_seq1_and_seq2.fna	genome6~random_sequence_length_11003	0	0
").unwrap();
    }

    #[test]
    fn test_bwa_parameters() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--output-format",
                "sparse",
                "-r",
                "tests/data/2seqs.fasta",
                "--single",
                "tests/data/2seqs.fasta"]).succeeds().stdout().contains(
                "Sample	Contig	Mean\n\
                 2seqs.fasta/2seqs.fasta	seq1	1\n\
                 2seqs.fasta/2seqs.fasta	seq2	1\n").unwrap();

        Assert::main_binary()
            .with_args(&[
                "contig",
                "--output-format",
                "sparse",
                "-p",
                "bwa-mem",
                "-r",
                "tests/data/2seqs.fasta",
                "--bwa-parameters",
                "'-k 5000'", // seed length longer than both sequences
                "--single",
                "tests/data/2seqs.fasta"]).succeeds().stdout().contains(
                "Sample	Contig	Mean\n\
                 2seqs.fasta/2seqs.fasta	seq1	0\n\
                 2seqs.fasta/2seqs.fasta	seq2	0\n").unwrap();
    }

    #[test]
    fn test_metabat_include_supplementary() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-m",
                "metabat",
                "-b",
                "tests/data/k141_7.reheadered.bam", // includes a supplementary alignment
            ]).succeeds().stdout().contains(
                "contigName	contigLen	totalAvgDepth	k141_7.reheadered.bam	k141_7.reheadered.bam-var
k141_7	350	0.69	0.69	2.0843215").unwrap();
    }


    #[test]
    fn test_metabat_97_of_100_bases_should_fail() {
        // 100 bases match, 3 mismatch. Seems metabat uses > 0.97 not >= 0.97
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-m",
                "metabat",
                "-b",
                "tests/data/k141_2005182.head11.bam", // includes a supplementary alignment
            ]).succeeds().stdout().contains(
                "contigName	contigLen	totalAvgDepth	k141_2005182.head11.bam	k141_2005182.head11.bam-var
k141_2005182	225	1.9333333	1.9333333	0.063063025").unwrap();
    }

    #[test]
    fn test_deletions_count_towards_perc_id() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-m",
                "metabat",
                "-b",
                "tests/data/k141_109815.stray_read.bam", // If D not counted, perc ID filters out
            ]).succeeds().stdout().contains(
        "contigName	contigLen	totalAvgDepth	k141_109815.stray_read.bam	k141_109815.stray_read.bam-var
k141_109815	362	0.6273585	0.6273585	0.23488776").unwrap();
    }

    #[test]
    fn test_no_zeroes_missing_column_bug_genomes_and_contigs() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--output-format",
                "sparse",
                "-m",
                "trimmed_mean",
                "mean",
                "covered_fraction",
                "--no-zeros",
                "-r",
                "tests/data/2seqs.fasta",
                "-c",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--genome-fasta-files",
                "tests/data/genomes_dir/seq1.fna",
                "tests/data/genomes_dir/seq2.fna",
                "--trim-max",
                "0.01",
                "--trim-min",
                "0.009"
            ]).succeeds().stdout().contains(
                "2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz	seq1	0	1.4117647	0.669
2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz	seq2	0	1.2435294	0.849"
            ).unwrap();
    }

    #[test]
    fn test_no_zeroes_missing_column_bug_separator() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--output-format",
                "sparse",
                "-m",
                "trimmed_mean",
                "mean",
                "covered_fraction",
                "--no-zeros",
                "-r",
                "tests/data/2seqs.fasta",
                "-c",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "-s",
                "q",
                "--trim-max",
                "0.01",
                "--trim-min",
                "0.009"
            ]).succeeds().stdout().contains(
                "2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz	se	0	1.3276471	0.759"
            ).unwrap();
    }

    #[test]
    fn test_autoconcatenation() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-m",
                "trimmed_mean",
                "mean",
                "covered_fraction",
                "--no-zeros",
                "-c",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--genome-fasta-files",
                "tests/data/genomes_dir/seq1.fna",
                "tests/data/genomes_dir/seq2.fna",
                "--output-format",
                "sparse",
                "--trim-max",
                "0.01",
                "--trim-min",
                "0.009"
            ]).succeeds().stdout().contains(
                "reads_for_seq1_and_seq2.1.fq.gz	seq1	0	1.4117647	0.669\n\
                 reads_for_seq1_and_seq2.1.fq.gz	seq2	0	1.2435294	0.849\n"
            ).unwrap();
    }

    #[test]
    fn test_single_genome_bug_dense() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--single-genome",
                "-r",
                "tests/data/dense_interleaved_single_genome_bug/ref.fna",
                "--interleaved",
                "tests/data/dense_interleaved_single_genome_bug/reads_interleaved.fna",
                "tests/data/dense_interleaved_single_genome_bug/reads_interleaved2.fna",
                "--output-format",
                "dense"
            ]).succeeds().stdout().is(
                "Genome\tref.fna/reads_interleaved.fna Relative Abundance (%)\t\
                 ref.fna/reads_interleaved2.fna Relative Abundance (%)\n\
                 unmapped\t25\t33.333332\n\
                 genome1\t75\t66.66667\n").unwrap();
    }

    #[test]
    fn test_filter_unmapped_not_inverse() {
        let tf1: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t1 = tf1.path().to_str().unwrap();
        let tf2: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t2 = tf2.path().to_str().unwrap();
        Assert::main_binary()
            .with_args(&[
                "filter",
                "--min-read-aligned-length",
                "1",
                "-b",
                "tests/data/dense_interleaved_single_genome_bug/ref.fna.reads_interleaved.fna.bam",
                "tests/data/dense_interleaved_single_genome_bug/ref.fna.reads_interleaved2.fna.bam",
                "-o",
                t1,
                t2]).succeeds().unwrap();
        Assert::command(&["samtools","view",t1])
            .stdout().contains("random_sequence_length_1000")
            .stdout().doesnt_contain("random_sequence_length_100\t").unwrap();
        Assert::command(&["samtools","view",t2])
            .stdout().contains("random_sequence_length_1000")
            .stdout().doesnt_contain("random_sequence_length_100\t").unwrap();
    }

    #[test]
    fn test_filter_unmapped_inverse(){
        let tf1: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t1 = tf1.path().to_str().unwrap();
        Assert::main_binary()
            .with_args(&[
                "filter",
                "--inverse",
                "-b",
                "tests/data/dense_interleaved_single_genome_bug/ref.fna.r1.fna.bam",
                "-o",
                t1]).succeeds().unwrap();
        Assert::command(&["samtools","view",t1])
            .stdout().doesnt_contain("random_sequence_length_1000")
            .stdout().contains("seq4\t77")
            .stdout().contains("seq4\t141")
            .unwrap();
    }

    #[test]
    fn test_filter_unmapped_inverse_improper_pairs(){
        let tf1: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t1 = tf1.path().to_str().unwrap();
        let tf2: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t2 = tf2.path().to_str().unwrap();
        Assert::main_binary()
            .with_args(&[
                "filter",
                "--inverse",
                "-b",
                // all mappings are improper pairs since too few reads were aligned
                "tests/data/dense_interleaved_single_genome_bug/ref.fna.reads_interleaved.fna.bam",
                "tests/data/dense_interleaved_single_genome_bug/ref.fna.reads_interleaved2.fna.bam",
                "-o",
                t1,
                t2]).succeeds().unwrap();
        Assert::command(&["samtools","view",t1])
            .stdout().contains("random_sequence_length_1000")
            .stdout().contains("random_sequence_length_100\t77")
            .stdout().contains("random_sequence_length_100\t141")
            .unwrap();
        Assert::command(&["samtools","view",t2])
            .stdout().contains("random_sequence_length_1000")
            .stdout().contains("random_sequence_length_100\t77")
            .stdout().contains("random_sequence_length_100\t141")
            .unwrap();
     }

    #[test]
    fn test_caches_when_reference_not_specified() {
        let td = tempfile::TempDir::new().unwrap();
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--genome-fasta-directory",
                "tests/data/genomes_dir/",
                "--bam-file-cache-directory",
                td.path().to_str().unwrap()
            ]).succeeds().unwrap();
        assert!(td.path()
                .join("coverm-genome.reads_for_seq1_and_seq2.1.fq.gz.bam")
                .is_file());
    }

    #[test]
    fn test_sharding_no_exclusion_genome_separator() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--sharded",
                "-b",
                "tests/data/shard1.bam",
                "tests/data/shard2.bam",
                "-s",
                "~",
            ])
            .stdout().is("Genome	shard1|shard2 Relative Abundance (%)
unmapped	0
genome3	25.024881
genome4	25.022575
genome5	0
genome6	25.020271
genome1	24.932274
genome2	0
")
            .succeeds().unwrap()
    }

    #[test]
    fn test_sharding_no_exclusion_contig() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--sharded",
                "-b",
                "tests/data/shard1.bam",
                "tests/data/shard2.bam",
            ])
            .stdout().is("Contig	shard1|shard2 Mean
genome3~random_sequence_length_11001	0.110588886
genome4~random_sequence_length_11002	0.11057869
genome5~seq2	0
genome6~random_sequence_length_11003	0.11056851
genome1~random_sequence_length_11000	0.109861754
genome1~random_sequence_length_11010	0.110497236
genome2~seq1	0
")
            .succeeds().unwrap()
    }

    #[test]
    fn test_sharding_no_exclusion_bwa_contig() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-p",
                "bwa-mem",
                "--sharded",
                "-b",
                "tests/data/shard1.bam",
                "tests/data/shard2.bam",
            ])
            .stdout().is("Contig	shard1|shard2 Mean
genome3~random_sequence_length_11001	0.110588886
genome4~random_sequence_length_11002	0.11057869
genome5~seq2	0
genome6~random_sequence_length_11003	0.11056851
genome1~random_sequence_length_11000	0.109861754
genome1~random_sequence_length_11010	0.110497236
genome2~seq1	0
")
            .succeeds().unwrap()
    }

    #[test]
    fn test_sharding_exclusion_genome_separator() {
        let mut tf1: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        writeln!(tf1, "genome3").unwrap();
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--sharded",
                "-b",
                "tests/data/shard1.bam",
                "tests/data/shard2.bam",
                "-s",
                "~",
                "--exclude-genomes-from-deshard",
                tf1.path().to_str().unwrap()
            ])
            .stdout().is("Genome	shard1|shard2 Relative Abundance (%)
unmapped	19.999998
genome3	0
genome4	26.699606
genome5	0
genome6	26.697144
genome1	26.60325
genome2	0
")
            .succeeds().unwrap()
    }

    #[test]
    fn test_sharding_exclusion_genomes_fasta_files_definition() {
        let mut tf1: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        writeln!(tf1, "genome3").unwrap();
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--sharded",
                "-b",
                "tests/data/shard1.bam",
                "tests/data/shard2.bam",
                "--genome-fasta-files",
                "tests/data/genomes_dir_7seqs/genome1.fasta",
                "tests/data/genomes_dir_7seqs/genome2.fasta",
                "tests/data/genomes_dir_7seqs/genome3.fasta",
                "tests/data/genomes_dir_7seqs/genome4.fasta",
                "tests/data/genomes_dir_7seqs/genome5.fasta",
                "tests/data/genomes_dir_7seqs/genome6.fasta",
                "--exclude-genomes-from-deshard",
                tf1.path().to_str().unwrap()
            ])
            .stdout().is("Genome	shard1|shard2 Relative Abundance (%)
unmapped	19.999998
genome1	26.60325
genome2	0
genome3	0
genome4	26.699606
genome5	0
genome6	26.697144
")
            .succeeds().unwrap()
    }

    #[test]
    fn test_correct_number_of_reads_total_with_filtering_paired() {
        let mut tf1: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        writeln!(tf1, "genome3").unwrap();
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-c",
                "tests/data/7seqs.reads_for_7_plus5random.1.fa",
                "tests/data/7seqs.reads_for_7_plus5random.2.fa",
                "-r",
                "tests/data/7seqs.fna",
                "--proper-pairs-only",
                "--min-read-aligned-length-pair",
                "50"])
            .stderr().contains(
                "coverm::contig] In sample '7seqs.fna/7seqs.reads_for_7_plus5random.1.fa', \
                 found 40 reads mapped out of 50 total (80.00%)")
            .succeeds().unwrap()
    }

    #[test]
    fn test_correct_number_of_reads_total_with_filtering_single() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-c",
                "tests/data/7seqs.reads_for_7_plus5random.1.fa",
                "tests/data/7seqs.reads_for_7_plus5random.2.fa",
                "-r",
                "tests/data/7seqs.fna",
                "--min-read-aligned-length",
                "50"])
            .stderr().contains(
                "coverm::contig] In sample '7seqs.fna/7seqs.reads_for_7_plus5random.1.fa', \
                 found 40 reads mapped out of 50 total (80.00%)")
            .succeeds().unwrap()
    }

    #[test]
    fn test_sharded_contig_input_reads() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-c",
                "tests/data/7seqs.reads_for_7.1.fq",
                "tests/data/7seqs.reads_for_7.2.fq",
                "-r",
                "tests/data/shard1.fna",
                "tests/data/shard2.fna",
                "--sharded",
            ])
            .stdout().is(
                "Contig	shard1.fna|shard2.fna/7seqs.reads_for_7.1.fq|7seqs.reads_for_7.1.fq Mean\n\
                 genome3~random_sequence_length_11001	0.110588886\n\
                 genome4~random_sequence_length_11002	0.11057869\n\
                 genome5~seq2	0\n\
                 genome6~random_sequence_length_11003	0.11056851\n\
                 genome1~random_sequence_length_11000	0.109861754\n\
                 genome1~random_sequence_length_11010	0.110497236\n\
                 genome2~seq1	0\n")
            .succeeds().unwrap()
    }

    #[test]
    fn test_genome_definition() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--genome-definition",
                "tests/data/7seqs.definition",
                "-b",
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam"])
            .succeeds()
            .stdout().contains(
                "Genome	7seqs.reads_for_seq1_and_seq2 Relative Abundance (%)\n")
            .stdout().contains("genome2	53.167923\n")
            .stdout().contains("genome5	46.832077\n")
            .unwrap();
    }

    #[test]
    fn test_some_samples_zero_coverage_genome() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-s",
                "~",
                "-r",
                "tests/data/7seqs.fna",
                "-1",
                "tests/data/7seqs.reads_for_7.1.fq",
                "tests/data/random.fq",
                "-2",
                "tests/data/7seqs.reads_for_7.2.fq",
                "tests/data/random.fq"])
            .succeeds()
            .stdout().is("Genome	7seqs.fna/7seqs.reads_for_7.1.fq Relative Abundance (%)\t\
                          7seqs.fna/random.fq Relative Abundance (%)\n\
                          unmapped	0	100\n\
                          genome1	24.932272	NaN\n\
                          genome2	0	NaN\n\
                          genome3	25.02488	NaN\n\
                          genome4	25.022572	NaN\n\
                          genome5	0	NaN\n\
                          genome6	25.02027	NaN\n")
            .unwrap();
    }

    #[test]
    fn test_make_with_bwa() {
        let td = tempfile::TempDir::new().unwrap();
        Assert::main_binary()
            .with_args(&[
                "make",
                "--mapper",
                "bwa-mem",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--reference",
                "tests/data/7seqs.fna",
                "--output-directory",
                td.path().to_str().unwrap()
            ]).succeeds().unwrap();
        let bam = td.path()
                .join("7seqs.fna.reads_for_seq1_and_seq2.1.fq.gz.bam");
        assert!(bam.is_file());
        Assert::command(&["samtools","view","-H",bam.to_str().unwrap()])
            .stdout().contains("PN:bwa").unwrap();
    }

    #[test]
    fn test_make_with_minimap2() {
        let td = tempfile::TempDir::new().unwrap();
        Assert::main_binary()
            .with_args(&[
                "make",
                "--mapper",
                "minimap2-sr",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--reference",
                "tests/data/7seqs.fna",
                "--output-directory",
                td.path().to_str().unwrap()
            ]).succeeds().unwrap();
        let bam = td.path()
                .join("7seqs.fna.reads_for_seq1_and_seq2.1.fq.gz.bam");
        assert!(bam.is_file());
        Assert::command(&["samtools","view","-H",bam.to_str().unwrap()])
            .stdout().contains("PN:minimap2").unwrap();
    }

    #[test]
    fn test_contig_sparse_rpkm() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-m",
                "rpkm",
                "reads_per_base",
                "-b",
                "tests/data/7seqs.fnaVbad_read.bam",
                "-m",
                "length",
                "count",
                "--output-format",
                "sparse"
            ])
            .succeeds()
            .stdout().is("Sample	Contig	RPKM	Reads per base	Length	Read Count\n\
                7seqs.fnaVbad_read	genome1~random_sequence_length_11000	0	0	11000	0\n\
                7seqs.fnaVbad_read	genome1~random_sequence_length_11010	0	0	11010	0\n\
                7seqs.fnaVbad_read	genome2~seq1	500000	0.01	1000	10\n\
                7seqs.fnaVbad_read	genome3~random_sequence_length_11001	0	0	11001	0\n\
                7seqs.fnaVbad_read	genome4~random_sequence_length_11002	0	0	11002	0\n\
                7seqs.fnaVbad_read	genome5~seq2	500000	0.01	1000	10\n\
                7seqs.fnaVbad_read	genome6~random_sequence_length_11003	0	0	11003	0\n")
            .unwrap();
    }

    #[test]
    fn test_contig_dense_rpkm() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-m",
                "rpkm",
                "reads_per_base",
                "length",
                "count",
                "-b",
                "tests/data/7seqs.fnaVbad_read.bam",
            ])
            .succeeds()
            .stdout().is(
                "Contig	7seqs.fnaVbad_read RPKM	7seqs.fnaVbad_read Reads per base	7seqs.fnaVbad_read Length	7seqs.fnaVbad_read Read Count\n\
                genome1~random_sequence_length_11000	0	0	11000	0\n\
                genome1~random_sequence_length_11010	0	0	11010	0\n\
                genome2~seq1	500000	0.01	1000	10\n\
                genome3~random_sequence_length_11001	0	0	11001	0\n\
                genome4~random_sequence_length_11002	0	0	11002	0\n\
                genome5~seq2	500000	0.01	1000	10\n\
                genome6~random_sequence_length_11003	0	0	11003	0\n")
            .unwrap();
    }

    #[test]
    fn test_single_genome_dense_rpkm() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--single-genome",
                "-m",
                "rpkm",
                "reads_per_base",
                "length",
                "count",
                "--min-covered-fraction",
                "0",
                "-b",
                "tests/data/7seqs.fnaVbad_read.bam",
            ])
            .succeeds()
            .stdout().is(
                "Genome	7seqs.fnaVbad_read RPKM	7seqs.fnaVbad_read Reads per base	7seqs.fnaVbad_read Length	7seqs.fnaVbad_read Read Count\n\
                genome1	17538.936	0.00035077872	57016	20\n")
            .unwrap();
    }

    #[test]
    fn test_single_genome_rpkm_min_covered_fraction() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--single-genome",
                "-m",
                "rpkm",
                "-b",
                "tests/data/7seqs.fnaVbad_read.bam",
            ])
            .succeeds()
            .stdout().is(
                "Genome	7seqs.fnaVbad_read RPKM\n\
                genome1	0\n")
            .unwrap();
    }

    #[test]
    fn test_ont_single_sample() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-m",
                "rpkm",
                "mean",
                "count",
                "-p",
                "minimap2-ont",
                "--single",
                "tests/data/ont.reads.fq.gz",
                "-r",
                "tests/data/ont.ref.fna",
            ])
            .succeeds()
            .stdout().is(
                "Contig	ont.ref.fna/ont.reads.fq.gz RPKM	ont.ref.fna/ont.reads.fq.gz Mean	ont.ref.fna/ont.reads.fq.gz Read Count\n\
                ctg4	1887.6234	0.02457587	6\n\
                ctg5	717.27264	0.0041760243	3\n\
                ctg6	126.346375	0.0021053297	1\n")
            .unwrap();
    }

    #[test]
    fn test_ont_two_samples() {
        // Also tests that -t is being set when minimap2 indexing
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-m",
                "rpkm",
                "mean",
                "count",
                "-p",
                "minimap2-ont",
                "--single",
                "tests/data/ont.reads.fq.gz",
                "tests/data/ont.reads.fq.gz",
                "-r",
                "tests/data/ont.ref.fna",
                "-t",
                "2",
                "-v",
            ])
            .succeeds()
            .stdout()
            .is(
                "Contig	ont.ref.fna/ont.reads.fq.gz RPKM	ont.ref.fna/ont.reads.fq.gz Mean	ont.ref.fna/ont.reads.fq.gz Read Count	ont.ref.fna/ont.reads.fq.gz RPKM	ont.ref.fna/ont.reads.fq.gz Mean	ont.ref.fna/ont.reads.fq.gz Read Count\n\
                ctg4	1887.6234	0.02457587	6	1887.6234	0.02457587	6\n\
                ctg5	717.27264	0.0041760243	3	717.27264	0.0041760243	3\n\
                ctg6	126.346375	0.0021053297	1	126.346375	0.0021053297	1\n")
            .stderr()
            .contains(
                    "Running DB indexing command: \"minimap2\" \"-x\" \"map-ont\" \"-t\" \"2\" \"-d\"")
            .unwrap();
    }

    #[test]
    fn test_pb_two_samples() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-m",
                "rpkm",
                "mean",
                "count",
                "-p",
                "minimap2-pb",
                "--single",
                "tests/data/ont.reads.fq.gz",
                "tests/data/ont.reads.fq.gz",
                "-r",
                "tests/data/ont.ref.fna",
                "-v",
            ])
            .succeeds()
            .stdout()
            .is(
                "Contig	ont.ref.fna/ont.reads.fq.gz RPKM	ont.ref.fna/ont.reads.fq.gz Mean	ont.ref.fna/ont.reads.fq.gz Read Count	ont.ref.fna/ont.reads.fq.gz RPKM	ont.ref.fna/ont.reads.fq.gz Mean	ont.ref.fna/ont.reads.fq.gz Read Count\n\
                ctg4	1797.7366	0.02386453	4	1797.7366	0.02386453	4\n\
                ctg5	683.11676	0.0035182887	2	683.11676	0.0035182887	2\n\
                ctg6	180.49483	0.0021053297	1	180.49483	0.0021053297	1\n")
            .stderr()
            .contains(
                    "Running DB indexing command: \"minimap2\" \"-x\" \"map-pb\" \"-t\" \"1\" \"-d\"")
            .unwrap();
    }

    #[test]
    fn test_minimap2_no_preset_with_params() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-m",
                "rpkm",
                "mean",
                "count",
                "-p",
                "minimap2-no-preset",
                "--single",
                "tests/data/ont.reads.fq.gz",
                "-r",
                "tests/data/ont.ref.fna",
                "-v",
                "--minimap2-parameters",
                "-A 20"
            ])
            .succeeds()
            .stdout()
            .is(
                "Contig	ont.ref.fna/ont.reads.fq.gz RPKM	ont.ref.fna/ont.reads.fq.gz Mean	ont.ref.fna/ont.reads.fq.gz Read Count\n\
                ctg4	2202.2273	0.081332035	7\n\
                ctg5	478.18173	0.010246328	2\n\
                ctg6	126.346375	0.0021596688	1\n")
            .stderr()
            .contains("-A 20 -t 1")
            .unwrap();
    }

    #[test]
    fn test_pregenerated_minimap2_index() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-m",
                "rpkm",
                "mean",
                "count",
                "--coupled",
                "tests/data/bad_read.1.fq",
                "tests/data/bad_read.2.fq",
                "-r",
                "tests/data/7seqs.fna.mmi",
                "-v",
                "--minimap2-reference-is-index",
            ])
            .succeeds()
            .stdout()
            .is(
                "Contig	7seqs.fna.mmi/bad_read.1.fq RPKM	7seqs.fna.mmi/bad_read.1.fq Mean	7seqs.fna.mmi/bad_read.1.fq Read Count
genome1~random_sequence_length_11000	0	0	0
genome1~random_sequence_length_11010	0	0	0
genome2~seq1	500000	1.6764706	10
genome3~random_sequence_length_11001	0	0	0
genome4~random_sequence_length_11002	0	0	0
genome5~seq2	500000	1.6764706	10
genome6~random_sequence_length_11003	0	0	0
")
            .stderr()
            .contains(
                    "Minimap2 uses mapping parameters defined when the index was created, not parameters defined when mapping")
            .unwrap();
    }

    #[test]
    fn test_genome_narrowing_rochelles_bug() {
        // This was an issue where samtools sort was running out of memory, so
        // starting to read the BAM failed, but coverm failed ungracefully.
        // Fixed in 4d02762a8b9f830c0b081c86b44a642ce9617c48. Test commented out
        // for posterity.

        // Assert::main_binary()
        //     .with_args(&[
        //         "genome",
        //         "-v",
        //         "-t",
        //         "2000", // many threads means sort asks for too much mem
        //         "--coupled",
        //         "tests/data/bad_read.1.fq",
        //         "tests/data/bad_read.2.fq",
        //         "tests/data/bad_read.1.fq",
        //         "tests/data/bad_read.2.fq",
        //         "-s",
        //         "q",
        //         "--min-read-aligned-length",
        //         "70",
        //         "--min-read-percent-identity",
        //         "0.97",
        //         "--min-covered-fraction",
        //         "0",
        //         "-m",
        //         "count",
        //         "--bam-file-cache-directory",
        //         "/tmp/bam_files",
        //         "-r",
        //         "tests/data/2seqs.fasta",
        //     ]).succeeds().unwrap();
    }
}


// TODO: Add mismatching bases test
// TODO: Filter fails when reference sequences are duplicated?
// TODO: Filter should spit things out if no thresholds are specified.
