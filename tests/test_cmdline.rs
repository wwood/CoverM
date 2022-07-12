extern crate assert_cli;
extern crate tempfile;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;
    extern crate tempfile;
    use std;
    use std::io::Read;
    use std::io::Write;
    use std::str;

    fn assert_equal_table(expected: &str, observed: &str) -> bool {
        // assert the first lines are the same
        let mut expected_lines = expected.lines();
        let mut observed_lines = observed.lines();
        assert_eq!(expected_lines.next(), observed_lines.next());

        // assert the rest of the lines are equal after sorting
        let mut expected_contents: Vec<_> = expected_lines.collect();
        let mut observed_contents: Vec<_> = observed_lines.collect();
        expected_contents.sort();
        observed_contents.sort();
        assert_eq!(expected_contents.join("\n"), observed_contents.join("\n"));

        true
    }

    #[test]
    fn test_filter_all_reads() {
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        Assert::main_binary()
            .with_args(&["filter", "-b", "tests/data/2seqs.bad_read.1.bam", "-o", t])
            .succeeds()
            .unwrap();
        Assert::command(&["samtools", "view", t])
            .stdout()
            .contains("1\t99\tseq1")
            .unwrap();
    }

    #[test]
    fn test_filter_filter_out() {
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
                "--proper-pairs-only",
            ])
            .succeeds()
            .unwrap();
        Assert::command(&["samtools", "view", t])
            .stdout()
            .doesnt_contain("1\t99\tseq1")
            .unwrap();
    }

    #[test]
    fn test_contig_tempdir_index_creation() {
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
            .stdout()
            .contains(
                format!(
                    "{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11000\t0
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11010\t0
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome2~seq1\t1.2
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome3~random_sequence_length_11001\t0
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome4~random_sequence_length_11002\t0
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome5~seq2\t1.2
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome6~random_sequence_length_11003\t0",
                    t, t, t, t, t, t, t
                )
                .as_str(),
            )
            .unwrap();
    }

    #[test]
    #[ignore] // known failure, cannot currently take multiple references
    fn test_contig_multiple_references() {
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
            .stdout()
            .contains(
                "7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11000\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11010\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome2~seq1\t1.2
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome3~random_sequence_length_11001\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome4~random_sequence_length_11002\t0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome5~seq2\t1.2
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome6~random_sequence_length_11003\t0
2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz\tseq1\t1.2
2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz\tseq2\t1.2",
            )
            .unwrap();
    }

    #[test]
    #[ignore] // cannot currently take multiple references
    fn test_coupled_reads_input() {
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
            .stdout()
            .contains(
                "7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11000\t0
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
2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz\tseq2\t1.2",
            )
            .unwrap();
    }

    #[test]
    fn test_unfiltered_interleaved_input() {
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
            .stdout()
            .contains(
                "2seqs.fasta/bad_reads.interleaved.fq\tseq1\t0.899\n\
                 2seqs.fasta/bad_reads.interleaved.fq\tseq2\t0",
            )
            .unwrap();
    }

    #[test]
    fn test_filtered_interleaved_input() {
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
            .stdout()
            .contains(
                "2seqs.fasta/bad_reads.all.interleaved.fa\tseq1\t1.2
2seqs.fasta/bad_reads.all.interleaved.fa\tseq2\t1.5",
            )
            .unwrap();
    }

    #[test]
    fn test_single_reads_input() {
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
            .stdout()
            .contains(
                "2seqs.fasta/bad_reads.interleaved.fq\tseq1\t0.899\n\
                 2seqs.fasta/bad_reads.interleaved.fq\tseq2\t0",
            )
            .unwrap();
    }

    #[test]
    fn test_genome_coupled_read_input_argparsing() {
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
                "-s",
                "~",
            ])
            .succeeds()
            .unwrap();
    }

    #[test]
    fn test_contig_coupled_read_input_argparsing() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--reference",
                "tests/data/7seqs.fna",
            ])
            .succeeds()
            .unwrap();
    }

    #[test]
    fn test_cache_bam_files() {
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
                td.path().to_str().unwrap(),
            ])
            .succeeds()
            .stdout()
            .contains(
                "Sample	Contig	Mean
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome1~random_sequence_length_11000	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome1~random_sequence_length_11010	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome2~seq1	1.4117647
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome3~random_sequence_length_11001	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome4~random_sequence_length_11002	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome5~seq2	1.2435294
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome6~random_sequence_length_11003	0",
            )
            .unwrap();
        assert!(td
            .path()
            .join("7seqs.fna.reads_for_seq1_and_seq2.1.fq.gz.bam")
            .is_file());
    }

    #[test]
    fn test_non_existant_cache_bam_files() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--reference",
                "tests/data/7seqs.fna",
                "--bam-file-cache-directory",
                "/no/no/no/165",
            ])
            .fails()
            .unwrap();
    }

    #[test]
    fn test_unwriteable_cache_bam_files() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--reference",
                "tests/data/7seqs.fna",
                "--bam-file-cache-directory",
                "/",
            ])
            .fails()
            .unwrap();
    }

    #[test]
    fn test_make() {
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
                td.path().to_str().unwrap(),
            ])
            .succeeds()
            .unwrap();
        assert!(td
            .path()
            .join("7seqs.fna.reads_for_seq1_and_seq2.1.fq.gz.bam")
            .is_file());
    }

    #[test]
    fn test_make_with_mkdir() {
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
                format!("{}/unmade_directory", td.path().to_str().unwrap()).as_str(),
            ])
            .succeeds()
            .unwrap();
        assert!(td
            .path()
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
                "-s",
                "~",
            ])
            .succeeds()
            .stdout()
            .contains(
                "Sample	Genome	Relative Abundance (%)
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	unmapped	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome1	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome2	53.16792
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome3	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome4	0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome5	46.832077
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz	genome6	0",
            )
            .unwrap();
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
                "~",
            ])
            .succeeds()
            .stdout()
            .contains(
                "Sample	Genome	Relative Abundance (%)	Mean
7seqs.reads_for_seq1_and_seq2	unmapped	0	NA
7seqs.reads_for_seq1_and_seq2	genome1	0	0
7seqs.reads_for_seq1_and_seq2	genome2	53.16792	1.4117647
7seqs.reads_for_seq1_and_seq2	genome3	0	0
7seqs.reads_for_seq1_and_seq2	genome4	0	0
7seqs.reads_for_seq1_and_seq2	genome5	46.832077	1.2435294
7seqs.reads_for_seq1_and_seq2	genome6	0	0",
            )
            .unwrap();
    }

    #[test]
    fn test_contig_dense_output_simple() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-b",
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
                "--output-format",
                "dense",
            ])
            .succeeds()
            .stdout()
            .contains(
                "Contig	7seqs.reads_for_seq1_and_seq2 Mean
genome1~random_sequence_length_11000	0
genome1~random_sequence_length_11010	0
genome2~seq1	1.4117647
genome3~random_sequence_length_11001	0
genome4~random_sequence_length_11002	0
genome5~seq2	1.2435294
genome6~random_sequence_length_11003	0",
            )
            .unwrap();
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
                "dense",
            ])
            .succeeds()
            .stdout()
            .contains(
                "Genome	7seqs.reads_for_seq1_and_seq2 Relative Abundance (%)
unmapped	0
genome1	0
genome2	53.167923
genome3	0
genome4	0
genome5	46.832077
genome6	0",
            )
            .unwrap();
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
                "~",
            ])
            .succeeds()
            .stdout()
            .contains(
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
",
            )
            .unwrap();
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
                "tests/data/reads_for_seq1_and_seq2.fna",
            ])
            .succeeds()
            .stdout()
            .contains(
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
",
            )
            .unwrap();
    }

    #[test]
    fn test_bwa_parameters() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--output-format",
                "sparse",
                "-p",
                "bwa-mem",
                "-r",
                "tests/data/2seqs.fasta",
                "--single",
                "tests/data/2seqs.fasta",
            ])
            .succeeds()
            .stdout()
            .contains(
                "Sample	Contig	Mean\n\
                 2seqs.fasta/2seqs.fasta	seq1	1\n\
                 2seqs.fasta/2seqs.fasta	seq2	1\n",
            )
            .unwrap();

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
                "tests/data/2seqs.fasta",
            ])
            .succeeds()
            .stdout()
            .contains(
                "Sample	Contig	Mean\n\
                 2seqs.fasta/2seqs.fasta	seq1	0\n\
                 2seqs.fasta/2seqs.fasta	seq2	0\n",
            )
            .unwrap();
    }

    #[test]
    fn test_bwa_mem2_parameters() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--output-format",
                "sparse",
                "-p",
                "bwa-mem2",
                "-r",
                "tests/data/2seqs.fasta",
                "--single",
                "tests/data/2seqs.fasta",
            ])
            .succeeds()
            .stdout()
            .contains(
                "Sample	Contig	Mean\n\
                 2seqs.fasta/2seqs.fasta	seq1	1\n\
                 2seqs.fasta/2seqs.fasta	seq2	1\n",
            )
            .unwrap();

        Assert::main_binary()
            .with_args(&[
                "contig",
                "--output-format",
                "sparse",
                "-p",
                "bwa-mem2",
                "-r",
                "tests/data/2seqs.fasta",
                "--bwa-parameters",
                "'-k 5000'", // seed length longer than both sequences
                "--single",
                "tests/data/2seqs.fasta",
            ])
            .succeeds()
            .stdout()
            .contains(
                "Sample	Contig	Mean\n\
                 2seqs.fasta/2seqs.fasta	seq1	0\n\
                 2seqs.fasta/2seqs.fasta	seq2	0\n",
            )
            .unwrap();
    }

    #[test]
    fn test_bwa_mem2_genome() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--output-format",
                "sparse",
                "-p",
                "bwa-mem2",
                "-r",
                "tests/data/2seqs.fasta",
                "--single",
                "tests/data/2seqs.fasta",
                "--single-genome",
            ])
            .succeeds()
            .stdout()
            .contains(
                "Sample	Genome	Relative Abundance (%)\n\
                 2seqs.fasta/2seqs.fasta	unmapped	0\n\
                 2seqs.fasta/2seqs.fasta	genome1	100\n",
            )
            .unwrap();
    }

    #[test]
    fn test_bwa_mem2_prefix_no_original() {
        // https://github.com/wwood/CoverM/issues/112
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--output-format",
                "sparse",
                "-p",
                "bwa-mem2",
                "-r",
                "tests/data/bwa_ref_without_original/2seqs.fasta",
                "--single",
                "tests/data/2seqs.fasta",
            ])
            .succeeds()
            .stdout()
            .contains(
                "Sample	Contig	Mean\n\
                 2seqs.fasta/2seqs.fasta	seq1	1\n\
                 2seqs.fasta/2seqs.fasta	seq2	1\n",
            )
            .unwrap();
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
            ])
            .succeeds()
            .stdout()
            .contains(
                "contigName	contigLen	totalAvgDepth	k141_7.reheadered.bam	k141_7.reheadered.bam-var
k141_7	350	0.69	0.69	2.0843",
            )
            .unwrap();
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
            ])
            .succeeds()
            .stdout()
            .contains(
                "contigName	contigLen	totalAvgDepth	k141_2005182.head11.bam	k141_2005182.head11.bam-var
k141_2005182	225	1.9333	1.9333	0.0631",
            )
            .unwrap();
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
k141_109815	362	0.6274	0.6274	0.2349").unwrap();
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
                "0.009",
            ])
            .succeeds()
            .stdout()
            .contains(
                "2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz	seq1	0	1.4117647	0.669
2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz	seq2	0	1.2435294	0.849",
            )
            .unwrap();
    }

    #[test]
    fn test_compressed_input_genome_fasta_files_with_reference() {
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
                "tests/data/genomes_dir_compressed/seq1.fna.gz",
                "tests/data/genomes_dir_compressed/seq2.fna.gz",
                "--trim-max",
                "0.01",
                "--trim-min",
                "0.009",
            ])
            .succeeds()
            .stdout()
            .contains(
                "2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz	seq1	0	1.4117647	0.669
2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz	seq2	0	1.2435294	0.849",
            )
            .unwrap();
    }

    #[test]
    fn test_compressed_input_genome_fasta_files_no_reference() {
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
                "-c",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--genome-fasta-files",
                "tests/data/genomes_dir_compressed/seq1.fna.gz",
                "tests/data/genomes_dir_compressed/seq2.fna.gz",
                "--trim-max",
                "0.01",
                "--trim-min",
                "0.009",
            ])
            .succeeds()
            .stdout()
            .contains(
                "reads_for_seq1_and_seq2.1.fq.gz	seq1	0	1.4117647	0.669
reads_for_seq1_and_seq2.1.fq.gz	seq2	0	1.2435294	0.849",
            )
            .unwrap();
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
                "0.009",
            ])
            .succeeds()
            .stdout()
            .contains("2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz	se	0	1.3276471	0.759")
            .unwrap();
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
                "0.009",
            ])
            .succeeds()
            .stdout()
            .contains(
                "reads_for_seq1_and_seq2.1.fq.gz	seq1	0	1.4117647	0.669\n\
                 reads_for_seq1_and_seq2.1.fq.gz	seq2	0	1.2435294	0.849\n",
            )
            .unwrap();
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
                "dense",
            ])
            .succeeds()
            .stdout()
            .is(
                "Genome\tref.fna/reads_interleaved.fna Relative Abundance (%)\t\
                 ref.fna/reads_interleaved2.fna Relative Abundance (%)\n\
                 unmapped\t25\t33.333332\n\
                 genome1\t75\t66.66667\n",
            )
            .unwrap();
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
                t2,
            ])
            .succeeds()
            .unwrap();
        Assert::command(&["samtools", "view", t1])
            .stdout()
            .contains("random_sequence_length_1000")
            .stdout()
            .doesnt_contain("random_sequence_length_100\t")
            .unwrap();
        Assert::command(&["samtools", "view", t2])
            .stdout()
            .contains("random_sequence_length_1000")
            .stdout()
            .doesnt_contain("random_sequence_length_100\t")
            .unwrap();
    }

    #[test]
    fn test_filter_unmapped_inverse() {
        let tf1: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t1 = tf1.path().to_str().unwrap();
        Assert::main_binary()
            .with_args(&[
                "filter",
                "--inverse",
                "-b",
                "tests/data/dense_interleaved_single_genome_bug/ref.fna.r1.fna.bam",
                "-o",
                t1,
            ])
            .succeeds()
            .unwrap();
        Assert::command(&["samtools", "view", t1])
            .stdout()
            .doesnt_contain("random_sequence_length_1000")
            .stdout()
            .contains("seq4\t77")
            .stdout()
            .contains("seq4\t141")
            .unwrap();
    }

    #[test]
    fn test_filter_unmapped_inverse_improper_pairs() {
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
                t2,
            ])
            .succeeds()
            .unwrap();
        Assert::command(&["samtools", "view", t1])
            .stdout()
            .contains("random_sequence_length_1000")
            .stdout()
            .contains("random_sequence_length_100\t77")
            .stdout()
            .contains("random_sequence_length_100\t141")
            .unwrap();
        Assert::command(&["samtools", "view", t2])
            .stdout()
            .contains("random_sequence_length_1000")
            .stdout()
            .contains("random_sequence_length_100\t77")
            .stdout()
            .contains("random_sequence_length_100\t141")
            .unwrap();
    }

    #[test]
    fn test_dot_in_extension() {
        // https://github.com/wwood/CoverM/issues/49
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--genome-fasta-directory",
                "tests/data/genomes_dir/",
                "-x",
                ".fna",
            ])
            .succeeds()
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
                td.path().to_str().unwrap(),
            ])
            .succeeds()
            .unwrap();
        assert!(td
            .path()
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
            .stdout()
            .is("Genome	shard1|shard2 Relative Abundance (%)
unmapped	0
genome3	25.024881
genome4	25.022575
genome5	0
genome6	25.020271
genome1	24.932274
genome2	0
")
            .succeeds()
            .unwrap()
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
            .stdout()
            .is("Contig	shard1|shard2 Mean
genome3~random_sequence_length_11001	0.110588886
genome4~random_sequence_length_11002	0.11057869
genome5~seq2	0
genome6~random_sequence_length_11003	0.11056851
genome1~random_sequence_length_11000	0.109861754
genome1~random_sequence_length_11010	0.110497236
genome2~seq1	0
")
            .succeeds()
            .unwrap()
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
            .stdout()
            .is("Contig	shard1|shard2 Mean
genome3~random_sequence_length_11001	0.110588886
genome4~random_sequence_length_11002	0.11057869
genome5~seq2	0
genome6~random_sequence_length_11003	0.11056851
genome1~random_sequence_length_11000	0.109861754
genome1~random_sequence_length_11010	0.110497236
genome2~seq1	0
")
            .succeeds()
            .unwrap()
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
                tf1.path().to_str().unwrap(),
            ])
            .stdout()
            .is("Genome	shard1|shard2 Relative Abundance (%)
unmapped	19.999998
genome3	0
genome4	26.699606
genome5	0
genome6	26.697144
genome1	26.60325
genome2	0
")
            .succeeds()
            .unwrap()
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
                tf1.path().to_str().unwrap(),
            ])
            .stdout()
            .is("Genome	shard1|shard2 Relative Abundance (%)
unmapped	19.999998
genome1	26.60325
genome2	0
genome3	0
genome4	26.699606
genome5	0
genome6	26.697144
")
            .succeeds()
            .unwrap()
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
                "50",
            ])
            .stderr()
            .contains(
                "coverm::contig] In sample '7seqs.fna/7seqs.reads_for_7_plus5random.1.fa', \
                 found 40 reads mapped out of 50 total (80.00%)",
            )
            .succeeds()
            .unwrap()
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
                "50",
            ])
            .stderr()
            .contains(
                "coverm::contig] In sample '7seqs.fna/7seqs.reads_for_7_plus5random.1.fa', \
                 found 40 reads mapped out of 50 total (80.00%)",
            )
            .succeeds()
            .unwrap()
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
            .stdout()
            .is(
                "Contig	shard1.fna|shard2.fna/7seqs.reads_for_7.1.fq|7seqs.reads_for_7.1.fq Mean\n\
                 genome3~random_sequence_length_11001	0.110588886\n\
                 genome4~random_sequence_length_11002	0.11057869\n\
                 genome5~seq2	0\n\
                 genome6~random_sequence_length_11003	0.11056851\n\
                 genome1~random_sequence_length_11000	0.109861754\n\
                 genome1~random_sequence_length_11010	0.110497236\n\
                 genome2~seq1	0\n",
            )
            .succeeds()
            .unwrap()
    }

    #[test]
    fn test_genome_definition_with_bam() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--genome-definition",
                "tests/data/7seqs.definition",
                "-b",
                "tests/data/7seqs.reads_for_seq1_and_seq2.bam",
            ])
            .succeeds()
            .stdout()
            .contains("Genome	7seqs.reads_for_seq1_and_seq2 Relative Abundance (%)\n")
            .stdout()
            .contains("genome2	53.167923\n")
            .stdout()
            .contains("genome5	46.832077\n")
            .unwrap();
    }

    #[test]
    fn test_genome_definition_with_reference() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--genome-definition",
                "tests/data/7seqs.definition",
                "-r",
                "tests/data/7seqs.fna",
                "--interleaved",
                "tests/data/reads_for_seq1_and_seq2.fna",
            ])
            .succeeds()
            .stdout()
            .contains("Genome	7seqs.fna/reads_for_seq1_and_seq2.fna Relative Abundance (%)\n")
            .stdout()
            .contains("genome2	53.167923\n")
            .stdout()
            .contains("genome5	46.832077\n")
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
                "tests/data/random.fq",
            ])
            .succeeds()
            .stdout()
            .is(
                "Genome	7seqs.fna/7seqs.reads_for_7.1.fq Relative Abundance (%)\t\
                          7seqs.fna/random.fq Relative Abundance (%)\n\
                          unmapped	0	100\n\
                          genome1	24.932272	NaN\n\
                          genome2	0	NaN\n\
                          genome3	25.02488	NaN\n\
                          genome4	25.022572	NaN\n\
                          genome5	0	NaN\n\
                          genome6	25.02027	NaN\n",
            )
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
                "tests/data/7seqs.fna.bwa1",
                "--output-directory",
                td.path().to_str().unwrap(),
            ])
            .succeeds()
            .unwrap();
        let bam = td
            .path()
            .join("7seqs.fna.bwa1.reads_for_seq1_and_seq2.1.fq.gz.bam");
        assert!(bam.is_file());
        Assert::command(&["samtools", "view", "-H", bam.to_str().unwrap()])
            .stdout()
            .contains("PN:bwa")
            .unwrap();
    }

    #[test]
    fn test_make_with_bwa_mem2() {
        let td = tempfile::TempDir::new().unwrap();
        Assert::main_binary()
            .with_args(&[
                "make",
                "--mapper",
                "bwa-mem2",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--reference",
                "tests/data/7seqs.fna",
                "--output-directory",
                td.path().to_str().unwrap(),
            ])
            .succeeds()
            .unwrap();
        let bam = td
            .path()
            .join("7seqs.fna.reads_for_seq1_and_seq2.1.fq.gz.bam");
        assert!(bam.is_file());
        Assert::command(&["samtools", "view", "-H", bam.to_str().unwrap()])
            .stdout()
            .contains("PN:bwa")
            .unwrap();
    }

    #[test]
    fn test_reference_specified_as_directory_genome() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--mapper",
                "bwa-mem",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--reference",
                "tests/data",
                "-x",
                "fna",
                "-s",
                "=",
            ])
            .fails()
            .stderr()
            .contains("should be a file, not e.g. a directory")
            .unwrap();
    }
    #[test]
    fn test_reference_not_existing_contig() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--mapper",
                "bwa-mem",
                "--coupled",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--reference",
                "testsblah",
            ])
            .fails()
            .stderr()
            .contains("does not appear to exist")
            .unwrap();
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
                td.path().to_str().unwrap(),
            ])
            .succeeds()
            .unwrap();
        let bam = td
            .path()
            .join("7seqs.fna.reads_for_seq1_and_seq2.1.fq.gz.bam");
        assert!(bam.is_file());
        Assert::command(&["samtools", "view", "-H", bam.to_str().unwrap()])
            .stdout()
            .contains("PN:minimap2")
            .unwrap();
    }

    #[test]
    fn test_contig_sparse_rpkm() {
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
                "--output-format",
                "sparse",
            ])
            .succeeds()
            .stdout()
            .is("Sample	Contig	RPKM	Reads per base	Length	Read Count\n\
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
            .stdout()
            .is("Genome	7seqs.fnaVbad_read RPKM\n\
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
                ctg4	1747.7994	0.024660854	5\n\
                ctg5	796.9696	0.0041760243	3\n\
                ctg6	140.38486	0.0021053297	1\n")
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
                ctg4	1747.7994	0.024660854	5	1747.7994	0.024660854	5\n\
                ctg5	796.9696	0.0041760243	3	796.9696	0.0041760243	3\n\
                ctg6	140.38486	0.0021053297	1	140.38486	0.0021053297	1\n")
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
    fn test_hifi_two_samples() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-m",
                "rpkm",
                "mean",
                "count",
                "-p",
                "minimap2-hifi",
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
                ctg4	1887.6234	0.020537598	3	1887.6234	0.020537598	3\n\
                ctg5	478.18173	0.0026476856	1	478.18173	0.0026476856	1\n\
                ctg6	252.69275	0.002099011	1	252.69275	0.002099011	1\n")
            .stderr()
            .contains(
                    "Running DB indexing command: \"minimap2\" \"-x\" \"map-hifi\" \"-t\" \"1\" \"-d\"")
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
                ctg4	2097.3594	0.067892104	6\n\
                ctg5	531.31305	0.010246328	2\n\
                ctg6	140.38486	0.0021596688	1\n")
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

    #[test]
    fn test_genome_all_methods() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--output-format",
                "sparse",
                "-b",
                "tests/data/7seqs.fnaVbad_read.bam",
                "--genome-fasta-directory",
                "tests/data/genomes_dir_7seqs/",
                "--genome-fasta-extension", "fasta",
                "-t", "5",
                "--methods", "covered_bases", "covered_fraction", "mean", "variance", "trimmed_mean", "rpkm", "relative_abundance", "length", 
                "--min-covered-fraction", "0"
            ])
            .succeeds()
            .stdout()
            .satisfies(|observed| assert_equal_table(
                "Sample	Genome	Covered Bases	Covered Fraction	Mean	Variance	Trimmed Mean	RPKM	Relative Abundance (%)	Length\n\
                7seqs.fnaVbad_read	unmapped	NA	NA	NA	NA	NA	NA	0	NA\n\
                7seqs.fnaVbad_read	genome2	899	0.899	1.6764706	0.51357985	1.6788511	500000	50	1000\n\
                7seqs.fnaVbad_read	genome6	0	0	0	0	0	0	0	11003\n\
                7seqs.fnaVbad_read	genome4	0	0	0	0	0	0	0	11002\n\
                7seqs.fnaVbad_read	genome3	0	0	0	0	0	0	0	11001\n\
                7seqs.fnaVbad_read	genome5	900	0.9	1.6764706	0.51357985	1.6788511	500000	50	1000\n\
                7seqs.fnaVbad_read	genome1	0	0	0	0	0	0	0	22010\n",
                observed
            ), "table incorrect")
            .unwrap();
    }

    #[test]
    fn test_dereplicate() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--genome-fasta-files",
                "tests/data/set1/1mbp.fna",
                "tests/data/set1/500kb.fna",
                "-t",
                "5",
                "--methods",
                "covered_fraction",
                "--min-covered-fraction",
                "0",
                "--dereplicate",
                "--single",
                "tests/data/set1/1read.actually_fasta.fq",
            ])
            .succeeds()
            .stdout()
            .is("Genome	1read.actually_fasta.fq Covered Fraction\n\
                1mbp	0.00232\n")
            .unwrap();
    }

    #[test]
    fn test_dereplicate_output_clusters() {
        let tf_clusters: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t_clusters = tf_clusters.path().to_str().unwrap();

        let tf_reps: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t_reps = tf_reps.path().to_str().unwrap();

        let td_symlink = tempfile::TempDir::new().unwrap();
        let td_copy = tempfile::TempDir::new().unwrap();

        Assert::main_binary()
            .with_args(&[
                "genome",
                "--genome-fasta-files",
                "tests/data/set1/1mbp.fna",
                "tests/data/set1/500kb.fna",
                "-t",
                "5",
                "--methods",
                "covered_fraction",
                "--min-covered-fraction",
                "0",
                "--dereplicate",
                "--single",
                "tests/data/set1/1read.actually_fasta.fq",
                "--dereplication-output-cluster-definition",
                t_clusters,
                "--dereplication-output-representative-fasta-directory",
                td_symlink.path().to_str().unwrap(),
                "--dereplication-output-representative-fasta-directory-copy",
                td_copy.path().to_str().unwrap(),
                "--dereplication-output-representative-list",
                t_reps,
            ])
            .succeeds()
            .stdout()
            .is("Genome	1read.actually_fasta.fq Covered Fraction\n\
                1mbp	0.00232\n")
            .unwrap();

        let mut s: String = "".to_string();
        std::fs::File::open(t_clusters)
            .unwrap()
            .read_to_string(&mut s)
            .unwrap();
        assert_eq!(
            "tests/data/set1/1mbp.fna	tests/data/set1/1mbp.fna\n\
                tests/data/set1/1mbp.fna	tests/data/set1/500kb.fna\n",
            s
        );

        let mut s2: String = "".to_string();
        std::fs::File::open(t_reps)
            .unwrap()
            .read_to_string(&mut s2)
            .unwrap();
        assert_eq!("tests/data/set1/1mbp.fna\n", s2);

        let out = td_symlink.path().join("1mbp.fna");
        assert!(out.exists());
        assert!(std::fs::symlink_metadata(out)
            .unwrap()
            .file_type()
            .is_symlink());
        assert!(!td_symlink.path().join("500kbp.fna").exists());

        let out2 = td_copy.path().join("1mbp.fna");
        assert!(out2.exists());
        assert!(!std::fs::symlink_metadata(out2)
            .unwrap()
            .file_type()
            .is_symlink());
        assert!(!td_symlink.path().join("500kbp.fna").exists());
    }

    #[test]
    fn test_dereplicate_checkm_ordering() {
        // 500kb specified first, should show up
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--genome-fasta-files",
                "tests/data/set1/500kb.fna",
                "tests/data/set1/1mbp.fna",
                "-t",
                "5",
                "--methods",
                "covered_fraction",
                "--min-covered-fraction",
                "0",
                "--dereplicate",
                "--single",
                "tests/data/set1/1read.actually_fasta.fq",
            ])
            .succeeds()
            .stdout()
            .is("Genome	1read.actually_fasta.fq Covered Fraction\n\
                500kb	0.00464\n")
            .unwrap();

        // 500kb specified first, but checkm specified so should be second
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--genome-fasta-files",
                "tests/data/set1/500kb.fna",
                "tests/data/set1/1mbp.fna",
                "-t",
                "5",
                "--methods",
                "covered_fraction",
                "--min-covered-fraction",
                "0",
                "--dereplicate",
                "--checkm-tab-table",
                "tests/data/set1/checkm.tsv",
                "--single",
                "tests/data/set1/1read.actually_fasta.fq",
            ])
            .succeeds()
            .stdout()
            .is("Genome	1read.actually_fasta.fq Covered Fraction\n\
                1mbp	0.00232\n")
            .unwrap();
    }

    #[test]
    fn test_dereplicate_genome_info_ordering() {
        // 500kb specified first, should show up
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--genome-fasta-files",
                "tests/data/set1/500kb.fna",
                "tests/data/set1/1mbp.fna",
                "-t",
                "5",
                "--methods",
                "covered_fraction",
                "--min-covered-fraction",
                "0",
                "--dereplicate",
                "--single",
                "tests/data/set1/1read.actually_fasta.fq",
            ])
            .succeeds()
            .stdout()
            .is("Genome	1read.actually_fasta.fq Covered Fraction\n\
                500kb	0.00464\n")
            .unwrap();

        // 500kb specified first, but checkm specified so should be second
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--genome-fasta-files",
                "tests/data/set1/500kb.fna",
                "tests/data/set1/1mbp.fna",
                "-t",
                "5",
                "--methods",
                "covered_fraction",
                "--min-covered-fraction",
                "0",
                "--dereplicate",
                "--genome-info",
                "tests/data/set1/genomeInfo.csv",
                "--single",
                "tests/data/set1/1read.actually_fasta.fq",
            ])
            .succeeds()
            .stdout()
            .is("Genome	1read.actually_fasta.fq Covered Fraction\n\
                1mbp	0.00232\n")
            .unwrap();
    }

    #[test]
    fn test_genome_fasta_list() {
        let mut tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();

        writeln!(tf, "tests/data/set1/500kb.fna").unwrap();
        writeln!(tf, "tests/data/set1/1mbp.fna").unwrap();
        tf.flush().unwrap();
        let t = tf.path().to_str().unwrap();

        Assert::main_binary()
            .with_args(&[
                "genome",
                "--genome-fasta-list",
                t,
                "-t",
                "5",
                "--methods",
                "covered_fraction",
                "--min-covered-fraction",
                "0",
                "--single",
                "tests/data/set1/1read.actually_fasta.fq",
            ])
            .succeeds()
            .stdout()
            .satisfies(
                |observed| {
                    assert_equal_table(
                        "Genome	1read.actually_fasta.fq Covered Fraction\n\
                500kb	0\n\
                1mbp	0.00232\n",
                        observed,
                    )
                },
                "table incorrect",
            )
            .unwrap();
    }

    #[test]
    fn test_contig_unsorted_bam_file() {
        Assert::main_binary()
            .with_args(&["contig", "-b", "tests/data/2seqs.bad_read.1.unsorted.bam"])
            .fails()
            .stderr()
            .contains("BAM file appears to be unsorted")
            .unwrap();
    }

    #[test]
    fn test_genome_separator_mode_unsorted_bam_file() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-s",
                "e",
                "-b",
                "tests/data/2seqs.bad_read.1.unsorted.bam",
            ])
            .fails()
            .stderr()
            .contains("BAM file appears to be unsorted")
            .unwrap();
    }

    #[test]
    fn test_genome_named_contigs_mode_unsorted_bam_file() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--genome-fasta-directory",
                "tests/data/genomes_dir",
                "-b",
                "tests/data/2seqs.bad_read.1.unsorted.bam",
            ])
            .fails()
            .stderr()
            .contains("BAM file appears to be unsorted")
            .unwrap();
    }

    #[test]
    fn test_no_zeros_bug1() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-c",
                "tests/data/rhys_bug/20120700_S3D.head100000.1.fq.gz",
                "tests/data/rhys_bug/20120700_S3D.head100000.2.fq.gz",
                "--genome-fasta-files",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.10.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.12.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.15.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.16.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.34.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.3.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.5.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.7.fna",
                "-t",
                "8",
                "-m",
                "mean",
                "covered_fraction",
                "--min-covered-fraction",
                "0.05",
                "--exclude-supplementary",
            ])
            .succeeds()
            .stdout()
            .satisfies(
                |observed| {
                    assert_equal_table(
                        "Genome	20120700_S3D.head100000.1.fq.gz Mean	20120700_S3D.head100000.1.fq.gz Covered Fraction\n\
                        73.20120700_S3D.10\t0.071023874\t0.06777273\n73.20120700_S3D.12\t0\t0\n73.20120700_S3D.15\t0\t0\n73.20120700_S3D.16\t0\t0\n73.20120700_S3D.3\t0\t0\n73.20120700_S3D.34\t0.06653676\t0.0630154\n73.20120700_S3D.5\t0.1341526\t0.123165175\n73.20120700_S3D.7\t0.100108385\t0.093486056\n\
                        ",
                        observed,
                    )
                },
                "table incorrect",
            )
            .unwrap();
    }

    #[test]
    fn test_no_zeros_bug2() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-c",
                "tests/data/rhys_bug/20120700_S3D.head100000.1.fq.gz",
                "tests/data/rhys_bug/20120700_S3D.head100000.2.fq.gz",
                "--genome-fasta-files",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.10.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.12.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.15.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.16.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.34.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.3.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.5.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.7.fna",
                "-t",
                "8",
                "--no-zeros",
                "-m",
                "mean",
                "covered_fraction",
                "--min-covered-fraction",
                "0.05",
                "--exclude-supplementary",
            ])
            .succeeds()
            .stdout()
            .satisfies(
                |observed| {
                    assert_equal_table(
                        "Genome	20120700_S3D.head100000.1.fq.gz Mean	20120700_S3D.head100000.1.fq.gz Covered Fraction\n\
                        73.20120700_S3D.10\t0.071023874\t0.06777273\n73.20120700_S3D.34\t0.06653676\t0.0630154\n73.20120700_S3D.5\t0.1341526\t0.123165175\n73.20120700_S3D.7\t0.100108385\t0.093486056\n",
                        observed,
                    )
                },
                "table incorrect",
            )
            .unwrap();
    }

    #[test]
    fn test_no_zeros_bug3() {
        Assert::main_binary()
        .with_args(&[
            "genome",
            "-c",
            "tests/data/rhys_bug/20120700_S3D.head100000.1.fq.gz",
            "tests/data/rhys_bug/20120700_S3D.head100000.2.fq.gz",
            "--genome-fasta-files",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.10.fna",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.12.fna",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.15.fna",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.16.fna",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.34.fna",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.3.fna",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.5.fna",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.7.fna",
            "-t",
            "8",
            "--no-zeros",
            "-m",
            "mean",
            "covered_fraction",
            "--min-covered-fraction",
            "0.03",
            "--exclude-supplementary",
        ])
        .succeeds()
        .stdout()
        .satisfies(
            |observed| {
                assert_equal_table(
                    "Genome	20120700_S3D.head100000.1.fq.gz Mean	20120700_S3D.head100000.1.fq.gz Covered Fraction\n\
                    73.20120700_S3D.10\t0.071023874\t0.06777273\n73.20120700_S3D.15\t0.03561887\t0.034370355\n73.20120700_S3D.16\t0.032864396\t0.031665392\n73.20120700_S3D.3\t0.036180563\t0.03499215\n73.20120700_S3D.34\t0.06653676\t0.0630154\n73.20120700_S3D.5\t0.1341526\t0.123165175\n73.20120700_S3D.7\t0.100108385\t0.093486056\n\
                    ",
                    observed,
                )
            },
            "table incorrect",
        )
        .unwrap();
    }

    #[test]
    fn test_no_zeros_bug4() {
        Assert::main_binary()
        .with_args(&[
            "genome",
            "-c",
            "tests/data/rhys_bug/20120700_S3D.head100000.1.fq.gz",
            "tests/data/rhys_bug/20120700_S3D.head100000.2.fq.gz",
            "--genome-fasta-files",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.10.fna",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.12.fna",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.15.fna",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.16.fna",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.34.fna",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.3.fna",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.5.fna",
            "tests/data/rhys_bug/genomes/73.20120700_S3D.7.fna",
            "-t",
            "8",
            //"--no-zeros",
            "-m",
            "mean",
            "covered_fraction",
            "--min-covered-fraction",
            "0.03",
            "--exclude-supplementary",
        ])
        .succeeds()
        .stdout()
        .satisfies(
            |observed| {
                assert_equal_table(
                    "Genome	20120700_S3D.head100000.1.fq.gz Mean	20120700_S3D.head100000.1.fq.gz Covered Fraction\n\
                    73.20120700_S3D.10\t0.071023874\t0.06777273\n73.20120700_S3D.12\t0\t0\n73.20120700_S3D.15\t0.03561887\t0.034370355\n73.20120700_S3D.16\t0.032864396\t0.031665392\n73.20120700_S3D.3\t0.036180563\t0.03499215\n73.20120700_S3D.34\t0.06653676\t0.0630154\n73.20120700_S3D.5\t0.1341526\t0.123165175\n73.20120700_S3D.7\t0.100108385\t0.093486056\n\
                    ",
                    observed,
                )
            },
            "table incorrect",
        )
        .unwrap();
    }

    #[test]
    fn test_contig_output_file() {
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();

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
                "-o",
                t,
            ])
            .succeeds()
            .stdout()
            .is("")
            .unwrap();

        let mut buf = vec![];
        std::fs::File::open(tf.path())
            .unwrap()
            .read_to_end(&mut buf)
            .unwrap();

        assert_eq!(
            "Sample\tContig\tMean\n\
            2seqs.fasta/bad_reads.interleaved.fq\tseq1\t0.899\n\
            2seqs.fasta/bad_reads.interleaved.fq\tseq2\t0\n",
            str::from_utf8(&buf).unwrap()
        )
    }

    #[test]
    fn test_autoconcatenation_with_clashing() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-m",
                "mean",
                "--genome-fasta-files",
                "tests/data/contig_name_clashing/genome1.fna",
                "tests/data/contig_name_clashing/genome2.fna",
                "tests/data/contig_name_clashing/genome3.fna",
                "--coupled",
                "tests/data/contig_name_clashing/reads_for_genome2.1.fa",
                "tests/data/contig_name_clashing/reads_for_genome2.2.fa",
            ])
            .succeeds()
            .stdout()
            .is("Genome	reads_for_genome2.1.fa Mean\n\
            genome1	0\n\
            genome2	4.8142858\n\
            genome3	0\n\
            ")
            .unwrap();
    }

    #[test]
    fn test_clashing_contig_names_with_reference_specified() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-m",
                "mean",
                "--genome-fasta-files",
                "tests/data/contig_name_clashing/genome1.fna",
                "tests/data/contig_name_clashing/genome2.fna",
                "tests/data/contig_name_clashing/genome3.fna",
                "--reference",
                "tests/data/contig_name_clashing/bad_reference.fna",
                "--coupled",
                "tests/data/contig_name_clashing/reads_for_genome2.1.fa",
                "tests/data/contig_name_clashing/reads_for_genome2.2.fa",
            ])
            .fails()
            .stderr()
            .contains(
                "The contig 'random_sequence_length_500_1' has been assigned to multiple genomes",
            )
            .unwrap();
    }

    #[test]
    fn test_tpm_contig_sparse() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--output-format",
                "sparse",
                "-m",
                "mean",
                "tpm",
                "-b",
                "tests/data/tpm_test.bam",
            ])
            .succeeds()
            .stdout()
            .is("Sample	Contig	Mean	TPM\n\
                tpm_test	genome1~random_sequence_length_11000	0	0\n\
                tpm_test	genome1~random_sequence_length_11010	0	0\n\
                tpm_test	genome2~seq1	1.5882353	900000.0357627869\n\
                tpm_test	genome3~random_sequence_length_11001	0	0\n\
                tpm_test	genome4~random_sequence_length_11002	0	0\n\
                tpm_test	genome5~seq2	0.14467005	99999.99403953552\n\
                tpm_test	genome6~random_sequence_length_11003	0	0\n")
            .unwrap();
    }

    #[test]
    fn test_tpm_contig_dense() {
        Assert::main_binary()
            .with_args(&[
                "contig",
                "-m",
                "mean",
                "tpm",
                "-b",
                "tests/data/tpm_test.bam",
            ])
            .succeeds()
            .stdout()
            .is("Contig	tpm_test Mean	tpm_test TPM\n\
                genome1~random_sequence_length_11000	0	0\n\
                genome1~random_sequence_length_11010	0	0\n\
                genome2~seq1	1.5882353	900000.06\n\
                genome3~random_sequence_length_11001	0	0\n\
                genome4~random_sequence_length_11002	0	0\n\
                genome5~seq2	0.14467005	99999.99\n\
                genome6~random_sequence_length_11003	0	0\n")
            .unwrap();
    }

    #[test]
    fn test_tpm_genome_sparse() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "--output-format",
                "sparse",
                "-m",
                "mean",
                "tpm",
                "-b",
                "tests/data/tpm_test.bam",
                "-s",
                "~",
                "--min-covered-fraction",
                "0",
            ])
            .succeeds()
            .stdout()
            .is("Sample	Genome	Mean	TPM\n\
                tpm_test	genome1	0	0\n\
                tpm_test	genome2	1.5882353	900000.0357627869\n\
                tpm_test	genome3	0	0\n\
                tpm_test	genome4	0	0\n\
                tpm_test	genome5	0.14467005	99999.99403953552\n\
                tpm_test	genome6	0	0\n")
            .unwrap();
    }

    #[test]
    fn test_tpm_genome_dense() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-m",
                "mean",
                "tpm",
                "-b",
                "tests/data/tpm_test.bam",
                "-s",
                "~",
                "--min-covered-fraction",
                "0",
            ])
            .succeeds()
            .stdout()
            .is("Genome	tpm_test Mean	tpm_test TPM\n\
            genome1	0	0\n\
            genome2	1.5882353	900000.06\n\
            genome3	0	0\n\
            genome4	0	0\n\
            genome5	0.14467005	99999.99\n\
            genome6	0	0\n")
            .unwrap();
    }

    #[test]
    fn test_mismatched_read_pairs() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-m",
                "mean",
                "tpm",
                "-c",
                "tests/data/bad_read.1.fa",
                "tests/data/7seqs.fna",
                "-r",
                "tests/data/7seqs.fna",
                "--single-genome",
            ])
            .fails()
            .stderr()
            .contains("Not continuing since when input file pairs have unequal numbers of reads")
            .unwrap();
    }

    #[test]
    fn test_single_genome_no_exclude_supplementary() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-m",
                "count",
                "-b",
                "tests/data/2seqs.bad_read.1.with_supplementary.bam",
                "--single-genome",
                "--min-covered-fraction",
                "0",
            ])
            .succeeds()
            .stdout()
            // We only count a supplementary alignment once in read count
            .is("Genome	2seqs.bad_read.1.with_supplementary Read Count\n\
            genome1	20\n")
            .unwrap();
    }

    #[test]
    fn test_genomes_and_contigs_with_supplementary() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-m",
                "mean",
                "covered_fraction",
                "count",
                "--genome-fasta-files",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.10.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.12.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.15.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.16.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.34.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.3.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.5.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.7.fna",
                "-c",
                "tests/data/rhys_bug/20120700_S3D.stray_read1.1.fq",
                "tests/data/rhys_bug/20120700_S3D.stray_read1.2.fq",
                "--min-covered-fraction",
                "0",
            ])
            .succeeds()
            .stdout()
            .satisfies(
                |observed| {
                    assert_equal_table(
                        "Genome	20120700_S3D.stray_read1.1.fq Mean	20120700_S3D.stray_read1.1.fq Covered Fraction	20120700_S3D.stray_read1.1.fq Read Count\n\
                        73.20120700_S3D.10	0.000022701126	0.00003879494	2\n\
                        73.20120700_S3D.12	0	0	0\n\
                        73.20120700_S3D.15	0	0	0\n\
                        73.20120700_S3D.16	0	0	0\n\
                        73.20120700_S3D.34	0	0	0\n\
                        73.20120700_S3D.3	0	0	0\n\
                        73.20120700_S3D.5	0.000043860742	0.000043714655	2\n\
                        73.20120700_S3D.7	0	0	0\n\
                        ",
                        observed,
                    )
                },
                "table incorrect",
            )
            .stderr()
            .contains("found 4 reads mapped out of 4 total (100.00%)")
            .unwrap();
    }

    #[test]
    fn test_genomes_and_contigs_without_supplementary() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-m",
                "mean",
                "covered_fraction",
                "count",
                "--genome-fasta-files",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.10.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.12.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.15.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.16.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.34.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.3.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.5.fna",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.7.fna",
                "-c",
                "tests/data/rhys_bug/20120700_S3D.stray_read1.1.fq",
                "tests/data/rhys_bug/20120700_S3D.stray_read1.2.fq",
                "--min-covered-fraction",
                "0",
                "--exclude-supplementary"
            ])
            .succeeds()
            .stdout()
            .satisfies(
                |observed| {
                    assert_equal_table(
                        "Genome	20120700_S3D.stray_read1.1.fq Mean	20120700_S3D.stray_read1.1.fq Covered Fraction	20120700_S3D.stray_read1.1.fq Read Count\n\
                        73.20120700_S3D.10	0.000008399416	0.000024585164	2\n\
                        73.20120700_S3D.12	0	0	0\n\
                        73.20120700_S3D.15	0	0	0\n\
                        73.20120700_S3D.16	0	0	0\n\
                        73.20120700_S3D.34	0	0	0\n\
                        73.20120700_S3D.3	0	0	0\n\
                        73.20120700_S3D.5	0.000043860742	0.000043714655	2\n\
                        73.20120700_S3D.7	0	0	0\n\
                        ",
                        observed,
                    )
                },
                "table incorrect",
            )
            .stderr()
            .contains("found 4 reads mapped out of 4 total (100.00%)")
            .unwrap();
    }

    #[test]
    fn test_genomes_separator_with_supplementary() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-m",
                "mean",
                "covered_fraction",
                "count",
                "-r",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.10.fna",
                "-s",
                "_",
                "-c",
                "tests/data/rhys_bug/20120700_S3D.stray_read1.1.fq",
                "tests/data/rhys_bug/20120700_S3D.stray_read1.2.fq",
                "--min-covered-fraction",
                "0",
            ])
            .succeeds()
            .stdout()
            .satisfies(
                |observed| {
                    assert_equal_table(
                        "Genome	73.20120700_S3D.10.fna/20120700_S3D.stray_read1.1.fq Mean	73.20120700_S3D.10.fna/20120700_S3D.stray_read1.1.fq Covered Fraction	73.20120700_S3D.10.fna/20120700_S3D.stray_read1.1.fq Read Count\n\
                        73.20120700	0.000022701126	0.00003879494	2\n\
                        ",
                        observed,
                    )
                },
                "table incorrect",
            )
            .stderr()
            .contains("found 2 reads mapped out of 4 total (50.00%)")
            .unwrap();
    }

    #[test]
    fn test_genomes_separator_without_supplementary() {
        Assert::main_binary()
            .with_args(&[
                "genome",
                "-m",
                "mean",
                "covered_fraction",
                "count",
                "-r",
                "tests/data/rhys_bug/genomes/73.20120700_S3D.10.fna",
                "-s",
                "_",
                "-c",
                "tests/data/rhys_bug/20120700_S3D.stray_read1.1.fq",
                "tests/data/rhys_bug/20120700_S3D.stray_read1.2.fq",
                "--min-covered-fraction",
                "0",
                "--exclude-supplementary",
            ])
            .succeeds()
            .stdout()
            .satisfies(
                |observed| {
                    assert_equal_table(
                        "Genome	73.20120700_S3D.10.fna/20120700_S3D.stray_read1.1.fq Mean	73.20120700_S3D.10.fna/20120700_S3D.stray_read1.1.fq Covered Fraction	73.20120700_S3D.10.fna/20120700_S3D.stray_read1.1.fq Read Count\n\
                        73.20120700	0.000008399416	0.000024585164	2\n\
                        ",
                        observed,
                    )
                },
                "table incorrect",
            )
            .stderr()
            .contains("found 2 reads mapped out of 4 total (50.00%)")
            .unwrap();
    }

    #[test]
    fn test_completion_generation() {
        Assert::main_binary()
            .with_args(&["shell-completion", "-o", "/dev/stdout", "--shell", "bash"])
            .succeeds()
            .stdout()
            .contains("genome")
            .unwrap()
    }
}

// TODO: Add mismatching bases test
// TODO: Filter fails when reference sequences are duplicated?
// TODO: Filter should spit things out if no thresholds are specified.
