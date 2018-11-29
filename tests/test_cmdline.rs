extern crate assert_cli;
extern crate tempfile;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;
    extern crate tempfile;
    use std;

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
                "--min-percent-identity",
                "0.99",
                "-b",
                "tests/data/2seqs.bad_read.1.bam",
                "-o",
                "/tmp/o"]).succeeds().unwrap();
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
                "--min-percent-identity",
                "0.95",
                "--contig-end-exclusion",
                "0",
                "-r",
                t_full,
                "-1",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "-2",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
            ])
            .succeeds()
            .stdout().contains(format!(
                "{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11000\t0.0
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11010\t0.0
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome2~seq1\t1.2
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome3~random_sequence_length_11001\t0.0
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome4~random_sequence_length_11002\t0.0
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome5~seq2\t1.2
{}/reads_for_seq1_and_seq2.1.fq.gz\tgenome6~random_sequence_length_11003\t0.0",
                t, t, t, t, t, t, t).as_str()).unwrap();
    }

    #[test]
    fn test_contig_multiple_references(){
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--min-percent-identity",
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
            .stdout().contains("7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11000\t0.0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11010\t0.0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome2~seq1\t1.2
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome3~random_sequence_length_11001\t0.0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome4~random_sequence_length_11002\t0.0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome5~seq2\t1.2
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome6~random_sequence_length_11003\t0.0
2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz\tseq1\t1.2
2seqs.fasta/reads_for_seq1_and_seq2.1.fq.gz\tseq2\t1.2").unwrap();
    }

    #[test]
    fn test_coupled_reads_input(){
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--min-percent-identity",
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
            .stdout().contains("7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11000\t0.0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11010\t0.0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome2~seq1\t1.2
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome3~random_sequence_length_11001\t0.0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome4~random_sequence_length_11002\t0.0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome5~seq2\t1.2
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome6~random_sequence_length_11003\t0.0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11000\t0.0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome1~random_sequence_length_11010\t0.0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome2~seq1\t1.2
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome3~random_sequence_length_11001\t0.0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome4~random_sequence_length_11002\t0.0
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome5~seq2\t1.2
7seqs.fna/reads_for_seq1_and_seq2.1.fq.gz\tgenome6~random_sequence_length_11003\t0.0
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
                "--interleaved",
                "tests/data/bad_reads.interleaved.fq",
                "--no-flag-filter",
            ])
            .succeeds()
            .stdout().contains("2seqs.fasta/bad_reads.interleaved.fq\tseq1\t0.899
2seqs.fasta/bad_reads.interleaved.fq\tseq2\t0.0").unwrap();
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
                "--min-aligned-length",
                "300",
                "--contig-end-exclusion",
                "0",
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
                "-r",
                "tests/data/2seqs.fasta",
                "--single",
                "tests/data/bad_reads.interleaved.fq",
                "--no-flag-filter",
                "--contig-end-exclusion",
                "0",
            ])
            .succeeds()
            .stdout().contains("2seqs.fasta/bad_reads.interleaved.fq\tseq1\t0.899
2seqs.fasta/bad_reads.interleaved.fq\tseq2\t0.0").unwrap();
    }

    #[test]
    fn test_genome_coupled_read_input_argparsing(){
        // Had trouble (because of a bug in clap? with this previously)
        Assert::main_binary()
            .with_args(&[
                "genome",
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
                "--reference",
                "tests/data/7seqs.fna",
                "--bam-file-cache-directory",
                td.path().to_str().unwrap()
            ]).succeeds().unwrap();
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
}
