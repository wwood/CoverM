extern crate assert_cli;

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
                "/tmp/o",
                ]).succeeds().unwrap();
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
                "-r",
                t_full,
                "-1",
                "tests/data/reads_for_seq1_and_seq2.1.fq.gz",
                "-2",
                "tests/data/reads_for_seq1_and_seq2.2.fq.gz",
                "--remove-headers",
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
    fn test_contig_multiple_references(){
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--remove-headers",
                "--min-percent-identity",
                "0.95",
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
    fn test_coupled_reads_input(){
        Assert::main_binary()
            .with_args(&[
                "contig",
                "--remove-headers",
                "--min-percent-identity",
                "0.95",
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
}
