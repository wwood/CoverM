#[macro_use]
extern crate assert_cli;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;

    #[test]
    fn test_filter_all_reads(){
        Assert::main_binary()
            .with_args(&[
                "filter",
                "-b",
                "test/data/2seqs.bad_read.1.bam",
                "-o",
                "/tmp/o"]).succeeds().unwrap();
        assert_cmd!(samtools view "/tmp/o").stdout().contains("1\t99\tseq1")
            .unwrap();
    }

    #[test]
    fn test_filter_filter_out(){
        Assert::main_binary()
            .with_args(&[
                "filter",
                "--min-percent-identity",
                "0.99",
                "-b",
                "test/data/2seqs.bad_read.1.bam",
                "-o",
                "/tmp/o"]).succeeds().unwrap();
        assert_cmd!(samtools view "/tmp/o").stdout().doesnt_contain("1\t99\tseq1")
            .unwrap();
    }
}
