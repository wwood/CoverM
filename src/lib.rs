pub mod contig;
pub mod genome;
pub mod mosdepth_genome_coverage_estimators;
pub mod genomes_and_contigs;
pub mod bam_generator;
pub mod filter;

extern crate bio;
#[macro_use]
extern crate log;

extern crate rust_htslib;
extern crate env_logger;
extern crate nix;
extern crate tempdir;

use std::path::Path;
use genomes_and_contigs::GenomesAndContigs;



pub fn read_genome_fasta_files(fasta_file_paths: &[&str]) -> GenomesAndContigs {
    let mut contig_to_genome = GenomesAndContigs::new();

    for file in fasta_file_paths {
        let path = Path::new(*file);
        let reader = bio::io::fasta::Reader::from_file(path)
            .expect(&format!("Unable to read fasta file {}", file));

        let genome_name = String::from(
            path.file_stem()
                .expect("Problem while determining file stem")
            .to_str()
                .expect("File name string conversion problem"));
        if contig_to_genome.genome_index(&genome_name).is_some() {
            panic!("The genome name {} was derived from >1 file", genome_name);
        }
        let genome_index = contig_to_genome.establish_genome(genome_name);
        for record in reader.records() {
            let contig = String::from(
                record
                    .expect(&format!("Failed to parse contig name in fasta file {:?}", path))
                    .id());
            contig_to_genome.insert(contig, genome_index);
        }
    }
    return contig_to_genome;
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contig_to_genome(){
        let mut contig_to_genome = GenomesAndContigs::new();
        let genome = String::from("genome0");
        let index = contig_to_genome.establish_genome(genome);
        contig_to_genome.insert(String::from("contig1"), index);
        assert_eq!(
            String::from("genome0"),
            *(contig_to_genome.genome_of_contig(&String::from("contig1")).unwrap()));
    }

    #[test]
    fn test_read_genome_fasta_files_one_genome(){
        let contig_to_genome = read_genome_fasta_files(&vec!["tests/data/genome1.fna"]);
        assert_eq!(String::from("genome1"), *contig_to_genome.genome_of_contig(&String::from("seq1")).unwrap());
        assert_eq!(String::from("genome1"), *contig_to_genome.genome_of_contig(&String::from("seq2")).unwrap());
    }
}
