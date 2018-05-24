pub mod contig;
pub mod genome;
pub mod mosdepth_genome_coverage_estimators;

extern crate bio;
#[macro_use]
extern crate log;

extern crate rust_htslib;
extern crate env_logger;

use std::collections::HashMap;
use std::path::Path;


#[derive(Debug)]
pub struct GenomesAndContigs {
    genomes: Vec<String>,
    contig_to_genome: HashMap<String, usize>
}


impl GenomesAndContigs {
    pub fn new() -> GenomesAndContigs {
        GenomesAndContigs {
            genomes: vec!(),
            contig_to_genome: HashMap::new()
        }
    }

    pub fn establish_genome(&mut self, genome_name: String) -> usize {
        let index = self.genomes.len();
        self.genomes.push(genome_name);
        return index
    }

    pub fn insert(&mut self, contig_name: String, genome_index: usize) {
        self.contig_to_genome.insert(contig_name, genome_index);
    }

    pub fn genome_index(&self, genome_name: &String) -> Option<usize> {
        match find_first(self.genomes.as_slice(), genome_name.to_string()) {
            Ok(index) => Some(index),
            Err(_) => None
        }
    }

    pub fn genome_of_contig(&self, contig_name: &String) -> Option<&String> {
        match self.contig_to_genome.get(contig_name) {
            Some(index) => Some(&self.genomes[*index]),
            None => None
        }
    }
}
/// Finds the first occurence of element in a slice
fn find_first<T>(slice: &[T], element: T) -> Result<usize, &'static str>
    where T: std::cmp::PartialEq<T> {

    let mut index: usize = 0;
    for el in slice {
        if *el == element {
            return Ok(index)
            //let res: Result<usize, None> = Ok(index)
        }
        index += 1;
    }
    return Err("Element not found in slice")
}

pub fn read_genome_fasta_files<'a, 'b, 'c>(fasta_file_paths: &'a [&'b str]) -> GenomesAndContigs {
    let mut contig_to_genome = GenomesAndContigs::new();

    for file in fasta_file_paths {
        //let genome_path_box = String::from(*file);
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
        let contig_to_genome = read_genome_fasta_files(&vec!["test/data/genome1.fna"]);
        assert_eq!(String::from("genome1"), *contig_to_genome.genome_of_contig(&String::from("seq1")).unwrap());
        assert_eq!(String::from("genome1"), *contig_to_genome.genome_of_contig(&String::from("seq2")).unwrap());
    }
}
