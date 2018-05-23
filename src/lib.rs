pub mod contig;
pub mod genome;
pub mod mosdepth_genome_coverage_estimators;

extern crate bio;
#[macro_use]
extern crate log;

extern crate rust_htslib;
extern crate env_logger;

use std::collections::HashMap;
// use std::path::Path;
// use std::boxed::Box;
// //use std::ffi::OsStr;
// use std::collections::BTreeSet;
use std::option::Option::*;

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

    pub fn insert(&mut self, contig_name: String, genome_name: String) {
        let genome_index_option = self.genome_index(&genome_name);
        match genome_index_option {
            Some(index) => {
                self.contig_to_genome.insert(contig_name, index);
            },
            None => {
                let index = self.genomes.len();
                self.genomes.push(genome_name);
                self.contig_to_genome.insert(contig_name, index);
            }
        };
    }

    fn genome_index(&self, genome_name: &String) -> Option<usize> {
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

// pub fn read_genome_fasta_files<'a, 'b, 'c>(fasta_file_paths: &'a [&'b str]) -> GenomesAndContigs {
//     let mut contig_to_genome = GenomesAndContigs::new();

//     for file in fasta_file_paths {
//         //let genome_path_box = String::from(*file);
//         let path = Path::new(*file);
//         let reader = bio::io::fasta::Reader::from_file(path)
//             .expect(&format!("Unable to read fasta file {}", file));
//         // TODO: Detect if the same genome name comes from multiple files.
//         let genome_name = Box::new(
//             String::from(path.file_stem().expect("Problem while determining file stem")
//                 .to_str().expect("File name string conversion problem")
//             ));
//         for record in reader.records() {
//             let contig = String::from(record
//                                        .expect(&format!("Failed to parse contig name in fasta file {:?}", path))
//                                        .id());
//             match contig_to_genome.contains_key(&contig) {
//                 true => panic!(format!(
//                     "The contig name {} is contained in multiple genomes. Please rename it.", contig)),
//                 false => {
//                     contig_to_genome.insert(contig, Box::new(String::from(&genome_name)));
//                 }
//             }
//         }
//     }
//     return contig_to_genome;
// }


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_genome_fasta_files_one_genome(){
        let mut contig_to_genome = GenomesAndContigs::new();
        contig_to_genome.insert(String::from("contig1"), String::from("genome0"));
        assert_eq!(
            String::from("genome0"),
            *(contig_to_genome.genome_of_contig(&String::from("contig1")).unwrap()));
    }

    // #[test]
    // fn test_read_genome_fasta_files_one_genome(){
    //     let contig_to_genome = read_genome_fasta_files(&vec!["test/data/genome1.fna"]);
    //     assert_eq!(Box::new(String::from("genome1")), contig_to_genome[&String::from("seq1")]);
    //     assert_eq!(Box::new(String::from("genome1")), contig_to_genome[&String::from("seq2")]);
    // }
}
