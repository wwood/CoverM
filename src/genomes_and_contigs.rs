use std;
use std::collections::HashMap;

#[derive(Debug)]
pub struct GenomesAndContigs {
    pub genomes: Vec<String>,
    pub contig_to_genome: HashMap<String, usize>
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

    pub fn genome_index_of_contig(&self, contig_name: &String) -> Option<usize> {
        match self.contig_to_genome.get(contig_name) {
            Some(index) => Some(*index),
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
