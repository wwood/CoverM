use std::collections::HashMap;

use find_first;

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
        match self.contig_to_genome.get(&contig_name) {
            Some(previous_index) => {
                let genome_prev = &self.genomes[*previous_index];
                let genome_current = &self.genomes[genome_index];
                panic!("The contig '{}' has been assigned to multiple genomes, \
                        at least '{}' and '{}'. You may try not using \
                        --reference and let coverm generate a reference of \
                        concatenated contigs, or rename the contigs in your \
                        genome file(s).",
                       contig_name,
                       genome_prev,
                       genome_current);
            },
            None => {}
        }
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
