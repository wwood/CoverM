// NOTE: This file is shared as a direct copy of code between cockatoo and coverm

use std::collections::HashMap;
use std::process;

#[derive(Serialize, Deserialize, Debug)]
pub struct GenomesAndContigs {
    pub genomes: Vec<String>,
    pub contig_to_genome: HashMap<String, usize>,
}

impl GenomesAndContigs {
    pub fn new() -> GenomesAndContigs {
        GenomesAndContigs {
            genomes: vec![],
            contig_to_genome: HashMap::new(),
        }
    }

    pub fn establish_genome(&mut self, genome_name: String) -> usize {
        let index = self.genomes.len();
        self.genomes.push(genome_name);
        index
    }

    pub fn insert(&mut self, contig_name: String, genome_index: usize) {
        if let Some(previous_index) = self.contig_to_genome.get(&contig_name) {
            let genome_prev = &self.genomes[*previous_index];
            let genome_current = &self.genomes[genome_index];
            error!(
                "The contig '{}' has been assigned to multiple genomes, \
                    at least '{}' and '{}'. You may try not using \
                    --reference/--bam-files and let coverm generate a reference of \
                    concatenated contigs, or rename the contigs in your \
                    genome file(s).",
                contig_name, genome_prev, genome_current
            );
            process::exit(1);
        }
        self.contig_to_genome.insert(contig_name, genome_index);
    }

    pub fn genome_index(&self, genome_name: &String) -> Option<usize> {
        self.genomes
            .iter()
            .position(|genome| genome.eq(genome_name))
    }

    pub fn genome_of_contig(&self, contig_name: &String) -> Option<&String> {
        match self.contig_to_genome.get(contig_name) {
            Some(index) => Some(&self.genomes[*index]),
            None => None,
        }
    }

    pub fn genome_index_of_contig(&self, contig_name: &String) -> Option<usize> {
        self.contig_to_genome.get(contig_name).copied()
    }
}

impl Default for GenomesAndContigs {
    fn default() -> Self {
        Self::new()
    }
}
