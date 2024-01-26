use std::collections::HashSet;
use std::str;

use genomes_and_contigs::GenomesAndContigs;

pub enum GenomeExcluders<'a> {
    SeparatorGenomeExclusionFilter {
        split_char: u8,
        excluded_genomes: HashSet<&'a [u8]>,
    },
    NoExclusionGenomeFilter {},
}

pub trait GenomeExclusion {
    fn is_excluded(&self, contig_name: &[u8]) -> bool;
}

pub struct GenomesAndContigsExclusionFilter<'a> {
    pub genomes_and_contigs: &'a GenomesAndContigs,
    pub excluded_genomes: HashSet<&'a [u8]>,
}

impl<'a> GenomeExclusion for GenomesAndContigsExclusionFilter<'a> {
    fn is_excluded(&self, contig_name: &[u8]) -> bool {
        let contig_str = str::from_utf8(contig_name).unwrap().to_string();
        match self.genomes_and_contigs.genome_of_contig(&contig_str) {
            Some(g) => {
                if self.excluded_genomes.contains(&g.as_bytes()) {
                    debug!(
                        "Excluding contig '{}' as it is part of excluded genome '{}'",
                        str::from_utf8(contig_name).unwrap(),
                        g
                    );
                    true
                } else {
                    false
                }
            }
            None => false,
        }
    }
}

pub struct SeparatorGenomeExclusionFilter<'a> {
    pub split_char: u8,
    pub excluded_genomes: HashSet<&'a [u8]>,
}

impl<'a> GenomeExclusion for SeparatorGenomeExclusionFilter<'a> {
    fn is_excluded(&self, contig_name: &[u8]) -> bool {
        debug!(
            "contig name {:?}, separator {:?}",
            contig_name, self.split_char
        );
        let offset = contig_name.iter().position(|&c|c == self.split_char).unwrap_or_else(|| panic!("Contig name {} does not contain split symbol, so cannot determine which genome it belongs to",
                     self.split_char));
        let genome = &contig_name[0..offset];
        self.excluded_genomes.contains(genome)
    }
}

pub struct NoExclusionGenomeFilter {}
impl GenomeExclusion for NoExclusionGenomeFilter {
    fn is_excluded(&self, _contig_name: &[u8]) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genomes_and_contigs_exclusion_filter() {
        let mut contig_to_genome = GenomesAndContigs::new();
        let genome = String::from("genome0");
        let index = contig_to_genome.establish_genome(genome);
        contig_to_genome.insert(String::from("contig1"), index);
        contig_to_genome.insert(String::from("contig2"), index);

        let mut hashset: HashSet<&[u8]> = HashSet::new();
        hashset.insert(b"genome0");

        let ex = GenomesAndContigsExclusionFilter {
            genomes_and_contigs: &contig_to_genome,
            excluded_genomes: hashset,
        };

        assert!(ex.is_excluded(b"contig1"));
        assert!(ex.is_excluded(b"contig2"));
        assert!(!ex.is_excluded(b"contig20"));
    }

    #[test]
    fn test_separator_exclusion_filter() {
        let mut hashset: HashSet<&[u8]> = HashSet::new();
        hashset.insert(b"genomeYes");
        let ex = SeparatorGenomeExclusionFilter {
            split_char: b"="[0],
            excluded_genomes: hashset,
        };
        assert!(ex.is_excluded(b"genomeYes=contig1"));
        assert!(!ex.is_excluded(b"genomeNo=contig1"));
    }
}
