use std::collections::HashSet;
use std::string::String;

use genomes_and_contigs::GenomesAndContigs;
use find_first_iter;

pub trait GenomeExclusion {
    fn is_excluded(&self, contig_name: &String) -> bool;
}

pub struct GenomesAndContigsExclusionFilter<'a> {
    genomes_and_contigs: &'a GenomesAndContigs
}

impl<'a> GenomesAndContigsExclusionFilter<'a> {
    pub fn generate_from_genomes_and_contigs(genomes_and_contigs: &'a GenomesAndContigs)
                                             -> GenomesAndContigsExclusionFilter<'a> {
        GenomesAndContigsExclusionFilter {
            genomes_and_contigs: genomes_and_contigs
        }
    }
}

impl<'a> GenomeExclusion for GenomesAndContigsExclusionFilter<'a> {
    fn is_excluded(&self, contig_name: &String) -> bool {
        self.genomes_and_contigs.genome_index_of_contig(contig_name).is_some()
    }
}

pub struct SeparatorGenomeExclusionFilter {
    split_char: char,
    excluded_genomes: HashSet<String>
}

impl GenomeExclusion for SeparatorGenomeExclusionFilter {
    fn is_excluded(&self, contig_name: &String) -> bool {
        debug!("contig name {:?}, separator {:?}", contig_name, self.split_char);
        let offset = find_first_iter(contig_name.chars(), self.split_char).expect(
            &format!("Contig name {} does not contain split symbol, so cannot determine which genome it belongs to",
                     self.split_char));
        let genome = &contig_name[(0..offset)];
        return self.excluded_genomes.contains(genome);
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

        let ex = GenomesAndContigsExclusionFilter::generate_from_genomes_and_contigs(&contig_to_genome);

        assert_eq!(ex.is_excluded(&"contig1".to_string()), true);
        assert_eq!(ex.is_excluded(&"contig2".to_string()), true);
        assert_eq!(ex.is_excluded(&"contig20".to_string()), false);
    }

    #[test]
    fn test_separator_exclusion_filter() {
        let mut hashset: HashSet<String> = HashSet::new();
        hashset.insert("genomeYes".to_string());
        let ex = SeparatorGenomeExclusionFilter {
            split_char: '=',
            excluded_genomes: hashset
        };
        assert_eq!(true, ex.is_excluded(&"genomeYes=contig1".to_string()));
        assert_eq!(false, ex.is_excluded(&"genomeNo=contig1".to_string()));
    }
}
