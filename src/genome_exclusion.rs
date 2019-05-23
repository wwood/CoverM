use std::collections::HashSet;

//use genomes_and_contigs::GenomesAndContigs;
use find_first;

pub enum GenomeExcluders<'a> {
    SeparatorGenomeExclusionFilter {
        split_char: u8,
        excluded_genomes: HashSet<&'a [u8]>
    },
    NoExclusionGenomeFilter{}
}

pub trait GenomeExclusion {
    fn is_excluded(&self, contig_name: &[u8]) -> bool;
}

// pub struct GenomesAndContigsExclusionFilter<'a> {
//     genomes_and_contigs: &'a GenomesAndContigs
// }

// impl<'a> GenomeExclusion for GenomesAndContigsExclusionFilter<'a> {
//     fn is_excluded(&self, contig_name: &[u8]) -> bool {
//         self.genomes_and_contigs.genome_index_of_contig(contig_name).is_some()
//     }
// }

pub struct SeparatorGenomeExclusionFilter<'a> {
    pub split_char: u8,
    pub excluded_genomes: HashSet<&'a [u8]>,
}

impl<'a> GenomeExclusion for SeparatorGenomeExclusionFilter<'a> {
    fn is_excluded(&self, contig_name: &[u8]) -> bool {
        debug!("contig name {:?}, separator {:?}", contig_name, self.split_char);
        let offset = find_first(contig_name, self.split_char).expect(
            &format!("Contig name {} does not contain split symbol, so cannot determine which genome it belongs to",
                     self.split_char));
        let genome = &contig_name[(0..offset)];
        return self.excluded_genomes.contains(genome);
    }
}

pub struct NoExclusionGenomeFilter {}
impl GenomeExclusion for NoExclusionGenomeFilter {
    fn is_excluded(&self, _contig_name: &[u8]) -> bool { false }
}

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn test_genomes_and_contigs_exclusion_filter() {
    //     let mut contig_to_genome = GenomesAndContigs::new();
    //     let genome = String::from("genome0");
    //     let index = contig_to_genome.establish_genome(genome);
    //     contig_to_genome.insert(String::from("contig1"), index);
    //     contig_to_genome.insert(String::from("contig2"), index);

    //     let ex = GenomesAndContigsExclusionFilter {
    //         genomes_and_contigs: &contig_to_genome
    //     };

    //     assert_eq!(ex.is_excluded(&"contig1".to_string()), true);
    //     assert_eq!(ex.is_excluded(&"contig2".to_string()), true);
    //     assert_eq!(ex.is_excluded(&"contig20".to_string()), false);
    // }

    #[test]
    fn test_separator_exclusion_filter() {
        let mut hashset: HashSet<&[u8]> = HashSet::new();
        hashset.insert(b"genomeYes");
        let ex = SeparatorGenomeExclusionFilter {
            split_char: b"="[0],
            excluded_genomes: hashset
        };
        assert_eq!(true, ex.is_excluded(b"genomeYes=contig1"));
        assert_eq!(false, ex.is_excluded(b"genomeNo=contig1"));
    }
}
