// NOTE: This file is shared as a direct copy of code between cockatoo and coverm

use genomes_and_contigs::GenomesAndContigs;
use needletail::parse_fastx_file;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;
use std::process;

pub fn read_genome_fasta_files(
    fasta_file_paths: &Vec<&str>,
    use_full_sequence_name: bool,
) -> GenomesAndContigs {
    let mut contig_to_genome = GenomesAndContigs::new();

    // NOTE: A lot of this code is shared with mapping_index_maintenance.rs#generate_concatenated_fasta_file
    for file in fasta_file_paths {
        let path = Path::new(file);
        let mut reader =
            parse_fastx_file(path).unwrap_or_else(|_| panic!("Unable to read fasta file {}", file));

        // Remove .gz .bz .xz from file names if present
        let mut genome_name1 =
            String::from(path.to_str().expect("File name string conversion problem"));
        if let Some(i) = genome_name1.rfind(".gz") {
            genome_name1.truncate(i);
        } else if let Some(i) = genome_name1.rfind(".bz") {
            genome_name1.truncate(i);
        } else if let Some(i) = genome_name1.rfind(".xz") {
            genome_name1.truncate(i);
        }
        let path1 = Path::new(&genome_name1);

        let genome_name = String::from(
            path1
                .file_stem()
                .expect("Problem while determining file stem")
                .to_str()
                .expect("File name string conversion problem"),
        );
        if contig_to_genome.genome_index(&genome_name).is_some() {
            error!("The genome name {genome_name} was derived from >1 file");
            process::exit(1);
        }
        let genome_index = contig_to_genome.establish_genome(genome_name);
        while let Some(record) = reader.next() {
            let record_expected = record
                .unwrap_or_else(|_| panic!("Failed to parse record in fasta file {:?}", path));

            if record_expected.format() != needletail::parser::Format::Fasta {
                panic!(
                    "File {:?} is not a fasta file, but a {:?}",
                    path,
                    record_expected.format()
                );
            }

            let contig_name = String::from(
                std::str::from_utf8(record_expected.id())
                    .expect("UTF-8 conversion problem in contig name"),
            );
            let contig = match use_full_sequence_name {
                false => match contig_name.split_once(' ') {
                    Some((contig, _)) => contig.to_string(),
                    None => contig_name,
                },
                true => contig_name,
            };
            contig_to_genome.insert(contig, genome_index);
        }
    }
    contig_to_genome
}

pub fn read_genome_definition_file(definition_file_path: &str) -> GenomesAndContigs {
    let f = std::fs::File::open(definition_file_path).unwrap_or_else(|_| {
        panic!(
            "Unable to find/read genome definition file {}",
            definition_file_path
        )
    });
    let file = std::io::BufReader::new(&f);
    let mut contig_to_genome: HashMap<String, String> = HashMap::new();
    let mut genome_to_contig: HashMap<String, Vec<String>> = HashMap::new();
    // Maintain the same order as the input file.
    let mut genome_order: Vec<String> = vec![];

    for line_res in file.lines() {
        let line = line_res.expect("Read error on genome definition file");
        let v: Vec<&str> = line.split('\t').collect();
        if v.len() == 2 {
            let genome = v[0].trim();
            let contig = v[1]
                .split_ascii_whitespace()
                .next()
                .expect("Failed to split contig name by whitespace in genome definition file");
            if contig_to_genome.contains_key(contig) {
                if contig_to_genome[contig] != genome {
                    error!(
                        "The contig name '{contig}' was assigned to multiple genomes"
                    );
                    process::exit(1);
                }
            } else {
                contig_to_genome.insert(contig.to_string(), genome.to_string());
            }

            if genome_to_contig.contains_key(genome) {
                genome_to_contig
                    .get_mut(genome)
                    .unwrap()
                    .push(contig.to_string());
            } else {
                genome_to_contig.insert(genome.to_string(), vec![contig.to_string()]);
                genome_order.push(genome.to_string());
            }
        } else if v.is_empty() {
            continue;
        } else {
            error!(
                "The line \"{line}\" in the genome definition file is not a \
                    genome name and contig name separated by a tab"
            );
            process::exit(1);
        }
    }

    info!(
        "Found {} contigs assigned to {} different genomes from \
           the genome definition file",
        contig_to_genome.len(),
        genome_to_contig.len()
    );

    let mut gc = GenomesAndContigs::new();
    for genome in genome_order {
        let contigs = &genome_to_contig[&genome];
        let genome_index = gc.establish_genome(genome);
        for contig in contigs {
            gc.insert(contig.to_string(), genome_index);
        }
    }
    gc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contig_to_genome() {
        let mut contig_to_genome = GenomesAndContigs::new();
        let genome = String::from("genome0");
        let index = contig_to_genome.establish_genome(genome);
        contig_to_genome.insert(String::from("contig1"), index);
        assert_eq!(
            String::from("genome0"),
            *(contig_to_genome
                .genome_of_contig(&String::from("contig1"))
                .unwrap())
        );
    }

    #[test]
    fn test_read_genome_fasta_files_one_genome() {
        let contig_to_genome = read_genome_fasta_files(&vec!["tests/data/genome1.fna"], false);
        assert_eq!(
            String::from("genome1"),
            *contig_to_genome
                .genome_of_contig(&String::from("seq1"))
                .unwrap()
        );
        assert_eq!(
            String::from("genome1"),
            *contig_to_genome
                .genome_of_contig(&String::from("seq2"))
                .unwrap()
        );
    }

    #[test]
    fn test_read_genome_definition_file() {
        let contig_to_genome = read_genome_definition_file("tests/data/7seqs.definition");
        assert_eq!(
            Some(&String::from("genome4")),
            contig_to_genome
                .genome_of_contig(&String::from("genome4~random_sequence_length_11002"))
        );
        assert_eq!(6, contig_to_genome.genomes.len());
    }

    #[test]
    fn test_read_genome_definition_file_with_comments() {
        let contig_to_genome =
            read_genome_definition_file("tests/data/7seqs.definition_with_comments");
        assert_eq!(
            Some(&String::from("genome4")),
            contig_to_genome
                .genome_of_contig(&String::from("genome4~random_sequence_length_11002"))
        );
        assert_eq!(6, contig_to_genome.genomes.len());
    }
}
