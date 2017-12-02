extern crate bio;
#[macro_use]
extern crate log;

extern crate rust_htslib;
use self::rust_htslib::bam;
use self::rust_htslib::bam::Read as bamRead;
use self::rust_htslib::bam::HeaderView;

use std::collections::HashMap;
use std::str;

pub fn genome_coverage(bam_files: &Vec<&str>, fasta_files: &Vec<&str>, print_stream: &mut std::io::Write){
    // Read in fasta files into a map of contig -> genome
    let contig_to_genome = read_contig_to_genome_from_fasta_files(fasta_files);
    println!("contig_to_genome: {:?}", contig_to_genome);

    let mut process_contig = |bam_file, total_coverage, bases_seen, tid, header: &HeaderView, genome_to_total_coverage: &mut HashMap<String, u32>, contig_to_genome: &HashMap<String, String>| {
        debug!("For {} in bam {}, found {} total coverage", tid, bam_file, total_coverage);
        let total_bases = header.target_len(tid).unwrap() as f32;
        let target_names = header.target_names();
        let contig_name = target_names.get(tid as usize).unwrap();
        let contig_name_str = str::from_utf8(contig_name).unwrap();
        let genome: String = contig_to_genome.get(contig_name_str).unwrap().clone();
        let cov = genome_to_total_coverage.entry(genome).or_insert(0);
        *cov += total_coverage;
        writeln!(print_stream, "{}\t{}", bam_file, total_coverage as f32/total_bases).unwrap();
    };

    let mut genome_to_total_coverage: HashMap<String, u32> = HashMap::new();

    for bam_file in bam_files {
        debug!("Working on BAM file {}", bam_file);
        let bam = bam::Reader::from_path(bam_file).unwrap();
        let header = bam.header();

        let mut total_coverage = 0;
        let mut bases_covered = 0;
        let mut last_tid: u32 = 0;
        for p in bam.pileup() {
            let pileup = p.unwrap();

            let tid = pileup.tid();
            if tid != last_tid {
                process_contig(bam_file, total_coverage, bases_covered, tid, header, &mut genome_to_total_coverage, &contig_to_genome);
                last_tid = tid;
                total_coverage = 0;
                bases_covered = 0;
            }
            bases_covered += 1; //TODO: What happens with indels?
            total_coverage += pileup.depth();
        }
        process_contig(bam_file, total_coverage, bases_covered, last_tid, header, &mut genome_to_total_coverage, &contig_to_genome);
        println!("{:?}", genome_to_total_coverage);
    }

}

fn read_contig_to_genome_from_fasta_files(fasta_files: &Vec<&str>) -> HashMap<String, String> {
    let mut contig_to_genome: HashMap<String, String> = HashMap::new();

    for fasta_file in fasta_files {
        let fasta_reader = bio::io::fasta::Reader::from_file(fasta_file).expect("Error opening FASTA file"); //TODO: Say which fasta file failed
        let genome = fasta_file;
        for record in fasta_reader.records() {
            let contig: String = String::from(record.unwrap().id());
            if contig_to_genome.contains_key(&contig) {
                let found_genome = contig_to_genome.get(&contig).unwrap();
                if genome != found_genome {
                    panic!(
                        "The contig '{}' is found in at least two genomes: '{}' and '{}",
                        contig, found_genome, genome)
                }
            } else {
                let genome_string: String = genome.to_string();
                contig_to_genome.insert(contig, genome_string);
            };
        }
    }

    return contig_to_genome;
}
