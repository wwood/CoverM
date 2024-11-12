use crate::bam_generator::{complete_processes, name_stoit};
use crate::mapping_parameters::MappingParameters;
use crate::{coverage_printer, coverage_takers::*, genomes_and_contigs, OutputWriter};
// use crate::mosdepth_genome_coverage_estimators::*;

use csv;

struct MappingResultStruct {
    tmpfile: tempfile::NamedTempFile,
    stoit_name: String,
}

// The function to actually do the mapping with strobealign --aemb, returning
// the coverage results as a vector of MappingResultStructs.
fn strobealign_aemb_coverage_result(
    mapping_parameters: MappingParameters,
    threads: u16,
) -> Vec<MappingResultStruct> {
    let mut commands = vec![];
    let mut log_files = vec![];

    let mut mapping_results: Vec<MappingResultStruct> = vec![];
    for p1 in mapping_parameters {
        for p2 in p1 {
            debug!("Processing mapping parameters: {:?}", p2);
            let mapping_result = tempfile::Builder::new()
                .prefix("coverm-strobealign-aemb-mapping-result")
                .tempfile()
                .unwrap_or_else(|_| panic!("Failed to create strobealign-aemb result tempfile"));
            let mapping_log = tempfile::Builder::new()
                .prefix("coverm-strobealign-aemb-mapping-log")
                .tempfile()
                .unwrap_or_else(|_| panic!("Failed to create strobealign-aemb log tempfile"));
            let cmd = format!(
                "strobealign --aemb -t {} '{}' '{}' '{}' > '{}' 2> '{}'",
                threads,
                p2.reference,
                p2.read1,
                p2.read2.unwrap_or(""),
                mapping_result
                    .path()
                    .to_str()
                    .expect("Failed to get strobealign-aemb result path"),
                mapping_log
                    .path()
                    .to_str()
                    .expect("Failed to get strobealign-aemb log path"),
            );
            log_files.push(mapping_log);
            commands.push(cmd);
            mapping_results.push(MappingResultStruct {
                tmpfile: mapping_result,
                stoit_name: name_stoit(p2.reference, p2.read1, true),
            });
        }
    }

    // Run each command in serial, and complete each process
    // before starting the next one.
    for cmd_string in commands {
        debug!("Queuing cmd_string: {}", cmd_string);
        let mut cmd = std::process::Command::new("bash");
        cmd.arg("-c")
            .arg(&cmd_string)
            .stderr(std::process::Stdio::piped());
        let process = cmd.spawn().expect("Unable to execute bash");
        complete_processes(
            vec![process],
            vec![cmd_string],
            vec!["strobealign-aemb log file".to_string()],
            vec![log_files.remove(0)],
            None,
        );
    }
    mapping_results
}

pub fn strobealign_aemb_coverage_per_contig(
    mapping_parameters: MappingParameters,
    coverage_taker: &mut CoverageTakerType,
    // _coverage_estimators: &mut [CoverageEstimator],
    // print_zero_coverage_contigs: bool,
    threads: u16,
    coverage_printer: &mut coverage_printer::CoveragePrinter,
    print_stream: &mut OutputWriter,
) {
    let mapping_results = strobealign_aemb_coverage_result(mapping_parameters, threads);

    // Read the mapping results and pass them to the coverage taker
    for mapping_result_tempfile in mapping_results {
        debug!("Starting stoid {}", mapping_result_tempfile.stoit_name);
        coverage_taker.start_stoit(&mapping_result_tempfile.stoit_name);
        // Read the TSV format
        let file = std::fs::File::open(mapping_result_tempfile.tmpfile.path())
            .expect("Failed to open strobealign-aemb mapping result");
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(file);
        for (i, result) in rdr.records().enumerate() {
            let record = result.expect("Failed to read strobealign-aemb mapping result");
            if record.len() != 2 {
                panic!(
                    "Unexpected number of columns in strobealign-aemb mapping result line {}: {:?}",
                    i, record
                );
            }
            coverage_taker.start_entry(i, &record[0]);
            coverage_taker.add_single_coverage(
                record[1]
                    .parse()
                    .expect("Failed to parse strobealign-aemb coverage"),
            );
            coverage_taker.finish_entry();
        }
    }
    coverage_printer.finalise_printing(coverage_taker, print_stream, None, &vec![], None, None);
}

pub fn strobealign_aemb_coverage_per_genome(
    mapping_parameters: MappingParameters,
    genomes_and_contigs: &Option<genomes_and_contigs::GenomesAndContigs>,
    coverage_taker: &mut CoverageTakerType,
    threads: u16,
    coverage_printer: &mut coverage_printer::CoveragePrinter,
    print_stream: &mut OutputWriter,
) {
    let mapping_results = strobealign_aemb_coverage_result(mapping_parameters, threads);

    // Read the mapping results and pass them to the coverage taker
    for mapping_result_tempfile in mapping_results {
        debug!("Starting stoid {}", mapping_result_tempfile.stoit_name);
        coverage_taker.start_stoit(&mapping_result_tempfile.stoit_name);
        // Read the TSV format
        let file = std::fs::File::open(mapping_result_tempfile.tmpfile.path())
            .expect("Failed to open strobealign-aemb mapping result");
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(file);

        let num_genomes = match &genomes_and_contigs {
            Some(gc) => gc.genomes.len(),
            None => 1, // single genome
        };
        let mut per_genome_coverage_and_length: Vec<(f32, u64)> = vec![(0., 0); num_genomes];
        for (i, result) in rdr.records().enumerate() {
            let record = result.expect("Failed to read strobealign-aemb mapping result");
            if record.len() != 2 {
                panic!(
                    "Unexpected number of columns in strobealign-aemb mapping result line {}: {:?}",
                    i, record
                );
            }
            let genome_index = match genomes_and_contigs {
                Some(gc) => gc
                    .genome_index_of_contig(&record[0])
                    .expect("Failed to get genome index for strobealign-aemb mapping result"),
                None => 0,
            };
            debug!("Genome index: {}", genome_index);
            let coverage: f32 = record[1]
                .parse()
                .expect("Failed to parse strobealign-aemb coverage");
            let contig_length = 1; // FIXME: Get contig length from somewhere
            per_genome_coverage_and_length[genome_index] = (
                per_genome_coverage_and_length[genome_index].0 + coverage,
                per_genome_coverage_and_length[genome_index].1 + contig_length,
            );
        }

        for (genome_index, (coverage, length)) in per_genome_coverage_and_length.iter().enumerate()
        {
            coverage_taker.start_entry(
                genome_index,
                match genomes_and_contigs {
                    Some(gc) => &gc.genomes[genome_index],
                    None => "genome1", // as per genome.rs
                },
            );
            let coverage = match length {
                0 => 0.,
                _ => coverage / (*length as f32),
            };
            coverage_taker.add_single_coverage(coverage);
            coverage_taker.finish_entry();
        }
    }
    coverage_printer.finalise_printing(coverage_taker, print_stream, None, &vec![], None, None);
}
