use crate::bam_generator::complete_processes;
use crate::coverage_takers::*;
use crate::mapping_parameters::MappingParameters;
use crate::mosdepth_genome_coverage_estimators::*;
use crate::ReadsMapped;

pub fn strobealign_aemb_coverage<T: CoverageTaker>(
    mapping_parameters: MappingParameters,
    _coverage_taker: &mut T,
    _coverage_estimators: &mut [CoverageEstimator],
    // print_zero_coverage_contigs: bool,
    threads: u16,
) -> Vec<ReadsMapped> {
    let mut commands = vec![];
    let mut log_files = vec![];

    for p1 in mapping_parameters {
        for p2 in p1 {
            let mapping_result = tempfile::Builder::new()
                .prefix("coverm-strobealign-aemb-mapping-result")
                .tempfile()
                .unwrap_or_else(|_| panic!("Failed to create strobealign-aemb result tempfile"));
            let mapping_log = tempfile::Builder::new()
                .prefix("coverm-strobealign-aemb-mapping-log")
                .tempfile()
                .unwrap_or_else(|_| panic!("Failed to create strobealign-aemb log tempfile"));
            let cmd = format!(
                "strobealign --aemb -t {} {} {} {} > {} 2> {}",
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

    // Read the mapping results and pass them to the coverage taker
    panic!();
}
