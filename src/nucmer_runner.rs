use std::io::prelude::*;
use std::io::BufReader;
use finish_command_safely;
use csv;

#[derive(Debug, PartialEq)]
pub struct NucmerAlignment {
    pub ref_start: u32,
    pub ref_stop: u32,
    pub query_start: u32,
    pub query_stop: u32,
    pub identity: f32,
    pub ref_contig: String,
    pub query_contig: String,
}

fn run_nucmer_and_show_coords(
    ref_suffix: &str,
    ref_file: &str,
    query_file: &str)
-> Vec<NucmerAlignment>{

    let mut cmd = std::process::Command::new("bash");
    let real_cmd = format!(
        "set -e -o pipefail; nucmer --delta=/dev/stdout --load={} {} {} |show-coords -T /dev/stdin |tail -n+5",
        ref_suffix, ref_file, query_file);
    cmd.arg("-c").arg(real_cmd)
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());

    debug!("Running nucmer cmd: {:?}", &cmd);
    // Parsnp can fail at steps after MUMi part, and that's OK.
    let mut process = cmd.spawn().expect(&format!("Failed to spawn {}", "nucmer"));
    process.wait()
        .expect(&format!("Failed to glean exitstatus from failing {} process", "nucmer"));
    // if the stdout has certain words, we are in business.
    let nucmer_stdout = process.stdout.as_mut().unwrap();
    let nucmer_stdout_reader = BufReader::new(nucmer_stdout);
    let aligns = parse_show_coords(nucmer_stdout_reader);
    finish_command_safely(process, "nucmer");
    return aligns;
}


fn parse_show_coords<R: Read>(
    reader: R) -> Vec<NucmerAlignment> {

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(reader);

    let mut to_return = Vec::new();
    for record_res in rdr.records() {
        match record_res {
            Ok(record) => {
                let aln = NucmerAlignment {
                    ref_start: record[0].parse::<u32>().expect("Unexpected format of show-coords output"),
                    ref_stop: record[1].parse::<u32>().expect("Unexpected format of show-coords output"),
                    query_start: record[2].parse::<u32>().expect("Unexpected format of show-coords output"),
                    query_stop: record[3].parse::<u32>().expect("Unexpected format of show-coords output"),
                    identity: record[6].parse::<f32>().expect("Unexpected format of show-coords output"),
                    ref_contig: record[7].to_string(),
                    query_contig: record[8].to_string(),
                };
                //debug!("Found nucmer alignment: {:?}", &aln);
                to_return.push(aln);
            },
            Err(e) => {
                error!("Error parsing nucmer/show-coords output: {}", e);
                std::process::exit(1);
            }
        }
    }
    return to_return;
}


#[cfg(test)]
mod tests {
    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn test_nucmer_parser_hello_world() {
        init();
        let parsed = run_nucmer_and_show_coords(
            "tests/data/parsnp/1_first_group/73.20120800_S1D.21.fna",
            "tests/data/parsnp/1_first_group/73.20120800_S1D.21.fna",
            "tests/data/parsnp/1_first_group/73.20120600_E3D.30.fna",
        );
        assert_eq!(674-4, parsed.len());
        assert!(parsed.iter().find(|a| **a == NucmerAlignment {
            ref_start: 908,
            ref_stop: 1599,
            query_start: 21120,
            query_stop: 20429,
            identity: 88.58,
            ref_contig: "73.20120800_S1D.21_contig_20912".to_string(),
            query_contig: "73.20120600_E3D.30_contig_20710".to_string(),
        }).is_some());
    }
}
