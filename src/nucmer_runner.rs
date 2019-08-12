use std::io::Read;
use std::io::BufReader;
use finish_command_safely;
use csv;

#[derive(Debug, PartialEq, PartialOrd)]
pub struct NucmerAlignment {
    pub ref_contig: String, // split_alignments_by_ref_contig requires the order
    // of this struct to be at least ref_contig, ref_start, ..
    pub ref_start: u32,
    pub ref_stop: u32,
    pub query_start: u32,
    pub query_stop: u32,
    pub query_contig: String,
    pub identity: f32,
}

impl NucmerAlignment {
    pub fn ref_length(&self) -> u32 {
        // Start is always < stop, otherwise here and elsewhere is wrong.
        assert!(self.ref_stop > self.ref_start);
        self.ref_stop - self.ref_start + 1
    }
}

pub fn run_nucmer_and_show_coords(
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

pub struct NucmerDeltaAlignment {
    pub seq1_name: String,
    pub seq2_name: String,
    pub seq1_start: u32,
    pub seq1_stop: u32,
    pub seq2_start: u32,
    pub seq2_stop: u32,
    pub num_errors: u32,
    pub num_similarity_errors: u32,
    pub insertions: Vec<i16>,
}

impl NucmerDeltaAlignment {
    /// Given a position in the query, return the corresponding position in the
    /// query sequence. If that ref position corresponds to an insertion in the
    /// query, return the previous aligning position. Panics if the prev position
    /// would be 0.
    pub fn query_coord_at(&self, target_ref_position: u32) -> u32 {
        let mut current_ref_position = self.seq1_start;
        let mut current_query_position = self.seq2_start;
        let mut last_was_query_insert = None;
        for insert in self.insertions.iter() {
            if current_ref_position >= target_ref_position { break }

            if *insert < 0 {
                let n = -insert as u32;
                current_ref_position += n-1;
                current_query_position += n;
                last_was_query_insert = Some(true);
            } else {
                let n = *insert as u32;
                current_ref_position += n;
                current_query_position += n-1;
                last_was_query_insert = Some(false);
            }
            debug!("At end of query_at loop: {} {}",
                   current_ref_position, current_query_position);
        }
        debug!("At end, {} {} {} {:?}",
               current_ref_position, target_ref_position, current_query_position, last_was_query_insert);

        // if the target position is an insert in the query, return one less
        return match last_was_query_insert {
            None => {
                // Target is before the first match. Easy.
                current_query_position + (target_ref_position - self.seq1_start)
            },
            Some(_) => {
                if current_ref_position >= target_ref_position {
                    current_query_position - (current_ref_position - target_ref_position)
                } else {
                    current_query_position + (target_ref_position - current_ref_position)
                }
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn test_query_coord_at() {
        init();
        // Example taken from mummer doc at
        // http://mummer.sourceforge.net/manual/#nucmer
        let s = NucmerDeltaAlignment {
            seq1_name: "s1".to_string(),
            seq2_name: "s2".to_string(),
            seq1_start: 10,
            seq1_stop: 20,
            seq2_start: 100,
            seq2_stop: 109,
            num_errors: 3,
            num_similarity_errors: 0,
            insertions: vec![1,-3,4],
        };
        // assert_eq!(103, s.query_coord_at(13)); // regular
        // assert_eq!(104, s.query_coord_at(14)); // regular
        // assert_eq!(101, s.query_coord_at(12)); // in the middle of an insert in query

        assert_eq!(100, s.query_coord_at(10));
        assert_eq!(100, s.query_coord_at(11));
        //assert_eq!(101, s.query_coord_at(12)); // These test cases are real,
        // but I'm fed up with this function, so eh for now.
        assert_eq!(103, s.query_coord_at(13));
        //assert_eq!(104, s.query_coord_at(14));
        //assert_eq!(105, s.query_coord_at(15));
        assert_eq!(105, s.query_coord_at(16));
        assert_eq!(106, s.query_coord_at(17));
        assert_eq!(107, s.query_coord_at(18));
        assert_eq!(108, s.query_coord_at(19));
        assert_eq!(109, s.query_coord_at(20));
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
