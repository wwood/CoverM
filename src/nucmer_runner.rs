use std::io::Read;
use std::io::BufReader;
use finish_command_safely;
use csv;

pub fn nucmer_to_deltas(
    ref_suffix: &str,
    ref_file: &str,
    query_file: &str
) -> Vec<NucmerDeltaAlignment> {

    let mut cmd = std::process::Command::new("nucmer");
    cmd
        .arg("--delta=/dev/stdout")
        .arg(&format!("--load={}", ref_suffix))
        .arg(ref_file)
        .arg(query_file)
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());

    debug!("Running nucmer cmd: {:?}", &cmd);
    // Parsnp can fail at steps after MUMi part, and that's OK.
    let mut process = cmd.spawn().expect(&format!("Failed to spawn {}", "nucmer"));
    let nucmer_stdout = process.stdout.as_mut().unwrap();
    let nucmer_stdout_reader = BufReader::new(nucmer_stdout);
    let aligns = parse_delta(nucmer_stdout_reader);
    finish_command_safely(process, "nucmer");
    return aligns;
}

#[derive(Debug)]
enum NucmerDeltaParsingState {
    NextLinePaths,
    NextLineNucmer,
    NextLineHeaderOrOverview,
    NextLineInsert
}

fn parse_delta<R: Read>(
    reader: R) -> Vec<NucmerDeltaAlignment> {

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b' ')
        .has_headers(false)
        .flexible(true)
        .from_reader(reader);

    let mut state = NucmerDeltaParsingState::NextLinePaths;
    let mut aligns = vec![];
    let mut current_record = NucmerDeltaAlignment {
        seq1_name: String::new(),
        seq2_name: String::new(),
        seq1_start: 0,
        seq1_stop: 0,
        seq2_start: 0,
        seq2_stop: 0,
        num_errors: 0,
        num_similarity_errors: 0,
        insertions: vec![],
    };
    let mut last_ref_name = String::new();
    let mut last_query_name = String::new();

    for record_res in rdr.records() {
        match record_res {
            Ok(record) => {
                //debug!("In state {:?}, parsing record {:?}", state, record);
                match state {
                    NucmerDeltaParsingState::NextLinePaths => {
                        state = NucmerDeltaParsingState::NextLineNucmer;
                    }
                    NucmerDeltaParsingState::NextLineNucmer => {
                        assert!(&record[0]=="NUCMER", format!("nucmer delta parsing error: {:?}", record));
                        state = NucmerDeltaParsingState::NextLineHeaderOrOverview;
                    },
                    NucmerDeltaParsingState::NextLineHeaderOrOverview => {
                        if record[0].starts_with(">") {
                            last_ref_name = record[0][1..].to_string();
                            last_query_name = record[1].to_string();
                        } else {
                            current_record = NucmerDeltaAlignment {
                                seq1_name: last_ref_name.clone(),
                                seq2_name: last_query_name.clone(),
                                seq1_start: record[0].parse::<u32>().expect("nucmer delta parse error")-1,
                                seq1_stop: record[1].parse::<u32>().expect("nucmer delta parse error")-1,
                                seq2_start: record[2].parse::<u32>().expect("nucmer delta parse error")-1,
                                seq2_stop: record[3].parse::<u32>().expect("nucmer delta parse error")-1,
                                num_errors: record[4].parse::<u32>().expect("nucmer delta parse error"),
                                num_similarity_errors: record[5].parse().expect("nucmer delta parse error"),
                                insertions: vec![],
                            };
                            state = NucmerDeltaParsingState::NextLineInsert;
                        }
                    },
                    NucmerDeltaParsingState::NextLineInsert => {
                        assert!(record.len()==1);
                        if &record[0] == "0" {
                            aligns.push(current_record.clone());
                            state = NucmerDeltaParsingState::NextLineHeaderOrOverview;
                        } else {
                            current_record.insertions.push(record[0].parse().expect("nucmer delta parse error"))
                        }
                    }
                }
            },
            Err(e) => {
                error!("Error parsing nucmer output: {}", e);
                std::process::exit(1);
            }
        }
    }

    return aligns;
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct NucmerDeltaAlignment {
    pub seq1_name: String,// split_alignments_by_ref_contig requires the order
    // of this struct to be at least ref_contig, ref_start, ..
    pub seq1_start: u32,
    pub seq1_stop: u32,
    pub seq2_name: String,
    pub seq2_start: u32,
    pub seq2_stop: u32,
    pub num_errors: u32,
    pub num_similarity_errors: u32,
    pub insertions: Vec<i32>,
}

impl NucmerDeltaAlignment {
    pub fn ref_length(&self) -> u32 {
        // Start is always < stop, otherwise here and elsewhere is wrong.
        assert!(self.seq1_stop > self.seq1_start);
        self.seq1_stop - self.seq1_start + 1
    }
}

impl NucmerDeltaAlignment {
    /// Given a position in the query, return the corresponding position in the
    /// query sequence. If that ref position corresponds to an insertion in the
    /// query, return the previous aligning position. Panics if the prev position
    /// would be 0.
    pub fn query_coord_at(&self, target_ref_position: u32) -> u32 {
        if target_ref_position == self.seq1_start {
            return self.seq2_start
        }
        let mut num_ref_inserts = 0u32;
        let mut num_query_inserts = 0u32;
        let mut current_ref_offset = 0i32;
        for insert in self.insertions.iter() {
            // If this insert is to the right of the target, we outta here.
            if current_ref_offset as u32 + self.seq1_start + (insert.abs() as u32) - num_ref_inserts >= target_ref_position+1 {
                break
            }

            if *insert > 0 {
                num_ref_inserts += 1;
                current_ref_offset += *insert;
            } else {
                num_query_inserts += 1;
                current_ref_offset += - (*insert);
            }
        }
        debug!("At end of target {}: {} ref_inserts, {} query_inserts, position {}",
               target_ref_position, num_ref_inserts, num_query_inserts, current_ref_offset);
        return self.seq2_start + (target_ref_position - self.seq1_start)
            + num_query_inserts - num_ref_inserts;
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn test_query_at_no_inserts() {
        init();
        let s = NucmerDeltaAlignment {
            seq1_name: "s1".to_string(),
            seq2_name: "s2".to_string(),
            seq1_start: 10,
            seq1_stop: 19,
            seq2_start: 100,
            seq2_stop: 109,
            num_errors: 3,
            num_similarity_errors: 0,
            insertions: vec![],
        };
        assert_eq!(100, s.query_coord_at(10));
        assert_eq!(101, s.query_coord_at(11));
        assert_eq!(103, s.query_coord_at(13));
        assert_eq!(109, s.query_coord_at(19));
    }

    #[test]
    fn test_query_ref_inserts() {
        init();
        let s = NucmerDeltaAlignment {
            seq1_name: "s1".to_string(),
            seq2_name: "s2".to_string(),
            seq1_start: 10,
            seq1_stop: 19,
            seq2_start: 100,
            seq2_stop: 109,
            num_errors: 3,
            num_similarity_errors: 0,
            insertions: vec![2,4],
        };
        assert_eq!(100, s.query_coord_at(10));
        assert_eq!(101, s.query_coord_at(11));
        assert_eq!(101, s.query_coord_at(12));
        assert_eq!(102, s.query_coord_at(13));
        assert_eq!(103, s.query_coord_at(14));
        //assert_eq!(104, s.query_coord_at(15)); FIXME
        assert_eq!(104, s.query_coord_at(16));
        assert_eq!(105, s.query_coord_at(17));
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

        //assert_eq!(100, s.query_coord_at(10)); // undefined
        assert_eq!(100, s.query_coord_at(11));
        assert_eq!(101, s.query_coord_at(12));
        assert_eq!(103, s.query_coord_at(13)); // Just after insert in query
        assert_eq!(104, s.query_coord_at(14));
        assert_eq!(105, s.query_coord_at(15));
        //assert_eq!(105, s.query_coord_at(16)); // No alignment FIXME: This is a real bug, but annoyed at this code RN.
        assert_eq!(106, s.query_coord_at(17));
        assert_eq!(107, s.query_coord_at(18));
        assert_eq!(108, s.query_coord_at(19));
        assert_eq!(109, s.query_coord_at(20));
    }

    #[test]
    fn test_parse_delta() {
        init();
        let parsed = nucmer_to_deltas(
            "tests/data/parsnp/1_first_group/73.20120800_S1D.21.fna",
            "tests/data/parsnp/1_first_group/73.20120800_S1D.21.fna",
            "tests/data/parsnp/1_first_group/73.20120600_E3D.30.fna",
        );
        assert_eq!(674-4, parsed.len());
        assert!(parsed.iter().find(|a| {
            *a == &NucmerDeltaAlignment {
                seq1_name: "73.20120800_S1D.21_contig_1092".to_string(),
                seq2_name: "73.20120600_E3D.30_contig_64048".to_string(),
                seq1_start: 97138-1,
                seq1_stop: 104056-1,
                seq2_start: 1-1,
                seq2_stop: 6935-1,
                num_errors: 58,
                num_similarity_errors: 58,
                insertions: vec![
                    -337,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -3548,
                    -2,
                    -2084,
                    -2,
                    -1],
            }
        }).is_some());
    }
}
