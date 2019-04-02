use rand::prelude::*;

use rust_htslib::bam;
use rust_htslib::bam::Record;
use rust_htslib::bam::Read as BamRead;

use tempdir::TempDir;
use tempfile;

use bam_generator::*;

pub struct ReadSortedShardedBamReader {
    shard_bam_readers: Vec<bam::Reader>,
    previous_read_records: Option<Vec<bam::Record>>,
    next_record_to_return: Option<bam::Record>,
    winning_index: Option<usize>,
    tid_offsets: Vec<i32>,
}

impl ReadSortedShardedBamReader {
    // Return a Vec of read in BAM records, or None if there are no more.
    fn read_a_record_set(&mut self) -> Option<Vec<Record>> {
        let mut current_qname: Option<&[u8]> = None;
        let mut current_alignments: Vec<Record> = Vec::with_capacity(
            self.shard_bam_readers.len());
        let mut some_unfinished = false;
        let mut some_finished = false;
        for (i, reader) in self.shard_bam_readers.iter_mut().enumerate() {
            // loop until there is a primary alignment or the BAM file ends.
            loop {
                match reader.read(&mut current_alignments[i]) {
                    Ok(_) => {
                        let record = &current_alignments[i];
                        if !record.is_paired() {
                            panic!("This code can only handle paired-end \
                                    input (at the moment), sorry.")
                        }
                        if !record.is_secondary() && !record.is_supplementary() {
                            some_unfinished = true;
                            match current_qname {
                                None => {current_qname = Some(record.qname())},
                                Some(prev) => {
                                    if prev != record.qname() {
                                        panic!(
                                            "BAM files do not appear to be \
                                             properly sorted by read name. \
                                             Expected read name {:?} from a \
                                             previous reader but found {:?} \
                                             in the current.", prev, record.qname());
                                    }
                                }
                            }
                            break;
                        }
                        // else we have a non-primary alignment, so loop again.
                    },
                    Err(_) => {
                        debug!("BAM reader #{} appears to be finished", i);
                        some_finished = true;
                        break;
                    }
                }
            }
        }
        // check we are properly finished.
        if some_unfinished && some_finished {
            panic!("Unexpectedly one BAM file input finished while another had further reads")
        }
        return match some_unfinished {
            true => Some(current_alignments),
            false => None
        }
    }

    /// There is no way to copy a Record into a waiting one, as far as I know,
    /// so implement here. Details here copied from impl PartialEq for Record.
    fn clone_record_into(from: &Record, to: &mut Record) {
        to.set_pos(from.pos());
        to.set_bin(from.bin());
        to.set_mapq(from.mapq());
        to.set_flags(from.flags());
        to.set_mtid(from.mtid());
        to.set_mpos(from.mpos());
        to.set_insert_size(from.insert_size());
        //to.set_data(from.data()); // TODO: Needed?
    }
}

impl Iterator for ReadSortedShardedBamReader {
    type Item = bam::Record;

    fn next(&mut self) -> Option<bam::Record> {
        if self.next_record_to_return.is_some() {
            let mut to_return = bam::Record::new();
            {
                let record = self.next_record_to_return.as_ref().unwrap();
                ReadSortedShardedBamReader::clone_record_into(record, &mut to_return);
            }
            self.next_record_to_return = None;
            let tid_now = to_return.tid();
            to_return.set_tid(tid_now + self.tid_offsets[self.winning_index.unwrap()]);
            return Some(to_return);
        } else {
            if self.previous_read_records.is_none() {
                self.previous_read_records = self.read_a_record_set();
                if self.previous_read_records.is_none() {
                    // All finished all the input files.
                    return None;
                }
            }
            // Read the second set
            let second_read_alignments_opt = self.read_a_record_set();

            // If we get None, then croak
            if second_read_alignments_opt.is_none() {
                panic!("Unexpectedly was able to read a first read set, but not a second. Hmm.")
            }
            let second_read_alignments = second_read_alignments_opt.unwrap();

            // Decide which pair is the winner
            // Cannot use max_by_key() here since we want a random winner
            let mut max_score: Option<i64> = None;
            let mut winning_indices: Vec<usize> = vec![];
            for (i, aln1) in self.previous_read_records.unwrap().iter().enumerate() {
                let mut score: i64 = 0;
                score += aln1.aux(b"AS")
                    .expect(&format!(
                        "Record {:#?} (read1) unexpectedly did not have AS tag, which is needed for \
                         ranking pairs of alignments", aln1.qname())).integer();
                score += second_read_alignments[i].aux(b"AS")
                    .expect(&format!(
                        "Record {:#?} (read2) unexpectedly did not have AS tag, which is needed for \
                         ranking pairs of alignments", aln1.qname())).integer();
                if max_score.is_none() || score > max_score.unwrap() {
                    max_score = Some(score);
                    winning_indices = vec![i]
                } else if score == max_score.unwrap() {
                    winning_indices.push(i)
                }
                // Else a loser
            }
            let winning_index: usize = *winning_indices
                .choose(&mut thread_rng()).unwrap();

            // Set the next read to return
            self.winning_index = Some(winning_index);
            self.next_record_to_return = Some(second_read_alignments[winning_index]);

            let mut to_return = self.previous_read_records.unwrap()[winning_index];
            to_return.set_tid(to_return.tid() + self.tid_offsets[winning_index]);

            // Unset the current crop of first reads
            self.previous_read_records = None;

            // Return the first of the pair
            return Some(to_return);
        }
    }
}

struct ShardedBamReaderGenerator {
    stoit_name: String,
    read_sorted_bam_readers: Vec<bam::Reader>,
}

impl NamedBamReaderGenerator<ShardedBamReader> for ShardedBamReaderGenerator {
    fn start(self) -> ShardedBamReader {
        let mut new_header = bam::header::Header::new();
        let mut tid_offsets: Vec<i32> = vec!();

        // Read header info for each BAM file, and write to new BAM header,
        // adding tid offsets.
        let mut current_tid_offset: i32 = 0;
        for ref reader in self.read_sorted_bam_readers {
            let header = reader.header();
            let mut current_tid: u32 = 1; // tids start with 1.
            for name in header.target_names() {
                let length = header.target_len(current_tid)
                    .expect(&format!("Failed to get target length for TID {}", current_tid));
                // e.g. @SQ	SN:a62_bin.100.fna=k141_20475	LN:15123
                let mut current_record = bam::header::HeaderRecord::new(b"SQ");
                current_record.push_tag(
                    b"SN",
                    &std::str::from_utf8(name).unwrap());
                current_record.push_tag(b"LN", &length);
                new_header.push_record(&current_record);
            }

            tid_offsets.push(current_tid_offset);
            current_tid_offset += header.target_count() as i32;
        }

        // Instantiate and start Demux struct
        let demux = ReadSortedShardedBamReader {
            shard_bam_readers: self.read_sorted_bam_readers,
            tid_offsets: tid_offsets,
            previous_read_records: None,
            next_record_to_return: None,
            winning_index: None,
        };

        // Add samtools sort process
        let tmp_dir = TempDir::new("coverm_fifo")
            .expect("Unable to create samtools sort temporary directory");
        let sort_input_fifo_path = tmp_dir.path().join("sort_input.pipe");
        let sort_output_fifo_path = tmp_dir.path().join("sort_input.pipe");
        let sort_log_file = tempfile::NamedTempFile::new()
            .expect("Failed to create samtools sort log tempfile");

        // TODO: make threads available in function
        // TODO: Allow cache
        // TODO: sort in tempdir
        let sort_command_string = format!(
            "set -e -o pipefail; \
             samtools sort -l0 {} > {} 2> {}",
            sort_input_fifo_path.to_str()
                .expect("Failed to convert sort tempfile input path to str"),
            sort_output_fifo_path.to_str()
                .expect("Failed to convert sort tempfile output path to str"),
            sort_log_file.path().to_str()
                .expect("Failed to convert tempfile log path to str")
        );
        debug!("Running cmd_string: {}", sort_command_string);
        let mut cmd = std::process::Command::new("bash");
        cmd
            .arg("-c")
            .arg(&sort_command_string)
            .stderr(std::process::Stdio::piped());
        let sort_child = cmd.spawn().expect("Unable to execute bash");

        // Write all BAM records to the samtools sort input fifo. May as well
        // here because we have to wait later anyway. This could be put out into
        // a new thread but that is just overcomplicating things, for now.
        // TODO: Buffer the output here?
        let mut writer = bam::Writer::from_path(sort_input_fifo_path, &new_header)
            .expect("Failed to open BAM to write to samtools sort process");
        debug!("Writing records to samtools sort input FIFO..");
        for record in demux {
            writer.write(&record)
                .expect("Failed to write BAM record to samtools sort input fifo");
        }
        debug!("Finished writing records to samtools sort input FIFO.");

        return ShardedBamReader {
            stoit_name: self.stoit_name,
            bam_reader: bam::Reader::from_path(sort_output_fifo_path)
                .expect("Unable to open reader from samtools sort output FIFO"),
            tempdir: tmp_dir,
            sort_process: sort_child,
            sort_command_string: sort_command_string,
            sort_log_file_description: "samtools sort".to_string(),
            sort_log_file: sort_log_file,
        }
    }
}

struct ShardedBamReader {
    stoit_name: String,
    bam_reader: bam::Reader,
    tempdir: TempDir,
    sort_process: std::process::Child,
    sort_command_string: String,
    sort_log_file_description: String,
    sort_log_file: tempfile::NamedTempFile,
}

// Given a list of paths to different BAM files which are all mappings of the
// same read set to different references (all sorted by read name), generate a
// BAM reader that chooses the best place for each read to map to.
// pub fn generate_sharded_bam_reader_from_bam_files(
//     bam_paths: Vec<&str>) -> Vec<ShardedBamReaderGenerator> {
//     // open an output BAM file that gets put to samtools sort without -n
//     // For each BAM path,
//     // open the path
//     // record the number of references in the header

//     // write each header entry read in to a new reference entry in the output
//     // BAM, after adding the tid offset.
// }


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shard_hello_world() {
        let gen = ShardedBamReaderGenerator {
            stoit_name: "stoita".to_string(),
            read_sorted_bam_readers: vec![
                bam::Reader::from_path("test/data/2seqs.fastaVbad_read.bam").unwrap(),
                bam::Reader::from_path("test/data/7seqs.fnaVbad_read.bam").unwrap()
            ]
        };
        let _reader = gen.start();
    }
}
