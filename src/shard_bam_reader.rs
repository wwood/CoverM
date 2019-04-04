use std::str;

use rand::prelude::*;

use rust_htslib::bam;
use rust_htslib::bam::Record;
use rust_htslib::bam::Read as BamRead;

use tempdir::TempDir;
use tempfile;
use nix::unistd;
use nix::sys::stat;


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
        let mut current_alignments: Vec<Record> = Vec::with_capacity(
            self.shard_bam_readers.len());
        let mut _current_qname: Option<&[u8]> = None;
        let mut some_unfinished = false;
        let mut some_finished = false;
        for (i, reader) in self.shard_bam_readers.iter_mut().enumerate() {
            // loop until there is a primary alignment or the BAM file ends.
            loop {
                {
                    current_alignments.push(bam::Record::new());
                    let res = reader.read(&mut current_alignments[i]);
                    if res.is_err() {
                        debug!("BAM reader #{} appears to be finished", i);
                        some_finished = true;
                        break;
                    }
                }

                {
                    let record = &current_alignments[i];
                    if !record.is_paired() {
                        panic!("This code can only handle paired-end \
                                input (at the moment), sorry.")
                    }
                    if !record.is_secondary() && !record.is_supplementary() {
                        some_unfinished = true;
                        // TODO: uncomment below and fix borrow checker issues
                        // match current_qname {
                        //     None => {current_qname = Some(record.qname().clone())},
                        //     Some(prev) => {
                        //         if prev != record.qname() {
                        //             panic!(
                        //                 "BAM files do not appear to be \
                        //                  properly sorted by read name. \
                        //                  Expected read name {:?} from a \
                        //                  previous reader but found {:?} \
                        //                  in the current.", prev, record.qname());
                        //         }
                        //     }
                        // }
                        break;
                    }
                    // else we have a non-primary alignment, so loop again.
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

    // There is no way to copy a Record into a waiting one, as far as I know,
    // so implement here. Details here copied from impl PartialEq for Record.
    fn clone_record_into(from: &Record, to: &mut Record) {
        to.set_pos(from.pos());
        to.set_bin(from.bin());
        to.set_mapq(from.mapq());
        to.set_flags(from.flags());
        to.set_mtid(from.mtid());
        to.set_mpos(from.mpos());
        to.set_insert_size(from.insert_size());
        to.set_qname(from.qname());
        to.set_tid(from.tid());
    }
}

impl ReadSortedShardedBamReader {
    fn read(&mut self, to_return: &mut bam::Record) -> Result<(), bam::ReadError> {
        if self.next_record_to_return.is_some() {
            {
                let record = self.next_record_to_return.as_ref().unwrap();
                ReadSortedShardedBamReader::clone_record_into(record, to_return);
            }
            self.next_record_to_return = None;
            let tid_now = to_return.tid();
            to_return.set_tid(tid_now + self.tid_offsets[self.winning_index.unwrap()]);
            return Ok(());
        } else {
            if self.previous_read_records.is_none() {
                self.previous_read_records = self.read_a_record_set();
                if self.previous_read_records.is_none() {
                    // All finished all the input files.
                    return Err(bam::ReadError::NoMoreRecord);
                }
            }
            // Read the second set
            // If we get None, then croak
            let second_read_alignments = self.read_a_record_set()
                .expect("Unexpectedly was able to read a first read set, but not a second. Hmm.");

            // Decide which pair is the winner
            // Cannot use max_by_key() here since we want a random winner
            let mut max_score: Option<i64> = None;
            let mut winning_indices: Vec<usize> = vec![];
            match self.previous_read_records {
                None => unreachable!(),
                Some(ref previous_records) => {
                    for (i, ref aln1) in previous_records.iter().enumerate() {
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
                }
            };

            let winning_index: usize = *winning_indices
                .choose(&mut thread_rng()).unwrap();
            debug!("Choosing winning index {} from winner pool {:?}",
                   winning_index, winning_indices);

            // Set the next read to return
            self.winning_index = Some(winning_index);
            self.next_record_to_return = Some(
                second_read_alignments[winning_index].clone());


            match self.previous_read_records {
                None => unreachable!(),
                Some(ref prev) => {
                    debug!("previous tid {}", &prev[winning_index].tid());
                    ReadSortedShardedBamReader::clone_record_into(
                        &prev[winning_index], to_return)
                }
            };

            let tid_now = to_return.tid();
            to_return.set_tid(tid_now + self.tid_offsets[winning_index]);
            debug!("Reindexed TID is {}, tid_offsets are {:?} and tid_now is {}",
                   to_return.tid(), self.tid_offsets, tid_now);

            // Unset the current crop of first reads
            self.previous_read_records = None;

            // TODO: Remove this debug code
            // let mut reader2 = bam::Reader::from_path(&"2seqs.fastaVbad_read.bam").unwrap();
            // reader2.read(to_return);

            // Return the first of the pair
            return Ok(());
        }
    }
}

pub struct ShardedBamReaderGenerator {
    pub stoit_name: String,
    pub read_sorted_bam_readers: Vec<bam::Reader>,
}

impl NamedBamReaderGenerator<ShardedBamReader> for ShardedBamReaderGenerator {
    fn start(self) -> ShardedBamReader {
        let mut new_header = bam::header::Header::new();
        let mut tid_offsets: Vec<i32> = vec!();

        // Read header info for each BAM file, and write to new BAM header,
        // adding tid offsets.
        let mut current_tid_offset: i32 = 0;
        for ref reader in &self.read_sorted_bam_readers {
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

        let tmp_dir = TempDir::new("coverm_fifo")
            .expect("Unable to create samtools sort temporary directory");
        // let tmp_dir = std::path::Path::new("/tmp/miner");
        // std::fs::create_dir(tmp_dir).expect("Failed to make dummy dir");
        let sort_input_fifo_path = tmp_dir.path().join("sort_input.pipe");
        let sort_output_fifo_path = tmp_dir.path().join("sort_output.pipe");
        // let sort_input_fifo_path = tmp_dir.join("sort_input.pipe");
        // let sort_output_fifo_path = tmp_dir.join("sort_output.pipe");
        let sort_log_file = tempfile::NamedTempFile::new()
            .expect("Failed to create samtools sort log tempfile");

        // create new fifo and give read, write and execute rights to the owner.
        // This is required because we cannot open a Rust stream as a BAM file
        // with rust-htslib.
        debug!("Creating FIFOs ..");
        unistd::mkfifo(&sort_input_fifo_path, stat::Mode::S_IRWXU)
            .expect(&format!("Error creating named pipe {:?}", sort_input_fifo_path));
        unistd::mkfifo(&sort_output_fifo_path, stat::Mode::S_IRWXU)
            .expect(&format!("Error creating named pipe {:?}", sort_output_fifo_path));


        // Instantiate and start Demux struct
        debug!("Instantiating ReadSortedShardedBamReader ..");
        let mut demux = ReadSortedShardedBamReader {
            shard_bam_readers: self.read_sorted_bam_readers,
            tid_offsets: tid_offsets,
            previous_read_records: None,
            next_record_to_return: None,
            winning_index: None,
        };




        // Start reader in a different thread because it needs to be running
        // when the sort starts spitting out results, but won't return from
        // instantiation until it has read the header (I think).
        let sort_output_fifo_path2 = sort_output_fifo_path.clone();
        let sorted_reader_join_handle: std::thread::JoinHandle<bam::Reader> = std::thread::spawn(
            move || {
                debug!("Starting to open sorted read BAM ..");
                let reader = bam::Reader::from_path(&sort_output_fifo_path2)
                    .expect("Unable to open reader from samtools sort output FIFO");
                debug!("Finished opening reader to samtools sort output FIFO..");
                reader
            }
        );

        // TODO: make threads available in function
        // TODO: Allow cache
        // TODO: sort in tempdir
        let sort_command_string = format!(
            "set -e -o pipefail; \
             samtools sort -l0 {} -o {}",
            sort_input_fifo_path.to_str()
                .expect("Failed to convert sort tempfile input path to str"),
            sort_output_fifo_path.to_str()
                .expect("Failed to convert sort tempfile output path to str"),
            //sort_log_file.path().to_str()
             //   .expect("Failed to convert tempfile log path to str")
        );
        debug!("Running cmd_string: {}", sort_command_string);
        let mut cmd = std::process::Command::new("bash");
        cmd
            .arg("-c")
            .arg(&sort_command_string);
            //.stderr(std::process::Stdio::piped())
        let sort_child = cmd.spawn().expect("Unable to execute bash");


        // Write all BAM records to the samtools sort input fifo. May as well
        // here because we have to wait later anyway. This could be put out into
        // a new thread but that is just overcomplicating things, for now.
        // TODO: Buffer the output here?
        {
            // Write in a scope so writer drops before we start reading from the sort.
            debug!("Opening BAM writer to {}", &sort_input_fifo_path.to_str().unwrap());
            let mut writer = bam::Writer::from_path(&sort_input_fifo_path, &new_header)
                .expect("Failed to open BAM to write to samtools sort process");
            debug!("Writing records to samtools sort input FIFO..");
            let mut record = bam::Record::new();
            while demux.read(&mut record).is_ok() {
                debug!("Writing tid {} for qname {}", record.tid(), str::from_utf8(record.qname()).unwrap());
                writer.write(&record)
                    .expect("Failed to write BAM record to samtools sort input fifo");
            }
            debug!("Finished writing records to samtools sort input FIFO.");
        }

        let reader = sorted_reader_join_handle.join()
            .expect("sorted reader thread failed");

        return ShardedBamReader {
            stoit_name: self.stoit_name,
            bam_reader: reader,
            //tempdir: tmp_dir,  // FIXME: remove after debug
            sort_process: sort_child,
            sort_command_string: sort_command_string,
            sort_log_file_description: "samtools sort".to_string(),
            sort_log_file: sort_log_file,
        }
    }
}

pub struct ShardedBamReader {
    pub stoit_name: String,
    pub bam_reader: bam::Reader,
    //tempdir: TempDir,
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
                bam::Reader::from_path("tests/data/2seqs.fastaVbad_read.bam").unwrap(),
                bam::Reader::from_path("tests/data/7seqs.fnaVbad_read.bam").unwrap()
            ]
        };
        let mut reader = gen.start();
        assert_eq!("stoita".to_string(), reader.stoit_name);
        let mut r = bam::Record::new();
        reader.bam_reader.read(&mut r);
        println!("{}",str::from_utf8(r.qname()).unwrap());
        println!("{}",r.tid());

        assert_eq!(1,2);
    }
}
