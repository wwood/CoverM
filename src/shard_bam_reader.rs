use rand::prelude::*;


pub struct ReadSortedShardedBamReader {
    stoit_name: String,
    shard_bam_readers: Vec<bam::Reader>,
    previous_read_records: Option<Vec<Record>>,
    next_record_to_return: Option<Record>,
    winning_index: Option<u32>,
}

impl ReadSortedShardedBamReader {
    // Return a Vec of read in BAM records, or None if there are no more.
    fn read_a_record_set(&mut self) -> Option<Vec<Record>> {
        let mut current_qname: Option<&str> = None;
        let mut current_alignments: Vec<Record> = Vec::with_capacity(
            self.shard_bam_readers.len());
        let mut some_unfinished = false;
        let mut some_finished = false;
        for (i, reader) in self.shard_bam_reads.enumerate() {
            // loop until there is a primary alignment or the BAM file ends.
            loop {
                match reader.read(current_alignments[i]) {
                    Ok(_) => {
                        let record = current_alignments[i];
                        if !record.is_paired() {
                            panic!("This code can only handle paired-end \
                                    input (at the moment), sorry.")
                        }
                        if !record.is_secondary() && !record.is_supplementary() {
                            some_unfinished = true;
                            match current_qname {
                                None => {current_qname = Some(record.qname)},
                                Some(prev) => {
                                    if prev != record.qname {
                                        panic!(
                                            "BAM files do not appear to be \
                                             properly sorted by read name. \
                                             Expected read name {} from a \
                                             previous reader but found {} \
                                             in the current.", prev, record.qname);
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
        to.set_insert_size(from.insert_size);
        to.set_data(from.data());
    }
}

impl NamedBamReader for ReadSortedShardedBamReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }
    
    fn read(&mut self, record: &mut bam::record::Record) -> Result<(), bam::ReadError> {
        match self.next_record_to_return {
            Some(r) => {
                self.clone_record_into(r, record);
                record.set_tid(record.tid + self.tid_offsets[self.winning_index.unwrap()]);
                self.next_record_to_return = None;
            },
            None => {
                if self.previous_read_records.is_none() {
                    self.previous_read_records = self.read_a_record_set();
                    if self.previous_read_reads.is_none() {
                        // All finished all the input files.
                        return Err(bam::ReadError);
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
                let mut max_score: Option<u32> = None;
                let mut winning_indices: Vec<usize> = vec![];
                for (i, aln1) in self.previous_record_reads.unwrap().iter().enumerate() {
                    let mut score: u32 = 0;
                    score += match aln1.aux(b"AS")
                        .expect(format!(
                            "Record {} (read1) unexpectedly did not have AS tag, which is needed for \
                             ranking pairs of alignments", aln1.qname));
                    score += match second_read_alignments[i].aux(b"AS")
                        .expect(format!(
                            "Record {} (read2) unexpectedly did not have AS tag, which is needed for \
                             ranking pairs of alignments", aln1.qname));
                    if max_score.is_none() || score > max_score.unwrap() {
                        max_score = Some(score);
                        winning_indices = vec![i]
                    } else if score == max_score.unwrap() {
                        winning_indices.push(i)
                    }
                    // Else a loser
                }
                let winning_index = winning_indices.choose(&mut thread_rng());

                // Set the next read to return
                self.winning_index = Some(winning_index);
                self.next_record_to_return = Some(second_read_alignments[winning_index]);

                let mut to_return = self.previous_record_reads[winning_index];
                to_return.set_tid(to_return.tid + self.tid_offsets[winning_index]);

                // Unset the current crop of first reads
                self.previous_record_reads = None;

                // Return the first of the pair
                self.clone_record_into(to_return, record);
                record.set_tid(record.tid + self.tid_offsets[winning_indices]);
            }
        }
    }

    fn header(&self) -> &bam::HeaderView {
        self.sorted_bam_reader.header()
    }
    fn finish(self) {;}

    fn num_detected_primary_alignments(&self) -> u64 {
        return self.num_detected_primary_alignments
    }
}

struct ShardedBamReaderGenerator {
    stoit_name: String,
    tempdir: TempDir,
    fifo_path: std::path::PathBuf,
    pre_processes: Vec<std::process::Command>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
}

// Given a list of paths to different BAM files which are all mappings of the
// same read set to different references (all sorted by read name), generate a
// BAM reader that chooses the best place for each read to map to.
pub fn generate_sharded_bam_reader(
    bam_paths: Vec<&str>) -> ShardedBamReader {
    // open an output BAM file that gets put to samtools sort without -n
    // For each BAM path,
    // open the path
    // record the number of references in the header

    // write each header entry read in to a new reference entry in the output
    // BAM, after adding the tid offset.
}
