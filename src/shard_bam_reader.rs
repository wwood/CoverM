use std::str;
use std;
use std::process;

use rand::prelude::*;

use rust_htslib::bam;
use rust_htslib::bam::Record;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::Read as BamRead;
use rust_htslib::bam::errors::Result as HtslibResult;

use mapping_parameters::ReadFormat;

use tempdir::TempDir;
use tempfile;
use nix::unistd;
use nix::sys::stat;

use bam_generator::*;
use bam_generator::complete_processes;
use genome_exclusion::*;

use std::slice;
use std::ffi;

// Like Header#names() except just return one name
fn bam_header_target_name(header: &bam::HeaderView, tid: usize) -> &[u8] {
    let names = unsafe {
        slice::from_raw_parts(header.inner().target_name, header.target_count() as usize)
    };
    return unsafe {
        ffi::CStr::from_ptr(names[tid]).to_bytes()
    }
}

pub struct ReadSortedShardedBamReader<'a, T>
where T: GenomeExclusion {
    shard_bam_readers: Vec<bam::Reader>,
    previous_read_records: Option<Vec<bam::Record>>,
    next_record_to_return: Option<bam::Record>,
    winning_index: Option<usize>,
    tid_offsets: Vec<i32>,
    genome_exclusion: &'a T,
}

impl<'a, T> ReadSortedShardedBamReader<'a, T>
where T: GenomeExclusion {
    // Return a Vec of read in BAM records, or None if there are no more.
    fn read_a_record_set(&mut self) -> Option<Vec<Record>> {
        let mut current_alignments: Vec<Record> = Vec::with_capacity(
            self.shard_bam_readers.len());
        let mut current_qname: Option<String> = None;
        let mut some_unfinished = false;
        let mut some_finished = false;
        let mut record;
        let mut res;
        let mut current_alignment;
        for (i, reader) in self.shard_bam_readers.iter_mut().enumerate() {
            // loop until there is a primary alignment or the BAM file ends.
            loop {
                {
                    current_alignment = bam::Record::new();
                    res = reader.read(&mut current_alignment);
                    if res.expect("Failure to read from a shard BAM file") == false {
                        debug!("BAM reader #{} appears to be finished", i);
                        some_finished = true;
                        break;
                    }
                }

                {
                    record = current_alignment.clone();
                    if !record.is_paired() {
                        error!("This code can only handle paired-end \
                                input (at the moment), sorry. Found \
                                record {:?}",
                                record);
                        process::exit(1);
                    }
                    if !record.is_secondary() && !record.is_supplementary() {
                        some_unfinished = true;
                         match current_qname.clone() {
                             None => {
                                 current_qname = Some(
                                     str::from_utf8(record.qname())
                                         .unwrap()
                                         .to_string())},
                             Some(prev) => {
                                 if prev != str::from_utf8(record.qname())
                                     .unwrap()
                                     .to_string() {
                                         error!(
                                             "BAM files do not appear to be \
                                              properly sorted by read name. \
                                              Expected read name {:?} from a \
                                              previous reader but found {:?} \
                                              in the current.",
                                             prev,
                                             str::from_utf8(
                                                 record.qname()).unwrap().to_string());
                                         process::exit(1);
                                 }
                             }
                         }
                        current_alignments.push(current_alignment);
                        break;
                    }
                    // else we have a non-primary alignment, so loop again.
                }
            }
        }
        // check we are properly finished.
        if some_unfinished && some_finished {
            error!("Unexpectedly one BAM file input finished while another had further reads");
            process::exit(1);
        }
        debug!("Read records {:?}", current_alignments);
        return match some_unfinished {
            true => Some(current_alignments),
            false => None
        }
    }

    // There is no way to copy a Record into a waiting one, as far as I know,
    // so implement here. Details here copied from impl PartialEq for Record.
    fn clone_record_into(from: &Record, to: &mut Record) {
        //Clone record using set() also add cigar and push aux tags NM
        to.set(from.qname(),
               Some(&CigarString::from_str(from.cigar().to_string().as_str()).unwrap()),
               from.seq().as_bytes().as_slice(),
               from.qual());
        to.set_pos(from.pos());
        to.set_bin(from.bin());
        to.set_mapq(from.mapq());
        to.set_flags(from.flags());
        to.set_mtid(from.mtid());
        to.set_mpos(from.mpos());
        to.set_insert_size(from.insert_size());
        to.set_tid(from.tid());
        match &from.aux("NM".as_bytes()) {
            Some(aux) => {
                to.push_aux("NM".as_bytes(), aux);
            },
            None => {
                if from.tid() >= 0 && from.cigar_len() != 0 {
                    error!("record {:?} with name {} had no NM tag",
                           from, str::from_utf8(from.qname()).unwrap());
                    process::exit(1);
                }
            },
        }
        debug!("from cigar {:?}, to cigar {:?}", from.cigar().to_string(), to.cigar().to_string())
    }
}

impl<'a, T> ReadSortedShardedBamReader<'a, T>
where T: GenomeExclusion {
    fn read(&mut self, to_return: &mut bam::Record) -> HtslibResult<bool> {
        if self.next_record_to_return.is_some() {
            {
                let record = self.next_record_to_return.as_ref().unwrap();
                ReadSortedShardedBamReader::<T>::clone_record_into(record, to_return);
            }
            self.next_record_to_return = None;
            let tid_now = to_return.tid();
            to_return.set_tid(tid_now + self.tid_offsets[self.winning_index.unwrap()]);
            return Ok(true);
        } else {
            if self.previous_read_records.is_none() {
                self.previous_read_records = self.read_a_record_set();
                if self.previous_read_records.is_none() {
                    // All finished all the input files.
                    return Ok(false);
                }
            }
            // Read the second set
            // If we get None, then croak
            let second_read_alignments = self.read_a_record_set()
                .expect("Unexpectedly was able to read a first read set, but not a second. Hmm.");

            debug!("Previous records {:?}", self.previous_read_records);
            debug!("Second read records {:?}", second_read_alignments);

            // Decide which pair is the winner
            // Cannot use max_by_key() here since we want a random winner
            let mut max_score: Option<i64> = None;
            let mut winning_indices: Vec<usize> = vec![];
            match self.previous_read_records {
                None => unreachable!(),
                Some(ref previous_records) => {
                    for (i, ref aln1) in previous_records.iter().enumerate() {
                        let tid = aln1.tid();
                        if tid < 0 || !self.genome_exclusion.is_excluded(
                            bam_header_target_name(&self.shard_bam_readers[i].header(), tid as usize)) {
                            let mut score: i64 = 0;
                            // Unlike BWA-MEM, Minimap2 does not have AS tags
                            // when the read is unmapped.
                            if !aln1.is_unmapped() {
                                score += aln1.aux(b"AS")
                                    .expect(&format!(
                                        "Record {:#?} (read1) unexpectedly did not have AS tag, which is needed for \
                                        ranking pairs of alignments", aln1.qname())).integer();
                            }
                            if !second_read_alignments[i].is_unmapped() {
                                score += second_read_alignments[i].aux(b"AS")
                                    .expect(&format!(
                                        "Record {:#?} (read2) unexpectedly did not have AS tag, which is needed for \
                                        ranking pairs of alignments", aln1.qname())).integer();
                            }
                            if max_score.is_none() || score > max_score.unwrap() {
                                max_score = Some(score);
                                winning_indices = vec![i]
                            } else if score == max_score.unwrap() {
                                winning_indices.push(i)
                            }
                            // Else a loser when there was a previous winner
                        }
                    }
                }
            };

            let winning_index: usize;
            if winning_indices.len() > 1 {
                winning_index = *winning_indices
                    .choose(&mut thread_rng()).unwrap();
            } else if winning_indices.len() == 1 {
                winning_index = winning_indices[0];
            } else {
                error!("CoverM cannot currently deal with reads that only map to excluded genomes");
                process::exit(1);
            }
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
                    ReadSortedShardedBamReader::<T>::clone_record_into(
                        &prev[winning_index], to_return)
                }
            };

            let tid_now = to_return.tid();
            to_return.set_tid(tid_now + self.tid_offsets[winning_index]);
            debug!("Reindexed TID is {}, tid_offsets are {:?} and tid_now is {}",
                   to_return.tid(), self.tid_offsets, tid_now);

            // Unset the current crop of first reads
            self.previous_read_records = None;

            // Return the first of the pair
            return Ok(true);
        }
    }
}

pub struct ShardedBamReaderGenerator<'a, T>
where T: GenomeExclusion {
    pub stoit_name: String,
    pub read_sorted_bam_readers: Vec<bam::Reader>,
    pub sort_threads: i32,
    pub genome_exclusion: &'a T,
}

impl<'a, T> NamedBamReaderGenerator<ShardedBamReader> for ShardedBamReaderGenerator<'a, T>
where T: GenomeExclusion {
    fn start(self) -> ShardedBamReader {
        let mut new_header = bam::header::Header::new();
        let mut tid_offsets: Vec<i32> = vec!();

        // Read header info for each BAM file, and write to new BAM header,
        // adding tid offsets.
        let mut current_tid_offset: i32 = 0;
        for reader in self.read_sorted_bam_readers.iter() {
            let header = reader.header();
            let mut current_tid: u32 = 0; // I think TID counting should start at 0, 1 returns wrong lengths.
            let names = header.target_names();
            for name in names.iter() {
                let length = header.target_len(current_tid)
                    .expect(&format!("Failed to get target length for TID {}", current_tid));
                // e.g. @SQ	SN:a62_bin.100.fna=k141_20475	LN:15123
                let mut current_record = bam::header::HeaderRecord::new(b"SQ");
                current_record.push_tag(
                    b"SN",
                    &std::str::from_utf8(name).unwrap());
                current_record.push_tag(b"LN", &length);
                new_header.push_record(&current_record);
                current_tid += 1; // increment current tid
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
        let sort_temp_fifo_path = tmp_dir.path().join("sort_temp.pipe");
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
            genome_exclusion: self.genome_exclusion,
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
        // Threads now available in function
        // TODO: Allow cache
        // Temp files now sorted and stored in tempdir
        let sort_command_string = format!(
            "set -e -o pipefail; \
             samtools sort -l0 {} -o {} -T {} -@ {}",
            sort_input_fifo_path.to_str()
                .expect("Failed to convert sort tempfile input path to str"),
            sort_output_fifo_path.to_str()
                .expect("Failed to convert sort tempfile output path to str"),
            sort_temp_fifo_path.to_str()
                .expect("Failed to convert sort tempfile tempfile path to str"),
            &self.sort_threads,
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
        // let mut output_buffer = bam::buffer::RecordBuffer::new(reader); Like this?
        {
            // Write in a scope so writer drops before we start reading from the sort.
            debug!("Opening BAM writer to {}", &sort_input_fifo_path.to_str().unwrap());
            let mut writer = bam::Writer::from_path(
                &sort_input_fifo_path, 
                &new_header,
                rust_htslib::bam::Format::BAM)
                .expect("Failed to open BAM to write to samtools sort process");
            // Do not compress since these records are just read back in again -
            // compression would be wasteful of CPU (and IO with samtools sort)
            writer.set_compression_level(bam::CompressionLevel::Uncompressed)
                .expect("Failure to set BAM writer compression level - programming bug?");
            debug!("Writing records to samtools sort input FIFO..");
            let mut record = bam::Record::new();
            while demux.read(&mut record) == Ok(true) {
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
            tempdir: tmp_dir,
            sort_process: sort_child,
            sort_command_string: sort_command_string,
            sort_log_file_description: "samtools sort".to_string(),
            sort_log_file: sort_log_file,
            num_detected_primary_alignments: 0,
        }
    }
}

pub struct ShardedBamReader {
    stoit_name: String,
    bam_reader: bam::Reader,
    tempdir: TempDir,
    sort_process: std::process::Child,
    sort_command_string: String,
    sort_log_file_description: String,
    sort_log_file: tempfile::NamedTempFile,
    num_detected_primary_alignments: u64,
}

impl NamedBamReader for ShardedBamReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }
    fn read(&mut self, record: &mut bam::record::Record) -> HtslibResult<bool> {
        let res = self.bam_reader.read(record);
        if res == Ok(true) && !record.is_secondary() && !record.is_supplementary() {
            self.num_detected_primary_alignments += 1;
        }
        return res;
    }
    fn header(&self) -> &bam::HeaderView {
        &self.bam_reader.header()
    }
    fn finish(self) {
        complete_processes(
            vec![self.sort_process],
            vec![self.sort_command_string],
            vec![self.sort_log_file_description],
            vec![self.sort_log_file],
            Some(self.tempdir));
    }

    fn set_threads(&mut self, n_threads: usize) {
        self.bam_reader.set_threads(n_threads).unwrap();
    }
    fn num_detected_primary_alignments(&self) -> u64 {
        return self.num_detected_primary_alignments
    }
}

// Given a list of paths to different BAM files which are all mappings of the
// same read set to different references (all sorted by read name), generate a
// BAM reader that chooses the best place for each read to map to.
pub fn generate_sharded_bam_reader_from_bam_files<'a, T>(
    bam_paths: Vec<&str>, sort_threads: i32, genome_exclusion: &'a T)
    -> Vec<ShardedBamReaderGenerator<'a, T>>
where T: GenomeExclusion {
    // open an output BAM file that gets put to samtools sort without -n
    // For each BAM path,
    // open the path
    // record the number of references in the header
    // write each header entry read in to a new reference entry in the output
    // BAM, after adding the tid offset.
    let bam_readers = bam_paths.iter().map(
        |f| {
            debug!("Opening BAM {} ..", f);
            bam::Reader::from_path(&f)
                .expect(&format!("Unable to open bam file {}", f))
        }
    ).collect();
    let stoit_name = bam_paths.iter().map(
        |f| std::path::Path::new(f).file_stem().unwrap().to_str().expect(
            "failure to convert bam file name to stoit name - UTF8 error maybe?").to_string())
        .fold(
            None,
            { |acc, s| match acc {
                None => Some(s),
                Some(prev) => Some(format!("{}|{}",prev,s))
            }}
    ).unwrap();
    debug!("Opened all input BAM files");
    let gen = ShardedBamReaderGenerator {
        stoit_name: stoit_name,
        read_sorted_bam_readers: bam_readers,
        sort_threads: sort_threads,
        genome_exclusion: genome_exclusion,
    };
    return vec!(gen);

}

pub fn generate_named_sharded_bam_readers_from_reads(
    mapping_program: MappingProgram,
    reference: &str,
    read1_path: &str,
    read2_path: Option<&str>,
    read_format: ReadFormat,
    threads: u16,
    cached_bam_file: Option<&str>,
    discard_unmapped: bool,
    mapping_options: Option<&str>) -> bam::Reader {

    let tmp_dir = TempDir::new("coverm_fifo")
        .expect("Unable to create temporary directory");

    let fifo_path = tmp_dir.path().join("foo.pipe");

    // create new fifo and give read, write and execute rights to the owner.
    // This is required because we cannot open a Rust stream as a BAM file with
    // rust-htslib.
    unistd::mkfifo(&fifo_path, stat::Mode::S_IRWXU)
        .expect(&format!("Error creating named pipe {:?}", fifo_path));

    let mapping_log = tempfile::NamedTempFile::new()
        .expect(&format!("Failed to create {:?} log tempfile", mapping_program));
    let samtools2_log = tempfile::NamedTempFile::new()
        .expect("Failed to create second samtools log tempfile");
    // tempfile does not need to be created but easier to create than get around
    // borrow checker.
    let samtools_view_cache_log = tempfile::NamedTempFile::new()
        .expect("Failed to create cache samtools view log tempfile");

    let cached_bam_file_args = match cached_bam_file {
        Some(path) => {
            format!(
                "|tee {:?} |samtools view {} -t {} -b -o '{}' 2>{}",
                // tee
                fifo_path,
                // samtools view
                // cannot discard unmapped for sharded bam files
                match discard_unmapped {
                    _ => "",
                },
                threads,
                path,
                samtools_view_cache_log.path().to_str()
                    .expect("Failed to convert tempfile path to str"))
        },
        None => format!("> {:?}", fifo_path)
    };

    let mapping_command = build_mapping_command(
        mapping_program, 
        read_format,
        threads,
        read1_path,
        reference,
        read2_path,
        mapping_options);

    let bwa_sort_prefix = tempfile::Builder::new()
        .prefix("coverm-make-samtools-sort")
        .tempfile_in(tmp_dir.path())
        .expect("Failed to create tempfile as samtools sort prefix");
    let cmd_string = format!(
        "set -e -o pipefail; \
         {} 2>{} \
         | samtools sort -n -T '{}' -l0 -@ {} 2>{} \
         {}",
        // Mapping
        mapping_command,
        mapping_log.path().to_str().expect("Failed to convert tempfile path to str"),
        // samtools
        bwa_sort_prefix.path().to_str()
            .expect("Failed to convert bwa_sort_prefix tempfile to str"),
        threads - 1,
        samtools2_log.path().to_str().expect("Failed to convert tempfile path to str"),
        // Caching (or not)
        cached_bam_file_args);
    debug!("Queuing cmd_string: {}", cmd_string);
    let mut cmd = std::process::Command::new("bash");
    cmd
        .arg("-c")
        .arg(&cmd_string)
        .stderr(std::process::Stdio::piped());

    let mut log_descriptions = vec![
        format!("{:?}",mapping_program).to_string(),
        "samtools sort".to_string()];
    let mut log_files = vec![mapping_log, samtools2_log];
    if cached_bam_file.is_some() {
        log_descriptions.push("samtools view for cache".to_string());
        log_files.push(samtools_view_cache_log);
    }

    debug!("Starting mapping processes");
    let pre_processes = vec![cmd];
    let command_strings = vec![format!("bash -c {}", cmd_string)];
    let mut processes = vec![];
    let mut i = 0;
    for mut preprocess in pre_processes {
        debug!("Running mapping command: {}", command_strings[i]);
        i += 1;
        processes.push(preprocess
            .spawn()
            .expect("Unable to execute bash"));
    }
    return bam::Reader::from_path(&fifo_path)
        .expect(&format!("Unable to open bam file {:?}", &fifo_path))
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shard_hello_world() {
        //This test needs to be revisited. It seems busted.
        let gen = ShardedBamReaderGenerator {
            stoit_name: "stoiter".to_string(),
            read_sorted_bam_readers: vec![
                bam::Reader::from_path("tests/data/2seqs.fastaVbad_read.bam").unwrap(),
                bam::Reader::from_path("tests/data/7seqs.fnaVbad_read.bam").unwrap()
            ],
            sort_threads: 1,
            genome_exclusion: &NoExclusionGenomeFilter{},
        };
        let mut reader = gen.start();
        assert_eq!("stoiter".to_string(), reader.stoit_name);
        let mut r = bam::Record::new();
        reader.bam_reader.read(&mut r).unwrap();
        println!("{}",str::from_utf8(r.qname()).unwrap());
        println!("{}",r.tid());

//        assert_eq!(1,2); // Not sure about this, might be because I changed the the starting TID to be 0 rather than 1
    }
}
