use std;
use std::io::Read;

use filter::*;
use bwa_index_maintenance::BwaIndexStruct;
use mapping_parameters::ReadFormat;
use FlagFilter;

use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

use nix::unistd;
use nix::sys::stat;
use tempdir::TempDir;
use tempfile;

pub trait NamedBamReader {
    // Name of the stoit
    fn name(&self) -> &str;

    // Read a record into record parameter
    fn read(&mut self, record: &mut bam::record::Record) -> Result<(), bam::ReadError>;

    // Return the bam header of the final BAM file
    fn header(&self) -> &bam::HeaderView;

    fn finish(self);

    // Number of reads that were detected
    fn num_detected_primary_alignments(&self) -> u64;
}

pub trait NamedBamReaderGenerator<T> {
    // For readers that map, start the process of mapping
    fn start(self) -> T;
}

pub struct BamFileNamedReader {
    stoit_name: String,
    bam_reader: bam::Reader,
    num_detected_primary_alignments: u64,
}

impl NamedBamReader for BamFileNamedReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }
    fn read(&mut self, record: &mut bam::record::Record) -> Result<(), bam::ReadError> {
        let res = self.bam_reader.read(record);
        if res.is_ok() && !record.is_secondary() && !record.is_supplementary() {
            self.num_detected_primary_alignments += 1;
        }
        return res;
    }
    fn header(&self) -> &bam::HeaderView {
        self.bam_reader.header()
    }
    fn finish(self) {;}

    fn num_detected_primary_alignments(&self) -> u64 {
        return self.num_detected_primary_alignments
    }
}

impl NamedBamReaderGenerator<BamFileNamedReader> for BamFileNamedReader {
    fn start(self) -> BamFileNamedReader {
        BamFileNamedReader {
            stoit_name: self.stoit_name,
            bam_reader: self.bam_reader,
            num_detected_primary_alignments: 0,
        }
    }
}

pub struct StreamingNamedBamReader {
    stoit_name: String,
    bam_reader: bam::Reader,
    tempdir: TempDir,
    processes: Vec<std::process::Child>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
    num_detected_primary_alignments: u64,
}

pub struct StreamingNamedBamReaderGenerator {
    stoit_name: String,
    tempdir: TempDir,
    fifo_path: std::path::PathBuf,
    pre_processes: Vec<std::process::Command>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
}

impl NamedBamReaderGenerator<StreamingNamedBamReader> for StreamingNamedBamReaderGenerator {
    fn start(self) -> StreamingNamedBamReader {
        debug!("Starting mapping processes");
        let mut processes = vec![];
        let mut i = 0;
        for mut preprocess in self.pre_processes {
            debug!("Running mapping command: {}", self.command_strings[i]);
            i += 1;
            processes.push(preprocess
                .spawn()
                .expect("Unable to execute bash"));
        }
        let bam_reader = bam::Reader::from_path(&self.fifo_path)
            .expect(&format!("Unable to find BAM file {:?}", self.fifo_path));
        return StreamingNamedBamReader {
            stoit_name: self.stoit_name,
            bam_reader: bam_reader,
            tempdir: self.tempdir,
            processes: processes,
            command_strings: self.command_strings,
            log_file_descriptions: self.log_file_descriptions,
            log_files: self.log_files,
            num_detected_primary_alignments: 0,
        }
    }
}

pub fn complete_processes(
    processes: Vec<std::process::Child>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
    tempdir: Option<TempDir>) {

    for mut process in processes {
        let es = process.wait().expect("Failed to glean exitstatus from mapping process");
        if !es.success() {
            error!("Error when running mapping process: {:?}", command_strings);
            let mut err = String::new();
            process.stderr.expect("Failed to grab stderr from failed mapping process")
                .read_to_string(&mut err).expect("Failed to read stderr into string");
            error!("The overall STDERR was: {:?}", err);
            for (description, tf) in log_file_descriptions.into_iter().zip(
                log_files.into_iter()) {
                let mut contents = String::new();
                tf.into_file().read_to_string(&mut contents)
                    .expect(&format!("Failed to read log file for {}", description));
                error!("The STDERR for the {:} part was: {}",
                       description, contents);
            }
            panic!("Cannot continue since mapping failed.");
        }
    }
    // Close tempdir explicitly. Maybe not needed.
    match tempdir {
        Some(td) => {
            td.close().expect("Failed to close tempdir");
        },
        None => {}
    }
}

impl NamedBamReader for StreamingNamedBamReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }
    fn read(&mut self, record: &mut bam::record::Record) -> Result<(), bam::ReadError> {
        let res = self.bam_reader.read(record);
        if res.is_ok() && !record.is_secondary() && !record.is_supplementary() {
            self.num_detected_primary_alignments += 1;
        }
        return res;
    }
    fn header(&self) -> &bam::HeaderView {
        self.bam_reader.header()
    }
    fn finish(self) {
        complete_processes(
            self.processes,
            self.command_strings,
            self.log_file_descriptions,
            self.log_files,
            Some(self.tempdir))
    }
    fn num_detected_primary_alignments(&self) -> u64 {
        return self.num_detected_primary_alignments
    }
}

pub fn generate_named_bam_readers_from_bam_files(
    bam_paths: Vec<&str>) -> Vec<BamFileNamedReader>{

    bam_paths.iter().map(
        |path|

       BamFileNamedReader {
           stoit_name: std::path::Path::new(path).file_stem().unwrap().to_str().expect(
               "failure to convert bam file name to stoit name - UTF8 error maybe?").to_string(),
           bam_reader: bam::Reader::from_path(path).expect(
               &format!("Unable to find BAM file {}", path)),
           num_detected_primary_alignments: 0,
       }
    ).collect()
}

pub fn generate_named_bam_readers_from_reads(
    reference: &str,
    read1_path: &str,
    read2_path: Option<&str>,
    read_format: ReadFormat,
    threads: u16,
    cached_bam_file: Option<&str>,
    discard_unmapped: bool,
    bwa_options: Option<&str>,
    include_reference_in_stoit_name: bool) -> StreamingNamedBamReaderGenerator {

    let tmp_dir = TempDir::new("coverm_fifo")
        .expect("Unable to create temporary directory");
    let fifo_path = tmp_dir.path().join("foo.pipe");

    // create new fifo and give read, write and execute rights to the owner.
    // This is required because we cannot open a Rust stream as a BAM file with
    // rust-htslib.
    unistd::mkfifo(&fifo_path, stat::Mode::S_IRWXU)
        .expect(&format!("Error creating named pipe {:?}", fifo_path));

    let bwa_log = tempfile::NamedTempFile::new()
        .expect("Failed to create BWA log tempfile");
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
                match discard_unmapped { true => "-F4", false => "" },
                threads,
                path,
                samtools_view_cache_log.path().to_str()
                    .expect("Failed to convert tempfile path to str"))
        },
        None => format!("> {:?}", fifo_path)
    };
    let bwa_read_params = match read_format {
        ReadFormat::Interleaved => format!("-p '{}'", read1_path),
        ReadFormat::Coupled => format!("'{}' '{}'", read1_path, read2_path.unwrap()),
        ReadFormat::Single => format!("'{}'", read1_path),
    };
    let bwa_sort_prefix = tempfile::Builder::new()
        .prefix("coverm-make-samtools-sort")
        .tempfile_in(tmp_dir.path())
        .expect("Failed to create tempfile as samtools sort prefix");
    let cmd_string = format!(
        "set -e -o pipefail; \
         bwa mem {} -t {} '{}' {} 2>{} \
         | samtools sort -T '{}' -l0 -@ {} 2>{} \
         {}",
        // BWA
        bwa_options.unwrap_or(""), threads, reference, bwa_read_params,
        bwa_log.path().to_str().expect("Failed to convert tempfile path to str"),
        // samtools
        bwa_sort_prefix.path().to_str()
            .expect("Failed to convert bwa_sort_prefix tempfile to str"),
        threads-1,
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
        "BWA".to_string(),
        "samtools sort".to_string()];
    let mut log_files = vec![bwa_log, samtools2_log];
    if cached_bam_file.is_some() {
        log_descriptions.push("samtools view for cache".to_string());
        log_files.push(samtools_view_cache_log);
    }

    let stoit_name = match include_reference_in_stoit_name {
        true => std::path::Path::new(reference).file_name()
            .expect("Unable to convert reference to file name").to_str()
            .expect("Unable to covert file name into str").to_string()+"/",
        false => "".to_string()
    } + &std::path::Path::new(read1_path).file_name()
        .expect("Unable to convert read1 name to file name").to_str()
        .expect("Unable to covert file name into str").to_string();

    return StreamingNamedBamReaderGenerator {
        stoit_name: stoit_name,
        tempdir: tmp_dir,
        fifo_path: fifo_path,
        pre_processes: vec![cmd],
        command_strings: vec![format!("bash -c {}", cmd_string)],
        log_file_descriptions: log_descriptions,
        log_files: log_files,
    }
}


pub struct FilteredBamReader {
    stoit_name: String,
    filtered_stream: ReferenceSortedBamFilter,
}

impl NamedBamReader for FilteredBamReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }
    fn read(&mut self, mut record: &mut bam::record::Record) -> Result<(), bam::ReadError> {
        self.filtered_stream.read(&mut record)
    }
    fn header(&self) -> &bam::HeaderView {
        &self.filtered_stream.reader.header()
    }
    fn finish(self) {;}
    fn num_detected_primary_alignments(&self) -> u64 {
        return self.filtered_stream.num_detected_primary_alignments
    }
}

impl NamedBamReaderGenerator<FilteredBamReader> for FilteredBamReader {
    fn start(self) -> FilteredBamReader {
        FilteredBamReader {
            stoit_name: self.stoit_name,
            filtered_stream: self.filtered_stream,
        }
    }
}

pub fn generate_filtered_bam_readers_from_bam_files(
    bam_paths: Vec<&str>,
    flag_filters: FlagFilter,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32) -> Vec<FilteredBamReader>{

    let mut generators: Vec<FilteredBamReader> = vec![];

    for path in bam_paths {
        let filtered: FilteredBamReader;
        let stoit_name = std::path::Path::new(path).file_stem().unwrap().to_str().expect(
            "failure to convert bam file name to stoit name - UTF8 error maybe?").to_string();
        let reader = bam::Reader::from_path(path).expect(
            &format!("Unable to find BAM file {}", path));

        filtered = FilteredBamReader {
            stoit_name: stoit_name,
            filtered_stream: ReferenceSortedBamFilter::new(
                reader,
                flag_filters.clone(),
                min_aligned_length_single,
                min_percent_identity_single,
                min_aligned_percent_single,
                min_aligned_length_pair,
                min_percent_identity_pair,
                min_aligned_percent_pair,
                true),
        };

        generators.push(
            filtered
        )
    }

    return generators;
}














pub struct StreamingFilteredNamedBamReader {
    stoit_name: String,
    filtered_stream: ReferenceSortedBamFilter,
    tempdir: TempDir,
    processes: Vec<std::process::Child>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
}

pub struct StreamingFilteredNamedBamReaderGenerator {
    stoit_name: String,
    tempdir: TempDir,
    fifo_path: std::path::PathBuf,
    pre_processes: Vec<std::process::Command>,
    command_strings: Vec<String>,
    flag_filters: FlagFilter,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
}

impl NamedBamReaderGenerator<StreamingFilteredNamedBamReader> for StreamingFilteredNamedBamReaderGenerator {
    fn start(self) -> StreamingFilteredNamedBamReader {
        debug!("Starting mapping processes");
        let mut processes = vec![];
        for mut preprocess in self.pre_processes {
            processes.push(preprocess
                           .spawn()
                           .expect("Unable to execute bash"));
        }
        let bam_reader = bam::Reader::from_path(&self.fifo_path)
            .expect(&format!("Unable to find BAM file {:?}", self.fifo_path));
        let filtered_stream = ReferenceSortedBamFilter::new(
            bam_reader,
            self.flag_filters,
            self.min_aligned_length_single,
            self.min_percent_identity_single,
            self.min_aligned_percent_single,
            self.min_aligned_length_pair,
            self.min_percent_identity_pair,
            self.min_aligned_percent_pair,
            true);
        return StreamingFilteredNamedBamReader {
            stoit_name: self.stoit_name,
            filtered_stream: filtered_stream,
            tempdir: self.tempdir,
            processes: processes,
            command_strings: self.command_strings,
            log_file_descriptions: self.log_file_descriptions,
            log_files: self.log_files,
        }
    }
}

impl NamedBamReader for StreamingFilteredNamedBamReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }
    fn read(&mut self, record: &mut bam::record::Record) -> Result<(), bam::ReadError> {
        self.filtered_stream.read(record)
    }
    fn header(&self) -> &bam::HeaderView {
        self.filtered_stream.reader.header()
    }
    fn finish(self) {
        complete_processes(
            self.processes,
            self.command_strings,
            self.log_file_descriptions,
            self.log_files,
            Some(self.tempdir))
    }
    fn num_detected_primary_alignments(&self) -> u64 {
        return self.filtered_stream.num_detected_primary_alignments
    }
}


pub fn generate_filtered_named_bam_readers_from_reads(
    reference: &str,
    read1_path: &str,
    read2_path: Option<&str>,
    read_format: ReadFormat,
    threads: u16,
    cached_bam_file: Option<&str>,
    flag_filters: FlagFilter,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
    bwa_options: Option<&str>,
    discard_unmapped: bool,
    include_reference_in_stoit_name: bool)
    -> StreamingFilteredNamedBamReaderGenerator {

    let streaming = generate_named_bam_readers_from_reads(
        reference, read1_path, read2_path, read_format, threads,
        cached_bam_file, discard_unmapped, bwa_options,
        include_reference_in_stoit_name);
    return StreamingFilteredNamedBamReaderGenerator {
        stoit_name: streaming.stoit_name,
        tempdir: streaming.tempdir,
        fifo_path: streaming.fifo_path,
        pre_processes: streaming.pre_processes,
        command_strings: streaming.command_strings,
        log_file_descriptions: streaming.log_file_descriptions,
        log_files: streaming.log_files,
        flag_filters: flag_filters,
        min_aligned_length_single: min_aligned_length_single,
        min_percent_identity_single: min_percent_identity_single,
        min_aligned_percent_single: min_aligned_percent_single,
        min_aligned_length_pair: min_aligned_length_pair,
        min_percent_identity_pair: min_percent_identity_pair,
        min_aligned_percent_pair: min_aligned_percent_pair,
    }
}




pub struct BamGeneratorSet<T> {
    pub generators: Vec<T>,
    pub index: Box<BwaIndexStruct>,
}



pub struct NamedBamMaker {
    stoit_name: String,
    processes: Vec<std::process::Child>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
}

pub struct NamedBamMakerGenerator {
    stoit_name: String,
    pre_processes: Vec<std::process::Command>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
}

pub fn generate_bam_maker_generator_from_reads(
    reference: &str,
    read1_path: &str,
    read2_path: Option<&str>,
    read_format: ReadFormat,
    threads: u16,
    cached_bam_file: &str,
    discard_unmapped: bool,
    bwa_options: Option<&str>) -> NamedBamMakerGenerator {

    let bwa_log = tempfile::NamedTempFile::new()
        .expect("Failed to create BWA log tempfile");
    let samtools2_log = tempfile::NamedTempFile::new()
        .expect("Failed to create second samtools log tempfile");
    // tempfile does not need to be created but easier to create than get around
    // borrow checker.
    let samtools_view_cache_log = tempfile::NamedTempFile::new()
        .expect("Failed to create cache samtools view log tempfile");

    let bwa_read_params = match read_format {
        ReadFormat::Interleaved => format!("-p '{}'", read1_path),
        ReadFormat::Coupled => format!("'{}' '{}'", read1_path, read2_path.unwrap()),
        ReadFormat::Single => format!("'{}'", read1_path),
    };
    let bwa_sort_prefix = tempfile::Builder::new()
        .prefix("coverm-make-samtools-sort")
        .tempfile()
        .expect("Failed to create tempfile as samtools sort prefix");
    let cmd_string = format!(
        "set -e -o pipefail; \
         bwa mem {} -t {} '{}' {} 2>{} \
         | samtools sort -T '{}' -l0 -@ {} 2>{} \
         | samtools view {} -b -t {} -o '{}' 2>{}",
        // BWA
        bwa_options.unwrap_or(""), threads, reference, bwa_read_params,
        bwa_log.path().to_str().expect("Failed to convert tempfile path to str"),
        // samtools
        bwa_sort_prefix.path().to_str()
            .expect("Failed to convert bwa_sort_prefix tempfile to str"),
        threads-1,
        samtools2_log.path().to_str().expect("Failed to convert tempfile path to str"),
        // samtools view
        match discard_unmapped { true => "-F4", false => ""},
        threads,
        cached_bam_file,
        samtools_view_cache_log.path().to_str()
            .expect("Failed to convert tempfile path to str"));
    debug!("Queuing cmd_string: {}", cmd_string);
    let mut cmd = std::process::Command::new("bash");
    cmd
        .arg("-c")
        .arg(&cmd_string)
        .stderr(std::process::Stdio::piped());

    let log_descriptions = vec![
        "BWA".to_string(),
        "samtools sort".to_string(),
        "samtools view for cache".to_string()];
    let log_files = vec![
        bwa_log,
        samtools2_log,
        samtools_view_cache_log];

    return NamedBamMakerGenerator {
        stoit_name: std::path::Path::new(reference).file_name()
            .expect("Unable to convert reference to file name").to_str()
            .expect("Unable to covert file name into str").to_string()+"/"+
            &std::path::Path::new(read1_path).file_name()
            .expect("Unable to convert read1 name to file name").to_str()
            .expect("Unable to covert file name into str").to_string(),
        pre_processes: vec![cmd],
        command_strings: vec![format!("bash -c {}", cmd_string)],
        log_file_descriptions: log_descriptions,
        log_files: log_files,
    }
}

impl NamedBamReaderGenerator<NamedBamMaker> for NamedBamMakerGenerator {
    fn start(self) -> NamedBamMaker {
        debug!("Starting mapping processes");
        let mut processes = vec![];
        let mut i = 0;
        for mut preprocess in self.pre_processes {
            debug!("Running mapping command: {}", self.command_strings[i]);
            i += 1;
            processes.push(preprocess
                           .spawn()
                           .expect("Unable to execute bash"));
        }
        return NamedBamMaker {
            stoit_name: self.stoit_name,
            processes: processes,
            command_strings: self.command_strings,
            log_file_descriptions: self.log_file_descriptions,
            log_files: self.log_files,
        }
    }
}



impl NamedBamMaker {
    pub fn name(&self) -> &str {
        &(self.stoit_name)
    }
    pub fn finish(self) {
        complete_processes(
            self.processes,
            self.command_strings,
            self.log_file_descriptions,
            self.log_files,
            None)
    }
}

