use std;
use std::io::Read;
use std::process;
use std::sync::atomic::{compiler_fence, Ordering};

use filter::*;
use mapping_index_maintenance::MappingIndex;
use mapping_parameters::ReadFormat;
use FlagFilter;

use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

use nix::sys::stat;
use nix::unistd;
use tempdir::TempDir;
use tempfile;

use rust_htslib::errors::Result as HtslibResult;

pub trait NamedBamReader {
    // Name of the stoit
    fn name(&self) -> &str;

    // Read a record into record parameter
    fn read(&mut self, record: &mut bam::record::Record) -> Option<HtslibResult<()>>;

    // Return the bam header of the final BAM file
    fn header(&self) -> &bam::HeaderView;

    fn finish(self);

    //set the number of threads for Bam reading
    fn set_threads(&mut self, n_threads: usize);

    // Number of reads that were detected
    fn num_detected_primary_alignments(&self) -> u64;
}

pub trait NamedBamReaderGenerator<T> {
    // For readers that map, start the process of mapping
    fn start(self) -> T;
}

#[derive(Debug, Clone, Copy)]
#[allow(non_camel_case_types)]
pub enum MappingProgram {
    BWA_MEM,
    BWA_MEM2,
    MINIMAP2_SR,
    MINIMAP2_ONT,
    MINIMAP2_PB,
    MINIMAP2_HIFI,
    MINIMAP2_NO_PRESET,
    STROBEALIGN,
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
    fn read(&mut self, record: &mut bam::record::Record) -> Option<HtslibResult<()>> {
        let res = self.bam_reader.read(record);
        if res == Some(Ok(())) && !record.is_secondary() && !record.is_supplementary() {
            self.num_detected_primary_alignments += 1;
        }
        res
    }
    fn header(&self) -> &bam::HeaderView {
        self.bam_reader.header()
    }
    fn finish(self) {}

    fn set_threads(&mut self, n_threads: usize) {
        if n_threads > 1 {
            self.bam_reader.set_threads(n_threads - 1).unwrap();
        }
    }

    fn num_detected_primary_alignments(&self) -> u64 {
        self.num_detected_primary_alignments
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
    minimap2_log_file_index: Option<usize>,
}

pub struct StreamingNamedBamReaderGenerator {
    stoit_name: String,
    tempdir: TempDir,
    fifo_path: std::path::PathBuf,
    pre_processes: Vec<std::process::Command>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
    minimap2_log_file_index: Option<usize>,
}

impl NamedBamReaderGenerator<StreamingNamedBamReader> for StreamingNamedBamReaderGenerator {
    fn start(self) -> StreamingNamedBamReader {
        debug!("Starting mapping processes");
        let mut processes = vec![];
        for (i, mut preprocess) in self.pre_processes.into_iter().enumerate() {
            debug!("Running mapping command: {}", self.command_strings[i]);
            processes.push(preprocess.spawn().expect("Unable to execute bash"));
        }
        let bam_reader = match bam::Reader::from_path(&self.fifo_path) {
            Ok(reader) => reader,
            Err(upstream_error) => {
                error!(
                    "Failed to correctly find or parse BAM file at {:?}: {}",
                    self.fifo_path, upstream_error
                );
                complete_processes(
                    processes,
                    self.command_strings,
                    self.log_file_descriptions,
                    self.log_files,
                    Some(self.tempdir),
                );
                panic!("Failure to find or parse BAM file, cannot continue");
            }
        };
        StreamingNamedBamReader {
            stoit_name: self.stoit_name,
            bam_reader,
            tempdir: self.tempdir,
            processes,
            command_strings: self.command_strings,
            log_file_descriptions: self.log_file_descriptions,
            log_files: self.log_files,
            num_detected_primary_alignments: 0,
            minimap2_log_file_index: self.minimap2_log_file_index,
        }
    }
}

pub fn name_stoit(
    index_path: &str,
    read1_path: &str,
    include_reference_in_stoit_name: bool,
) -> String {
    (match include_reference_in_stoit_name {
        true => (std::path::Path::new(&index_path)
            .file_name()
            .expect("Unable to convert reference to file name")
            .to_str()
            .expect("Unable to covert file name into str")
            .to_string()
            + "/")
            .to_string(),
        false => "".to_string(),
    }) + std::path::Path::new(read1_path)
        .file_name()
        .expect("Unable to convert read1 name to file name")
        .to_str()
        .expect("Unable to covert file name into str")
}

pub fn complete_processes(
    processes: Vec<std::process::Child>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
    tempdir: Option<TempDir>,
) {
    let mut failed_any = false;
    let mut overall_stderrs = vec![];
    for mut process in processes {
        let es = process
            .wait()
            .expect("Failed to glean exitstatus from mapping process");
        let failed = !es.success();
        if failed || log_enabled!(log::Level::Debug) {
            if failed {
                failed_any = true;
                error!("Error when running mapping process. Exitstatus was {:?}. Command run was: {:?}", es, command_strings);
            } else {
                debug!("Successfully finished process {:?}", process);
            }
            let mut err = String::new();
            process
                .stderr
                .expect("Failed to grab stderr from failed mapping process")
                .read_to_string(&mut err)
                .expect("Failed to read stderr into string");
            debug!("The overall STDERR was: {:?}", err);
            overall_stderrs.push(err);
        }
    }
    if failed_any || log_enabled!(log::Level::Debug) {
        for (description, tf) in log_file_descriptions.iter().zip(log_files.into_iter()) {
            let mut contents = String::new();
            tf.into_file()
                .read_to_string(&mut contents)
                .unwrap_or_else(|_| panic!("Failed to read log file for {}", description));
            if failed_any {
                error!("The STDERR for the {:} part was: {}", description, contents);
            } else {
                debug!("The STDERR for the {:} part was: {}", description, contents);
            }
        }
    }
    if failed_any {
        error!("Cannot continue since mapping failed.");
        process::exit(1);
    }

    // There's (maybe) a (difficult to reproduce) single-thread issue where the
    // tempdir gets dropped before the process is finished. Hopefully putting a
    // compiler fence here stops this.
    compiler_fence(Ordering::SeqCst);
    debug!("After fence, for tempdir {:?}", tempdir);
}

impl NamedBamReader for StreamingNamedBamReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }
    fn read(&mut self, record: &mut bam::record::Record) -> Option<HtslibResult<()>> {
        let res = self.bam_reader.read(record);
        if res == Some(Ok(())) && !record.is_secondary() && !record.is_supplementary() {
            self.num_detected_primary_alignments += 1;
        }
        res
    }
    fn header(&self) -> &bam::HeaderView {
        self.bam_reader.header()
    }
    fn finish(self) {
        // Check minimap2 didn't complain about unequal numbers of reads
        match self.minimap2_log_file_index {
            None => {}
            Some(log_file_index) => {
                let mut contents = String::new();
                std::fs::File::open(&self.log_files[log_file_index])
                    .expect("Failed to read minimap2 log file")
                    .read_to_string(&mut contents)
                    .expect("Failed to read minimap2 log file to string");
                if contents.contains("query files have different number of records") {
                    error!("The STDERR for the minimap2 part was: {}", contents);
                    error!(
                        "Not continuing since when input file pairs have \
                        unequal numbers of reads this usually means \
                        incorrect / corrupt files were specified"
                    );
                    process::exit(1);
                }
            }
        };

        complete_processes(
            self.processes,
            self.command_strings,
            self.log_file_descriptions,
            self.log_files,
            Some(self.tempdir),
        );
    }

    fn set_threads(&mut self, n_threads: usize) {
        if n_threads > 1 {
            self.bam_reader.set_threads(n_threads - 1).unwrap();
        }
    }

    fn num_detected_primary_alignments(&self) -> u64 {
        self.num_detected_primary_alignments
    }
}

pub fn generate_named_bam_readers_from_bam_files(bam_paths: Vec<&str>) -> Vec<BamFileNamedReader> {
    bam_paths
        .iter()
        .map(|path| BamFileNamedReader {
            stoit_name: std::path::Path::new(path)
                .file_stem()
                .unwrap()
                .to_str()
                .expect("failure to convert bam file name to stoit name - UTF8 error maybe?")
                .to_string(),
            bam_reader: bam::Reader::from_path(path)
                .unwrap_or_else(|_| panic!("Unable to find BAM file {}", path)),
            num_detected_primary_alignments: 0,
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
pub fn generate_named_bam_readers_from_reads(
    mapping_program: MappingProgram,
    index: &dyn MappingIndex,
    read1_path: &str,
    read2_path: Option<&str>,
    read_format: ReadFormat,
    threads: u16,
    cached_bam_file: Option<&str>,
    discard_unmapped: bool,
    mapping_options: Option<&str>,
    include_reference_in_stoit_name: bool,
) -> StreamingNamedBamReaderGenerator {
    let tmp_dir = TempDir::new("coverm_fifo").expect("Unable to create temporary directory");
    let fifo_path = tmp_dir.path().join("foo.pipe");

    // create new fifo and give read, write and execute rights to the owner.
    // This is required because we cannot open a Rust stream as a BAM file with
    // rust-htslib.
    unistd::mkfifo(&fifo_path, stat::Mode::S_IRWXU)
        .unwrap_or_else(|_| panic!("Error creating named pipe {:?}", fifo_path));

    let mapping_log = tempfile::Builder::new()
        .prefix("coverm-mapping-log")
        .tempfile()
        .unwrap_or_else(|_| panic!("Failed to create {:?} log tempfile", mapping_program));
    let samtools2_log = tempfile::Builder::new()
        .prefix("coverm-samtools2-log")
        .tempfile()
        .expect("Failed to create second samtools log tempfile");
    // tempfile does not need to be created but easier to create than get around
    // borrow checker.
    let samtools_view_cache_log = tempfile::Builder::new()
        .prefix("coverm-samtools-view-log")
        .tempfile()
        .expect("Failed to create cache samtools view log tempfile");

    let cached_bam_file_args = match cached_bam_file {
        Some(path) => {
            format!(
                "|tee {:?} |samtools view {} -@ {} -b -o '{}' 2>{}",
                // tee
                fifo_path,
                // samtools view
                match discard_unmapped {
                    true => "-F4",
                    false => "",
                },
                threads - 1,
                path,
                samtools_view_cache_log
                    .path()
                    .to_str()
                    .expect("Failed to convert tempfile path to str")
            )
        }
        None => format!("> {:?}", fifo_path),
    };

    let mapping_command = build_mapping_command(
        mapping_program,
        read_format,
        threads,
        read1_path,
        index,
        read2_path,
        mapping_options,
    );
    let bwa_sort_prefix = tempfile::Builder::new()
        .prefix("coverm-make-samtools-sort")
        .tempfile_in(tmp_dir.path())
        .expect("Failed to create tempfile as samtools sort prefix");
    let cmd_string = format!(
        "set -e -o pipefail; \
         {} 2>{} \
         | samtools sort -T '{}' -l0 -@ {} 2>{} \
         {}",
        // Mapping program
        mapping_command,
        mapping_log
            .path()
            .to_str()
            .expect("Failed to convert tempfile path to str"),
        // samtools
        bwa_sort_prefix
            .path()
            .to_str()
            .expect("Failed to convert bwa_sort_prefix tempfile to str"),
        threads - 1,
        samtools2_log
            .path()
            .to_str()
            .expect("Failed to convert tempfile path to str"),
        // Caching (or not)
        cached_bam_file_args
    );
    debug!("Queuing cmd_string: {}", cmd_string);
    let mut cmd = std::process::Command::new("bash");
    cmd.arg("-c")
        .arg(&cmd_string)
        .stderr(std::process::Stdio::piped());

    // Required because of https://github.com/wwood/CoverM/issues/58
    let minimap2_log_file_index = match mapping_program {
        MappingProgram::BWA_MEM | MappingProgram::BWA_MEM2 | MappingProgram::STROBEALIGN => None,
        // Required because of https://github.com/lh3/minimap2/issues/527
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_PB
        | MappingProgram::MINIMAP2_HIFI
        | MappingProgram::MINIMAP2_NO_PRESET => Some(0),
    };

    let mut log_descriptions = vec![
        format!("{:?}", mapping_program),
        "samtools sort".to_string(),
    ];
    let mut log_files = vec![mapping_log, samtools2_log];
    if cached_bam_file.is_some() {
        log_descriptions.push("samtools view for cache".to_string());
        log_files.push(samtools_view_cache_log);
    }

    let stoit_name = name_stoit(
        index.index_path(),
        read1_path,
        include_reference_in_stoit_name,
    );

    StreamingNamedBamReaderGenerator {
        stoit_name,
        tempdir: tmp_dir,
        fifo_path,
        pre_processes: vec![cmd],
        command_strings: vec![format!("bash -c \"{}\"", cmd_string)],
        log_file_descriptions: log_descriptions,
        log_files,
        minimap2_log_file_index,
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
    fn read(&mut self, record: &mut bam::record::Record) -> Option<HtslibResult<()>> {
        self.filtered_stream.read(record)
    }
    fn header(&self) -> &bam::HeaderView {
        self.filtered_stream.reader.header()
    }
    fn finish(self) {}
    fn set_threads(&mut self, n_threads: usize) {
        if n_threads > 1 {
            self.filtered_stream
                .reader
                .set_threads(n_threads - 1)
                .unwrap();
        }
    }
    fn num_detected_primary_alignments(&self) -> u64 {
        self.filtered_stream.num_detected_primary_alignments
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

#[allow(clippy::too_many_arguments)]
pub fn generate_filtered_bam_readers_from_bam_files(
    bam_paths: Vec<&str>,
    flag_filters: FlagFilter,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    min_mapq: u8,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
) -> Vec<FilteredBamReader> {
    let mut generators: Vec<FilteredBamReader> = vec![];

    for path in bam_paths {
        let filtered: FilteredBamReader;
        let stoit_name = std::path::Path::new(path)
            .file_stem()
            .unwrap()
            .to_str()
            .expect("failure to convert bam file name to stoit name - UTF8 error maybe?")
            .to_string();
        let reader = bam::Reader::from_path(path)
            .unwrap_or_else(|_| panic!("Unable to find BAM file {}", path));

        filtered = FilteredBamReader {
            stoit_name,
            filtered_stream: ReferenceSortedBamFilter::new(
                reader,
                flag_filters.clone(),
                min_aligned_length_single,
                min_percent_identity_single,
                min_aligned_percent_single,
                min_mapq,
                min_aligned_length_pair,
                min_percent_identity_pair,
                min_aligned_percent_pair,
                true,
            ),
        };

        generators.push(filtered)
    }

    generators
}

pub struct StreamingFilteredNamedBamReader {
    stoit_name: String,
    filtered_stream: ReferenceSortedBamFilter,
    tempdir: TempDir,
    processes: Vec<std::process::Child>,
    command_strings: Vec<String>,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
    output_filtered_bam_stream: Option<bam::Writer>,
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
    min_mapq: u8,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
    log_file_descriptions: Vec<String>,
    log_files: Vec<tempfile::NamedTempFile>,
    output_filtered_bam_path: Option<String>,
}

impl NamedBamReaderGenerator<StreamingFilteredNamedBamReader>
    for StreamingFilteredNamedBamReaderGenerator
{
    fn start(self) -> StreamingFilteredNamedBamReader {
        debug!("Starting mapping processes");
        let mut processes = vec![];
        for mut preprocess in self.pre_processes {
            processes.push(preprocess.spawn().expect("Unable to execute bash"));
        }
        let bam_reader = match bam::Reader::from_path(&self.fifo_path) {
            Ok(reader) => reader,
            Err(upstream_error) => {
                error!(
                    "Failed to correctly find or parse BAM file at {:?}: {}",
                    self.fifo_path, upstream_error
                );
                complete_processes(
                    processes,
                    self.command_strings,
                    self.log_file_descriptions,
                    self.log_files,
                    Some(self.tempdir),
                );
                panic!("Failure to find or parse BAM file, cannot continue");
            }
        };

        // If there's a output_filtered_bam_file, open it for writing
        let output_filtered_bam_stream = match self.output_filtered_bam_path {
            Some(path) => {
                let output_filtered_bam_stream = bam::Writer::from_path(
                    path,
                    &bam::header::Header::from_template(bam_reader.header()),
                    bam::Format::Bam,
                )
                .expect("Failed to open output filtered BAM file for writing");
                // TODO: Set compression level, num threads?
                Some(output_filtered_bam_stream)
            }
            None => None,
        };

        let filtered_stream = ReferenceSortedBamFilter::new(
            bam_reader,
            self.flag_filters,
            self.min_aligned_length_single,
            self.min_percent_identity_single,
            self.min_aligned_percent_single,
            self.min_mapq,
            self.min_aligned_length_pair,
            self.min_percent_identity_pair,
            self.min_aligned_percent_pair,
            true,
        );
        StreamingFilteredNamedBamReader {
            stoit_name: self.stoit_name,
            filtered_stream,
            tempdir: self.tempdir,
            processes,
            command_strings: self.command_strings,
            log_file_descriptions: self.log_file_descriptions,
            log_files: self.log_files,
            output_filtered_bam_stream,
        }
    }
}

impl NamedBamReader for StreamingFilteredNamedBamReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }
    fn read(&mut self, record: &mut bam::record::Record) -> Option<HtslibResult<()>> {
        // Write the record to the output filtered BAM stream if it exists
        if let Some(ref mut output_stream) = self.output_filtered_bam_stream {
            if !record.is_secondary() && !record.is_supplementary() {
                output_stream
                    .write(record)
                    .expect("Failed to write record to output filtered BAM stream");
            }
        }
        self.filtered_stream.read(record)
    }
    fn header(&self) -> &bam::HeaderView {
        self.filtered_stream.reader.header()
    }
    fn finish(self) {
        debug!(
            "Finishing StreamingFilteredNamedBamReader. Tempdir is {:?}",
            self.tempdir.path()
        );
        complete_processes(
            self.processes,
            self.command_strings,
            self.log_file_descriptions,
            self.log_files,
            Some(self.tempdir),
        )
    }

    fn set_threads(&mut self, n_threads: usize) {
        if n_threads > 1 {
            self.filtered_stream
                .reader
                .set_threads(n_threads - 1)
                .unwrap();
        }
    }
    fn num_detected_primary_alignments(&self) -> u64 {
        self.filtered_stream.num_detected_primary_alignments
    }
}

#[allow(clippy::too_many_arguments)]
pub fn generate_filtered_named_bam_readers_from_reads(
    mapping_program: MappingProgram,
    index: &dyn MappingIndex,
    read1_path: &str,
    read2_path: Option<&str>,
    read_format: ReadFormat,
    threads: u16,
    cached_bam_file: Option<&str>,
    output_filtered_bam_files: Vec<&str>,
    flag_filters: FlagFilter,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    min_mapq: u8,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
    bwa_options: Option<&str>,
    discard_unmapped: bool,
    include_reference_in_stoit_name: bool,
) -> StreamingFilteredNamedBamReaderGenerator {
    let streaming = generate_named_bam_readers_from_reads(
        mapping_program,
        index,
        read1_path,
        read2_path,
        read_format,
        threads,
        cached_bam_file,
        discard_unmapped,
        bwa_options,
        include_reference_in_stoit_name,
    );
    StreamingFilteredNamedBamReaderGenerator {
        stoit_name: streaming.stoit_name,
        tempdir: streaming.tempdir,
        fifo_path: streaming.fifo_path,
        pre_processes: streaming.pre_processes,
        command_strings: streaming.command_strings,
        log_file_descriptions: streaming.log_file_descriptions,
        log_files: streaming.log_files,
        flag_filters,
        min_aligned_length_single,
        min_percent_identity_single,
        min_aligned_percent_single,
        min_mapq,
        min_aligned_length_pair,
        min_percent_identity_pair,
        min_aligned_percent_pair,
        output_filtered_bam_file,
    }
}

pub struct BamGeneratorSet<T> {
    pub generators: Vec<T>,
    pub index: Box<dyn MappingIndex>,
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

#[allow(clippy::too_many_arguments)]
pub fn generate_bam_maker_generator_from_reads(
    mapping_program: MappingProgram,
    index: &dyn MappingIndex,
    read1_path: &str,
    read2_path: Option<&str>,
    read_format: ReadFormat,
    threads: u16,
    cached_bam_file: &str,
    discard_unmapped: bool,
    mapping_options: Option<&str>,
) -> NamedBamMakerGenerator {
    let mapping_log = tempfile::Builder::new()
        .prefix("coverm-mapping-log")
        .tempfile()
        .unwrap_or_else(|_| panic!("Failed to create {:?} log tempfile", mapping_program));
    let samtools2_log = tempfile::Builder::new()
        .prefix("coverm-samtools2-log")
        .tempfile()
        .expect("Failed to create second samtools log tempfile");
    // tempfile does not need to be created but easier to create than get around
    // borrow checker.
    let samtools_view_cache_log = tempfile::Builder::new()
        .prefix("coverm-samtools-view-log")
        .tempfile()
        .expect("Failed to create cache samtools view log tempfile");

    let mapping_command = build_mapping_command(
        mapping_program,
        read_format,
        threads,
        read1_path,
        index,
        read2_path,
        mapping_options,
    );

    let bwa_sort_prefix = tempfile::Builder::new()
        .prefix("coverm-make-samtools-sort")
        .tempfile()
        .expect("Failed to create tempfile as samtools sort prefix");
    let cmd_string = format!(
        "set -e -o pipefail; \
         {} 2>{} \
         | samtools sort -T '{}' -l0 -@ {} 2>{} \
         | samtools view {} -b -@ {} -o '{}' 2>{}",
        // Mapping program
        mapping_command,
        mapping_log
            .path()
            .to_str()
            .expect("Failed to convert tempfile path to str"),
        // samtools
        bwa_sort_prefix
            .path()
            .to_str()
            .expect("Failed to convert bwa_sort_prefix tempfile to str"),
        threads - 1,
        samtools2_log
            .path()
            .to_str()
            .expect("Failed to convert tempfile path to str"),
        // samtools view
        match discard_unmapped {
            true => "-F4",
            false => "",
        },
        threads - 1,
        cached_bam_file,
        samtools_view_cache_log
            .path()
            .to_str()
            .expect("Failed to convert tempfile path to str")
    );
    debug!("Queuing cmd_string: {}", cmd_string);
    let mut cmd = std::process::Command::new("bash");
    cmd.arg("-c")
        .arg(&cmd_string)
        .stderr(std::process::Stdio::piped());

    let log_descriptions = vec![
        format!("{:?}", mapping_program),
        "samtools sort".to_string(),
        "samtools view for cache".to_string(),
    ];
    let log_files = vec![mapping_log, samtools2_log, samtools_view_cache_log];

    return NamedBamMakerGenerator {
        stoit_name: name_stoit(index.index_path(), read1_path, true),
        pre_processes: vec![cmd],
        command_strings: vec![format!("bash -c \"{}\"", cmd_string)],
        log_file_descriptions: log_descriptions,
        log_files,
    };
}

impl NamedBamReaderGenerator<NamedBamMaker> for NamedBamMakerGenerator {
    fn start(self) -> NamedBamMaker {
        debug!("Starting mapping processes");
        let mut processes = vec![];
        for (i, mut preprocess) in self.pre_processes.into_iter().enumerate() {
            debug!("Running mapping command: {}", self.command_strings[i]);
            processes.push(preprocess.spawn().expect("Unable to execute bash"));
        }
        NamedBamMaker {
            stoit_name: self.stoit_name,
            processes,
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
            None,
        )
    }
}

pub fn build_mapping_command(
    mapping_program: MappingProgram,
    read_format: ReadFormat,
    threads: u16,
    read1_path: &str,
    reference: &dyn MappingIndex,
    read2_path: Option<&str>,
    mapping_options: Option<&str>,
) -> String {
    let read_params1 = match mapping_program {
        // minimap2 auto-detects interleaved based on read names
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_HIFI
        | MappingProgram::MINIMAP2_PB
        | MappingProgram::MINIMAP2_NO_PRESET => "",
        MappingProgram::BWA_MEM | MappingProgram::BWA_MEM2 => match read_format {
            ReadFormat::Interleaved => "-p",
            ReadFormat::Coupled | ReadFormat::Single => "",
        },
        MappingProgram::STROBEALIGN => match read_format {
            ReadFormat::Interleaved => "--interleaved",
            ReadFormat::Coupled | ReadFormat::Single => "",
        },
    };

    let read_params2 = match read_format {
        ReadFormat::Interleaved => format!("'{}'", read1_path),
        ReadFormat::Coupled => format!("'{}' '{}'", read1_path, read2_path.unwrap()),
        ReadFormat::Single => format!("'{}'", read1_path),
    };

    format!(
        "{} {} -t {} {} {} '{}' {}",
        match mapping_program {
            MappingProgram::BWA_MEM => "bwa mem".to_string(),
            MappingProgram::BWA_MEM2 => "bwa-mem2 mem".to_string(),
            MappingProgram::STROBEALIGN => "strobealign".to_string(),
            _ => {
                let split_prefix = tempfile::Builder::new()
                    .prefix("coverm-minimap2-split-index")
                    .tempfile()
                    .unwrap_or_else(|_| {
                        panic!(
                            "Failed to create {:?} minimap2 split_prefix file",
                            mapping_program
                        )
                    });
                format!(
                    "minimap2 --split-prefix {} -a {}",
                    split_prefix
                        .path()
                        .to_str()
                        .expect("Failed to convert split prefix tempfile path to str"),
                    match mapping_program {
                        MappingProgram::BWA_MEM
                        | MappingProgram::BWA_MEM2
                        | MappingProgram::STROBEALIGN => unreachable!(),
                        MappingProgram::MINIMAP2_SR => "-x sr",
                        MappingProgram::MINIMAP2_ONT => "-x map-ont",
                        MappingProgram::MINIMAP2_HIFI => "-x map-hifi",
                        MappingProgram::MINIMAP2_PB => "-x map-pb",
                        MappingProgram::MINIMAP2_NO_PRESET => "",
                    }
                )
            }
        },
        mapping_options.unwrap_or(""),
        threads,
        read_params1,
        reference.command_prefix(),
        reference.index_path(),
        read_params2
    )
}
