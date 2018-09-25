use std;
use std::io::Read;

use filter::*;
use bwa_index_maintenance::BwaIndexStruct;

use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

use nix::unistd;
use nix::sys::stat;
use tempdir::TempDir;

pub trait NamedBamReader {
    // Name of the stoit
    fn name(&self) -> &str;

    // Read a record into record parameter
    fn read(&mut self, record: &mut bam::record::Record) -> Result<(), bam::ReadError>;

    // Return the bam header of the final BAM file
    fn header(&self) -> &bam::HeaderView;

    fn finish(&self);
}

pub trait NamedBamReaderGenerator<T> {
    // For readers that map, start the process of mapping
    fn start(&self) -> T;
}

pub struct BamFileNamedReader<'a> {
    stoit_name: &'a String,
    bam_reader: &'a mut bam::Reader
}

impl<'a> NamedBamReader for BamFileNamedReader<'a> {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }
    fn read(&mut self, record: &mut bam::record::Record) -> Result<(), bam::ReadError> {
        self.bam_reader.read(record)
    }
    fn header(&self) -> &bam::HeaderView {
        self.bam_reader.header()
    }
    fn finish(&self) {;}
}

impl<'a> NamedBamReaderGenerator<BamFileNamedReader<'a>> for BamFileNamedReader<'a> {
    fn start(&self) -> BamFileNamedReader<'a> {
        BamFileNamedReader<'a> {
            stoit_name: self.stoit_name,
            bam_reader: self.bam_reader
        }
    }
}

pub struct StreamingNamedBamReader<'a> {
    stoit_name: &'a String,
    bam_reader: bam::Reader,
    tempdir: TempDir,
    processes: &'a Vec<std::process::Child>,
}

pub struct StreamingNamedBamReaderGenerator<'a> {
    stoit_name: &'a String,
    tempdir: TempDir,
    fifo_path: std::path::PathBuf,
    pre_processes: &'a Vec<std::process::Command>,
}

impl<'a> NamedBamReaderGenerator<StreamingNamedBamReader<'a>> for StreamingNamedBamReaderGenerator<'a> {
    fn start(&self) -> StreamingNamedBamReader<'a> {
        debug!("Starting mapping processes");
        let mut processes = vec![];
        for preprocess in self.pre_processes {
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
            processes: &processes,
        }
    }
}

impl<'a> NamedBamReader for StreamingNamedBamReader<'a> {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }
    fn read(&mut self, record: &mut bam::record::Record) -> Result<(), bam::ReadError> {
        self.bam_reader.read(record)
    }
    fn header(&self) -> &bam::HeaderView {
        self.bam_reader.header()
    }
    fn finish(&self) {
        for mut process in self.processes {
            let es = process.wait().expect("Failed to glean exitstatus from mapping process");
            if !es.success() {
                error!("Error when running mapping process.");
                let mut err = String::new();
                process.stderr.expect("Failed to grab stderr from failed mapping process")
                    .read_to_string(&mut err).expect("Failed to read stderr into string");
                error!("The STDERR was: {:?}", err);
                panic!("Cannot continue since mapping failed.");
            }
        }
        // Close tempdir explicitly. Maybe not needed.
        self.tempdir.close().expect("Failed to close tempdir");
    }
}

pub fn generate_named_bam_readers_from_bam_files(
    bam_paths: Vec<&str>) -> Vec<BamFileNamedReader>{

    bam_paths.iter().map(
        |path|

       BamFileNamedReader {
           stoit_name: &std::path::Path::new(path).file_stem().unwrap().to_str().expect(
               "failure to convert bam file name to stoit name - UTF8 error maybe?").to_string(),
           bam_reader: &bam::Reader::from_path(path).expect(
               &format!("Unable to find BAM file {}", path))
       }
    ).collect()
}



pub fn generate_named_bam_readers_from_read_couple<'a>(
    reference: &str,
    read1_path: &str,
    read2_path: &str,
    threads: u16) -> StreamingNamedBamReaderGenerator<'a> {

    let tmp_dir = TempDir::new("coverm_fifo")
        .expect("Unable to create temporary directory");
    let fifo_path = tmp_dir.path().join("foo.pipe");

    // create new fifo and give read, write and execute rights to the owner.
    // This is required because we cannot open a Rust stream as a BAM file with
    // rust-htslib.
    unistd::mkfifo(&fifo_path, stat::Mode::S_IRWXU)
        .expect(&format!("Error creating named pipe {:?}", fifo_path));

    let cmd_string = format!(
        "set -e -o pipefail; \
         bwa mem -t {} '{}' '{}' '{}' \
         | samtools view -Sub -F4 \
         | samtools sort -l0 -@ {} -o {:?}",
        threads, reference, read1_path, read2_path,
        threads-1, fifo_path);
    debug!("Executing with bash: {}", cmd_string);
    let mut cmd = std::process::Command::new("bash");
    cmd
        .arg("-c")
        .arg(cmd_string)
        .stderr(std::process::Stdio::piped());

    return StreamingNamedBamReaderGenerator {
        stoit_name: &std::path::Path::new(read1_path).file_name()
            .expect("Unable to convert read1 name fq to path").to_str()
            .expect("Unable to covert file name into str").to_string(),
        tempdir: tmp_dir,
        fifo_path: fifo_path,
        pre_processes: &vec![cmd],
    }
}

pub struct FilteredBamReader<'a> {
    stoit_name: &'a String,
    filtered_stream: &'a ReferenceSortedBamFilter
}

impl<'a> NamedBamReader for FilteredBamReader<'a> {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }
    fn read(&mut self, mut record: &mut bam::record::Record) -> Result<(), bam::ReadError> {
        self.filtered_stream.read(&mut record)
    }
    fn header(&self) -> &bam::HeaderView {
        &self.filtered_stream.reader.header()
    }
    fn finish(&self) {;}
}

impl<'a> NamedBamReaderGenerator<FilteredBamReader<'a>> for FilteredBamReader<'a> {
    fn start(&self) -> FilteredBamReader<'a> {
        FilteredBamReader {
            stoit_name: self.stoit_name,
            filtered_stream: self.filtered_stream
        }
    }
}

pub fn generate_filtered_bam_readers_from_bam_files(
    bam_paths: Vec<&str>,
    min_aligned_length: u32,
    min_percent_identity: f32) -> Vec<FilteredBamReader>{

    let mut generators: Vec<FilteredBamReader> = vec![];

    for path in bam_paths {
        let filtered: FilteredBamReader;
        let stoit_name = std::path::Path::new(path).file_stem().unwrap().to_str().expect(
            "failure to convert bam file name to stoit name - UTF8 error maybe?").to_string();
        let reader = bam::Reader::from_path(path).expect(
            &format!("Unable to find BAM file {}", path));

        filtered = FilteredBamReader {
                stoit_name: &stoit_name,
                filtered_stream: &ReferenceSortedBamFilter::new(
                    reader,
                    min_aligned_length,
                    min_percent_identity)
            };

        generators.push(
            filtered
        )
    }

    return generators;
}














pub struct StreamingFilteredNamedBamReader<'a> {
    stoit_name: &'a String,
    filtered_stream: ReferenceSortedBamFilter,
    tempdir: TempDir,
    processes: &'a Vec<std::process::Child>
}

pub struct StreamingFilteredNamedBamReaderGenerator<'a> {
    stoit_name: &'a String,
    tempdir: TempDir,
    fifo_path: std::path::PathBuf,
    pre_processes: &'a Vec<std::process::Command>,
    min_aligned_length: u32,
    min_percent_identity: f32
}

impl<'a> NamedBamReaderGenerator<StreamingFilteredNamedBamReader<'a>> for StreamingFilteredNamedBamReaderGenerator<'a> {
    fn start(&self) -> StreamingFilteredNamedBamReader<'a> {
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
            self.min_aligned_length,
            self.min_percent_identity);
        return StreamingFilteredNamedBamReader {
            stoit_name: self.stoit_name,
            filtered_stream: filtered_stream,
            tempdir: self.tempdir,
            processes: &processes
        }
    }
}

impl<'a> NamedBamReader for StreamingFilteredNamedBamReader<'a> {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }
    fn read(&mut self, record: &mut bam::record::Record) -> Result<(), bam::ReadError> {
        self.filtered_stream.read(record)
    }
    fn header(&self) -> &bam::HeaderView {
        self.filtered_stream.reader.header()
    }
    fn finish(&self) {
        for mut process in self.processes {
            let es = process.wait().expect("Failed to glean exitstatus from mapping process");
            if !es.success() {
                error!("Error when running mapping process.");
                let mut err = String::new();
                process.stderr.expect("Failed to grab stderr from failed mapping process")
                    .read_to_string(&mut err).expect("Failed to read stderr into string");
                error!("The STDERR was: {:?}", err);
                panic!("Cannot continue since mapping failed.");
            }
        }
        // Close tempdir explicitly. Maybe not needed.
        self.tempdir.close().expect("Failed to close tempdir");
    }
}

pub fn generate_filtered_named_bam_readers_from_read_couple<'a>(
    reference: &str,
    read1_path: &str,
    read2_path: &str,
    threads: u16,
    min_aligned_length: u32,
    min_percent_identity: f32) -> StreamingFilteredNamedBamReaderGenerator<'a> {

    let tmp_dir = TempDir::new("coverm_fifo")
        .expect("Unable to create temporary directory");
    let fifo_path = tmp_dir.path().join("foo.pipe");

    // create new fifo and give read, write and execute rights to the owner.
    // This is required because we cannot open a Rust stream as a BAM file with
    // rust-htslib.
    unistd::mkfifo(&fifo_path, stat::Mode::S_IRWXU)
        .expect(&format!("Error creating named pipe {:?}", fifo_path));

    let cmd_string = format!(
        "set -e -o pipefail; \
         bwa mem -t {} '{}' '{}' '{}' \
         | samtools view -Sub -F4 \
         | samtools sort -l0 -@ {} -o {:?}",
        threads, reference, read1_path, read2_path,
        threads-1, fifo_path);
    debug!("Executing with bash: {}", cmd_string);
    let mut cmd = std::process::Command::new("bash");
    cmd
        .arg("-c")
        .arg(cmd_string)
        .stderr(std::process::Stdio::piped());

    return StreamingFilteredNamedBamReaderGenerator {
        stoit_name: &std::path::Path::new(read1_path).file_name()
            .expect("Unable to convert read1 name fq to path").to_str()
            .expect("Unable to covert file name into str").to_string(),
        tempdir: tmp_dir,
        fifo_path: fifo_path,
        pre_processes: &vec![cmd],
        min_aligned_length: min_aligned_length,
        min_percent_identity: min_percent_identity
    }
}


pub struct BamGeneratorSet<T> {
    pub generators: Vec<T>,
    pub index: Box<BwaIndexStruct>,
}
