use std;
use rust_htslib::bam;
use rust_htslib::bam::Read;

use nix::unistd;
use nix::sys::stat;
use tempdir::TempDir;

pub trait NamedBamReader: Sized {
    fn name(&self) -> &str;

    fn read(&mut self, record: &mut bam::record::Record) -> Result<(), bam::ReadError>;

    fn header(&self) -> &bam::HeaderView;
}

pub struct BamFileNamedReader {
    stoit_name: String,
    bam_reader: bam::Reader
}

impl NamedBamReader for BamFileNamedReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }
    fn read(&mut self, record: &mut bam::record::Record) -> Result<(), bam::ReadError> {
        self.bam_reader.read(record)
    }
    fn header(&self) -> &bam::HeaderView {
        self.bam_reader.header()
    }
}

pub struct StreamingNamedBamReader {
    stoit_name: String,
    bam_reader: bam::Reader,
    pub processes: Vec<std::process::Child>
}

impl NamedBamReader for StreamingNamedBamReader {
    fn name(&self) -> &str {
        &(self.stoit_name)
    }
    fn read(&mut self, record: &mut bam::record::Record) -> Result<(), bam::ReadError> {
        self.bam_reader.read(record)
    }
    fn header(&self) -> &bam::HeaderView {
        self.bam_reader.header()
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
               &format!("Unable to find BAM file {}", path))
       }
    ).collect()
}



pub fn generate_named_bam_readers_from_read_couple(
    reference: &str,
    read1_path: &str,
    read2_path: &str) -> StreamingNamedBamReader {

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
         bwa mem '{}' '{}' '{}' \
         | samtools view -Sub -F4 \
         | samtools sort -o {:?}",
        reference, read1_path, read2_path,
        fifo_path);
    debug!("Executing with bash: {}", cmd_string);
    let cmd = std::process::Command::new("bash")
        .arg("-c")
        .arg(cmd_string)
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("Unable to execute bash");
    debug!("Spawned bash pipe");

    return StreamingNamedBamReader {
        stoit_name: std::path::Path::new(read1_path).file_name()
            .expect("Unable to convert read1 name fq to path").to_str()
            .expect("Unable to covert file name into str").to_string(),
        bam_reader: bam::Reader::from_path(&fifo_path)
            .expect(&format!("Unable to find BAM file {:?}", fifo_path)),
        processes: vec![cmd]
    }
}

