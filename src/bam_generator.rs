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
    bam_reader: bam::Reader// ,
    // processes: Vec<std::process::Child>
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

    let bwa_cmd = std::process::Command::new("bwa")
        .arg("mem")
        .arg("-t").arg("16")
        .arg(reference)
        .arg(read1_path)
        .arg(read2_path)
        .stdout(std::process::Stdio::piped())
        //.stderr(std::process::Stdio::piped()) // currently ignored
        .spawn()
        .expect("failed to execute bwa");
    debug!("Spawned BWA");
    // Convert to BAM
    let sam_to_bam = std::process::Command::new("samtools")
        .arg("view")
        .arg("-Sub")
        .stdin(bwa_cmd.stdout.expect("Unable to bind bwa STDOUT"))
        .stdout(std::process::Stdio::piped())
        //.stderr(std::process::Stdio::piped()) // currently ignored
        .spawn()
        .expect("failed to execute samtools view");
    debug!("Spawned sam to bam {:?}", sam_to_bam);
    // Reference-wise sort
    let bam_to_sorted = std::process::Command::new("samtools")
        .arg("sort")
        .stdin(sam_to_bam.stdout.expect("Unable to bind sam_to_bam stdout"))
        .arg("-o")
        .arg(&fifo_path)
        .spawn()
        .expect("Failed to execute samtools sort");
    debug!("Spawned samtools sort");

    return StreamingNamedBamReader {
        stoit_name: std::path::Path::new(read1_path).file_name()
            .expect("Unable to convert read1 name fq to path").to_str()
            .expect("Unable to covert file name into str").to_string(),
        bam_reader: bam::Reader::from_path(&fifo_path)
            .expect(&format!("Unable to find BAM file {:?}", fifo_path))// ,
        // processes: vec![bwa_cmd, sam_to_bam, bam_to_sorted]
    }
}

