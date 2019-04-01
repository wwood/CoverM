use std::fs::OpenOptions;
use std::io::prelude::*;
use std::io::BufWriter;
use std::path::Path;
use std::thread;

use fastq;
use fastq::Record;
use nix::unistd;
use nix::sys::stat;
use tempdir::TempDir;

#[allow(dead_code)]
struct OrderedReadBamGenerator {
}

#[allow(dead_code)]
impl OrderedReadBamGenerator {
    pub fn start_read_renumbering_writer(
        read_path: String
    ) -> std::thread::JoinHandle<()> { // TODO: Return the number of written reads.
        // Open tmpdir
        // Open pipe in tmpdir for writing
        let tmp_dir = TempDir::new("coverm_read_renaming_tmpdir")
            .expect("Unable to create temporary directory");
        let fifo_path = tmp_dir.path().join("read_pipe");
        // create new fifo and give read, write and execute rights to the owner.
        // This is required because we cannot open a Rust stream as a BAM file with
        // rust-htslib.
        unistd::mkfifo(&fifo_path, stat::Mode::S_IRWXU)
            .expect(&format!("Error creating named pipe 1 {:?}", fifo_path));
        // Open output streams to fifo
        //let f = OpenOptions::new().append(true).open(fifo_path);
        let f = OpenOptions::new().create(true).append(true).open("/tmp/blah")
            .expect("Unable to open renumbered read file");
        let mut buffer1 = BufWriter::new(f);

        // open and rename until exhausted in a new thread
        let child = thread::spawn(move || {
            let my_path = Some(Path::new(&read_path));
            let read_counter = fastq::parse_path(
                my_path,
                |parser| {
                    let mut read_counter: u64 = 0;
                    parser.each(|record| {
                        write!(buffer1,"@{}\n", read_counter).expect("Write error on numbered FASTQ writer");
                        buffer1.write(record.seq()).expect("Write error on numbered FASTQ writer");;
                        buffer1.write(b"\n").expect("Write error on numbered FASTQ writer");;
                        buffer1.write(b"+").expect("Write error on numbered FASTQ writer");;
                        buffer1.write(b"\n").expect("Write error on numbered FASTQ writer");;
                        buffer1.write(record.qual()).expect("Write error on numbered FASTQ writer");;
                        buffer1.write(b"\n").expect("Write error on numbered FASTQ writer");;
                        read_counter += 1;
                        true
                    }).expect("Invalid fastq file");

                    buffer1.flush().expect("Flush error on numbered FASTQ writer");
                }).expect("Invalid compression");
            read_counter
        });
        return child;
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_renumbering_hello_world() {
        OrderedReadBamGenerator::start_read_renumbering_writer("/tmp/1.fq".to_string()).join();
    }
}
