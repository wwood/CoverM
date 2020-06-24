use std::str;

#[macro_use]
extern crate log;
extern crate rust_htslib;
use rust_htslib::bam;
use rust_htslib::bam::Read;

extern crate coverm;
use coverm::genomes_and_contigs::find_first;

fn main() {
    // Open SAM format input
    let mut sam_input = bam::Reader::from_stdin().expect("Failed to read SAM input");
    let header = sam_input.header().clone();
    let target_names = header.target_names();

    // Open BAM output format
    let mut bam_output = bam::Writer::from_stdout(
        &bam::header::Header::from_template(&header),
        bam::Format::BAM,
    )
    .expect("Failed to start BAM output");
    bam_output
        .set_compression_level(bam::CompressionLevel::Uncompressed)
        .expect("Failed to set bam compression level");

    // for each record, is it the same as the last read name?
    let mut current_record_name: Vec<u8> = vec![];
    let mut num_mapped_reads: u64 = 0;
    let mut num_unmapped_reads: u64 = 0;

    let mut process_current_read = || {
        num_mapped_reads += 1;

        //FIXME implement
        // For each vector entry in
        // if Some, print as is except set as primary alignment
    };

    let num_reference_sets: usize = 4; //FIXME read in from ARGV
                                       // TODO: make into Vec<bool> and add another vec to hold data, to avoid realloc
    let mut current_best_hits: Vec<Option<bam::Record>> = vec![None; num_reference_sets];

    for record_res in sam_input.records() {
        let record = record_res.expect("Failed to parse SAM record");
        // If new read process the last read
        let qname = record.qname().to_vec();
        if current_record_name != qname {
            process_current_read();
            current_best_hits = vec![None; num_reference_sets];
            current_record_name = qname;
        }

        let tid = record.tid();
        if tid < 0 {
            // Unmapped read, ignore
            num_unmapped_reads += 1;
        } else {
            let reference_id = extract_reference(tid as u32, &target_names);

            //TODO: Check that there is no new best alignment
            let replace = current_best_hits[reference_id].is_none();
            if replace {
                current_best_hits[reference_id] = Some(record.clone());
            }
        }
    }
    // if current_record_name != vec!() {
    //     process_current_read();
    // }

    info!(
        "Found {} reads mapped, with {} mapped and {} unmapped",
        num_mapped_reads + num_unmapped_reads,
        num_mapped_reads,
        num_unmapped_reads
    );

    // Find the best hit for each assembly, and mark that as primary.
}

fn extract_reference<'a>(tid: u32, target_names: &'a Vec<&[u8]>) -> usize {
    let split_char = b'~';
    let target_name = target_names[tid as usize];
    trace!("target name {:?}, separator {:?}", target_name, split_char);
    let offset = find_first(target_name, split_char).expect(
        &format!("Contig name {} does not contain split symbol, so cannot determine which genome it belongs to",
                 str::from_utf8(target_name).unwrap()));
    str::from_utf8(&target_name[(0..offset)])
        .expect(&format!(
            "Failed to discern reference sequence set number from {:?}",
            target_name
        ))
        .parse()
        .expect(&format!(
            "Failed to parse reference sequence set number from {:?}",
            target_name
        ))
}
