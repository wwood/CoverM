use std::rc::Rc;
use std::str;

use rust_htslib::bam;
use std::collections::BTreeMap;
use rust_htslib::bam::Read;

pub struct ReferenceSortedBamFilter {
    first_set: BTreeMap<Rc<String>, Rc<bam::Record>>,
    current_reference: i32,
    known_next_read: Option<bam::Record>,
    pub reader: bam::Reader,
    min_aligned_length: u32,
    min_percent_identity: f32,
    pub num_detected_primary_alignments: u64,
}

impl ReferenceSortedBamFilter {
    pub fn new(
        reader: bam::Reader,
        min_aligned_length: u32,
        min_percent_identity: f32) -> ReferenceSortedBamFilter {

        ReferenceSortedBamFilter {
            first_set: BTreeMap::new(),
            current_reference: -1,
            known_next_read: None,
            reader: reader,
            min_aligned_length: min_aligned_length,
            min_percent_identity: min_percent_identity,
            num_detected_primary_alignments: 0,
        }
    }
}

impl ReferenceSortedBamFilter {
    pub fn read(&mut self, mut record: &mut bam::record::Record) -> Result<(), bam::ReadError> {
        if self.known_next_read.is_none() {
            while self.reader.read(&mut record).is_ok() {
                debug!("record: {:?}", record);

                debug!("passed flags, {} {} {}",
                       record.is_secondary(),
                       record.is_supplementary(),
                       !record.is_proper_pair());
                // TODO: make usage ensure flag_filtering when mapping
                if record.is_secondary() ||
                    record.is_supplementary() {
                        continue
                    }
                self.num_detected_primary_alignments += 1;
                if !record.is_proper_pair() {
                    continue
                }

                // if a new reference ID is encountered, instantiate a new first read set
                if record.tid() != self.current_reference {
                    self.current_reference = record.tid();
                    self.first_set = BTreeMap::new();
                }
                // if this is a first read
                if record.insert_size() > 0 {
                    if record.mtid() == self.current_reference {
                        // if tlen is +ve and < threshold
                        // add to first read set
                        let qname = String::from(str::from_utf8(record.qname())
                                                 .expect("UTF8 error in conversion of read name"));
                        debug!("Testing qname1 {}", qname);
                        self.first_set.insert(Rc::new(qname), Rc::new(record.clone()));
                        // continue the loop without returning as we need to see the second record
                    }
                    // pairs from different contigs are ignored.
                }
                else { // Second read in insert
                    let qname = String::from(str::from_utf8(record.qname())
                                             .expect("UTF8 error in conversion of read name"));

                    debug!("Testing qname2 {}", qname);
                    match self.first_set.remove(&qname) {
                        Some(record1) => {
                            if read_pair_passes_filter(
                                &record,
                                &record1,
                                self.min_aligned_length, self.min_percent_identity) {

                                debug!("Read pair passed QC");
                                self.known_next_read = Some(record.clone());
                                record.clone_from(
                                    &Rc::try_unwrap(record1).expect("Cannot get strong RC pointer"));
                                debug!("Returning..");
                                return Ok(())
                            } else {
                                debug!("Read pair did not pass QC");
                            }
                        },
                        // if pair is not in first read set, ignore it
                        None => {}
                    }
                }
            }

            // No more records, we are finished.
            return Err(bam::ReadError::NoMoreRecord)
        }


        else {
            record.clone_from(self.known_next_read.as_ref().unwrap());
            self.known_next_read = None;
            return Ok(())
        }
    }
}

fn read_pair_passes_filter(
    record1: &bam::Record,
    record2: &bam::Record,
    min_aligned_length: u32,
    min_percent_identity: f32) -> bool {

    let edit_distance1 = match record1.aux(b"NM") {
        Some(i) => i.integer(),
        None => {panic!("Alignment of read {:?} did not have an NM aux tag", record1.qname())}
    };
    let edit_distance2 = match record2.aux(b"NM") {
        Some(i) => i.integer(),
        None => {panic!("Alignment of read {:?} did not have an NM aux tag", record2.qname())}
    };

    let mut aligned_length1: u32 = 0;
    for cig in record1.cigar().iter() {
        match cig {
            bam::record::Cigar::Match(i) |
            bam::record::Cigar::Ins(i) => {
                aligned_length1 = aligned_length1 + i;
            },
            _ => {}
        }
    }
    let mut aligned_length2: u32 = 0;
    for cig in record2.cigar().iter() {
        match cig {
            bam::record::Cigar::Match(i) |
            bam::record::Cigar::Ins(i) => {
                aligned_length2 = aligned_length2 + i;
            },
            _ => {}
        }
    }

    let aligned = aligned_length1 + aligned_length2;
    debug!("num_bases {} {}, edit distances {} {}, perc {}",
             aligned_length1, aligned_length2, edit_distance1, edit_distance2,
             1.0 - ((edit_distance1 + edit_distance2) as f32 / aligned as f32));

    if aligned >= min_aligned_length &&
        1.0 - ((edit_distance1 + edit_distance2) as f32 / aligned as f32) >= min_percent_identity {
            return true
        } else {
            return false
        }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hello_world(){
        let mut reader = bam::Reader::from_path(
            &"tests/data/7seqs.reads_for_seq1_and_seq2.bam").unwrap();
        let mut sorted = ReferenceSortedBamFilter::new(
            reader, 90, 0.99);
        let queries = vec![
            "9",
            "9",
            "12",
            "12",
            "7",
            "7",
            "11",
            "11",
            "10",
            "10",
            "8",
            "8",
            "4",
            "4",
            "6",
            "6",
            "1",
            "1",
            "2",
            "2",
            "3",
            "3",
            "5",
            "5"];
        let mut record = bam::record::Record::new();
        for i in queries {
            println!("query: {}", i);
            sorted.read(&mut record).expect("");
            assert_eq!(i, str::from_utf8(record.qname()).unwrap());
        }
        assert!(sorted.read(&mut record).is_err())
    }

    #[test]
    fn test_one_bad_read(){
        let mut reader = bam::Reader::from_path(
            &"tests/data/2seqs.bad_read.1.bam").unwrap();
        let mut sorted = ReferenceSortedBamFilter::new(
            reader, 250, 0.99); // perc too high
        let queries = vec![
            "2",
            "2",
            "3",
            "3"];
        let mut record = bam::record::Record::new();
        for i in queries {
            println!("query: {}", i);
            sorted.read(&mut record).expect("");
            assert_eq!(i, str::from_utf8(record.qname()).unwrap());
        }

        let mut reader = bam::Reader::from_path(
            &"tests/data/2seqs.bad_read.1.bam").unwrap();
        let mut sorted = ReferenceSortedBamFilter::new(
            reader, 300, 0.98); // aligned length too high
        let queries = vec![
            "2",
            "2",
            "3",
            "3"];
        for i in queries {
            println!("query: {}", i);
            sorted.read(&mut record).expect("");
            assert_eq!(i, str::from_utf8(record.qname()).unwrap());
        }

        let mut reader = bam::Reader::from_path(
            &"tests/data/2seqs.bad_read.1.bam").unwrap();
        let mut sorted = ReferenceSortedBamFilter::new(
            reader, 299, 0.98);
        let queries = vec![
            "1",
            "1",
            "2",
            "2"];
        for i in queries {
            println!("query: {}", i);
            sorted.read(&mut record).expect("");
            assert_eq!(i, str::from_utf8(record.qname()).unwrap());
        }
    }
}
