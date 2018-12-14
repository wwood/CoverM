use std::rc::Rc;
use std::str;
use std::collections::BTreeMap;

use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Cigar;

pub struct ReferenceSortedBamFilter {
    first_set: BTreeMap<Rc<String>, Rc<bam::Record>>,
    current_reference: i32,
    known_next_read: Option<bam::Record>,
    pub reader: bam::Reader,
    filter_single_reads: bool,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32,
    filter_pairs: bool,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32,
    pub num_detected_primary_alignments: u64,
}

impl ReferenceSortedBamFilter {
    pub fn new(
        reader: bam::Reader,
        min_aligned_length_single: u32,
        min_percent_identity_single: f32,
        min_aligned_percent_single: f32,
        min_aligned_length_pair: u32,
        min_percent_identity_pair: f32,
        min_aligned_percent_pair: f32) -> ReferenceSortedBamFilter {

        let filtering_single =
            min_aligned_length_single > 0 ||
            min_percent_identity_single > 0.0 ||
            min_aligned_percent_single > 0.0;
        let filtering_pairs =
            min_aligned_length_pair > 0 ||
            min_percent_identity_pair > 0.0 ||
            min_aligned_percent_pair > 0.0;

        ReferenceSortedBamFilter {
            first_set: BTreeMap::new(),
            current_reference: -1,
            known_next_read: None,
            reader: reader,
            filter_single_reads: filtering_single,
            min_aligned_length_single: min_aligned_length_single,
            min_percent_identity_single: min_percent_identity_single,
            min_aligned_percent_single: min_aligned_percent_single,
            filter_pairs: filtering_pairs,
            min_aligned_length_pair: min_aligned_length_pair,
            min_percent_identity_pair: min_percent_identity_pair,
            min_aligned_percent_pair: min_aligned_percent_pair,
            num_detected_primary_alignments: 0,
        }
    }
}

impl ReferenceSortedBamFilter {
    pub fn read(&mut self, mut record: &mut bam::record::Record) -> Result<(), bam::ReadError> {
        // if doing only singles, remove filter them and return
        if self.filter_single_reads && !self.filter_pairs {
            loop {
                self.reader.read(&mut record)?;
                if !record.is_unmapped() &&
                    !record.is_secondary() &&
                    !record.is_supplementary() &&
                    single_read_passes_filter(
                        &record,
                        self.min_aligned_length_single,
                        self.min_percent_identity_single,
                        self.min_aligned_percent_single) {
                        return Ok(())
                    }
                // else this read shall not pass, try another
            }
        }

        // else doing pairs, so do as before except maybe filter out single reads too
        else {
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
                                // if filtering single and paired reads then
                                // both must pass QC, as well as the pair
                                // together.
                                if (!self.filter_single_reads ||
                                    (single_read_passes_filter(
                                        &record1,
                                        self.min_aligned_length_single,
                                        self.min_percent_identity_single,
                                        self.min_aligned_percent_single) &&
                                     single_read_passes_filter(
                                         &record,
                                         self.min_aligned_length_single,
                                         self.min_percent_identity_single,
                                         self.min_aligned_percent_single))) &&
                                    read_pair_passes_filter(
                                        &record,
                                        &record1,
                                        self.min_aligned_length_pair,
                                        self.min_percent_identity_pair,
                                        self.min_aligned_percent_pair) {

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
}

fn single_read_passes_filter(
    record: &bam::Record,
    min_aligned_length_single: u32,
    min_percent_identity_single: f32,
    min_aligned_percent_single: f32) -> bool {

    let edit_distance1 = match record.aux(b"NM") {
        Some(i) => i.integer(),
        None => {panic!("Alignment of read {:?} did not have an NM aux tag", record.qname())}
    };

    let mut aligned: u32 = 0;
    for cig in record.cigar().iter() {
        match cig {
            Cigar::Match(i) |
            Cigar::Ins(i) |
            Cigar::Diff(i) |
            Cigar::Equal(i) => {
                aligned += i;
            },
            _ => {}
        }
    }

    debug!("num_bases {}, distance {}, perc id {}, percent aligned {}",
           aligned, edit_distance1,
           1.0 - edit_distance1 as f32 / aligned as f32,
           aligned as f32 / record.seq().len() as f32);

    return aligned >= min_aligned_length_single &&
        aligned as f32 / record.seq().len() as f32 >= min_aligned_percent_single &&
        1.0 - edit_distance1 as f32 / aligned as f32 >= min_percent_identity_single
}

fn read_pair_passes_filter(
    record1: &bam::Record,
    record2: &bam::Record,
    min_aligned_length_pair: u32,
    min_percent_identity_pair: f32,
    min_aligned_percent_pair: f32) -> bool {

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
            Cigar::Match(i) |
            Cigar::Ins(i) |
            Cigar::Diff(i) |
            Cigar::Equal(i) => {
                aligned_length1 = aligned_length1 + i;
            },
            _ => {}
        }
    }
    let mut aligned_length2: u32 = 0;
    for cig in record2.cigar().iter() {
        match cig {
            Cigar::Match(i) |
            Cigar::Ins(i) |
            Cigar::Diff(i) |
            Cigar::Equal(i) => {
                aligned_length2 = aligned_length2 + i;
            },
            _ => {}
        }
    }

    let aligned = aligned_length1 + aligned_length2;
    debug!("num_bases {} {}, edit distances {} {}, perc id {}, percent aligned {}",
           aligned_length1, aligned_length2, edit_distance1, edit_distance2,
           1.0 - ((edit_distance1 + edit_distance2) as f32 / aligned as f32),
           aligned as f32 / ((record1.seq().len() + record2.seq().len()) as f32));

    return aligned >= min_aligned_length_pair &&
        aligned as f32 / (record1.seq().len() + record2.seq().len()) as f32 >= min_aligned_percent_pair &&
        1.0 - ((edit_distance1 + edit_distance2) as f32 / aligned as f32) >= min_percent_identity_pair
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hello_world(){
        let reader = bam::Reader::from_path(
            &"tests/data/7seqs.reads_for_seq1_and_seq2.bam").unwrap();
        let mut sorted = ReferenceSortedBamFilter::new(
            reader, 0,0.0,0.0, 90, 0.99, 0.0);
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
        let reader = bam::Reader::from_path(
            &"tests/data/2seqs.bad_read.1.bam").unwrap();
        let mut sorted = ReferenceSortedBamFilter::new(
            reader, 0,0.0,0.0, 250, 0.99, 0.0); // perc too high
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

        let reader = bam::Reader::from_path(
            &"tests/data/2seqs.bad_read.1.bam").unwrap();
        let mut sorted = ReferenceSortedBamFilter::new(
            reader, 0,0.0,0.0, 300, 0.98, 0.0); // aligned length too high
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

        let reader = bam::Reader::from_path(
            &"tests/data/2seqs.bad_read.1.with_extra.bam").unwrap();
        let mut sorted = ReferenceSortedBamFilter::new(
            reader, 0,0.0,0.0, 0, 0.98, 0.94); // aligned percent too high
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

        let reader = bam::Reader::from_path(
            &"tests/data/2seqs.bad_read.1.bam").unwrap();
        let mut sorted = ReferenceSortedBamFilter::new(
            reader, 0,0.0,0.0, 299, 0.98, 0.0);
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

    #[test]
    fn test_filter_single_reads(){
        let reader = bam::Reader::from_path(
            &"tests/data/2seqs.bad_read.1.bam").unwrap();
        let mut sorted = ReferenceSortedBamFilter::new(
            reader, 0, 0.99, 0.0, 0,0.0,0.0); // perc too high
        assert_eq!(true, sorted.filter_single_reads);
        assert_eq!(false, sorted.filter_pairs);
        let queries = vec!["2",
                           "3",
                           "4",
                           "1",];
        let mut record = bam::record::Record::new();
        println!("query: {:?}", queries);
        for i in queries {
            sorted.read(&mut record).expect("");
            assert_eq!(i, str::from_utf8(record.qname()).unwrap());
        }
    }

    #[test]
    fn test_filter_single_and_paired_reads(){
        let reader = bam::Reader::from_path(
            &"tests/data/2seqs.bad_read.1.bam").unwrap();
        let mut sorted = ReferenceSortedBamFilter::new(
            reader, 0, 0.95, 0.0, 300,0.0,0.0); // perc OK, but pair fails on length
        assert_eq!(true, sorted.filter_single_reads);
        assert_eq!(true, sorted.filter_pairs);
        let queries = vec!["2",
                           "2",
                           "3",
                           "3",
                           "4",
                           "4"];
        let mut record = bam::record::Record::new();
        println!("query: {:?}", queries);
        for i in queries {
            sorted.read(&mut record).expect("");
            assert_eq!(i, str::from_utf8(record.qname()).unwrap());
        }
    }
}
