use std::rc::Rc;
use std::str;

use rust_htslib::bam;
use std::collections::BTreeMap;
use rust_htslib::bam::Read;

struct ReferenceSortedBamFilter<'a> {
    first_set: BTreeMap<Rc<String>, bam::Record>,
    current_reference: i32,
    known_next_read: Option<bam::Record>,
    records: &'a mut bam::Records<'a, bam::Reader>,

}

impl<'a> ReferenceSortedBamFilter<'a> {
    fn new(records: &'a mut bam::Records<'a, bam::Reader>) -> ReferenceSortedBamFilter<'a> {
        ReferenceSortedBamFilter {
            first_set: BTreeMap::new(),
            current_reference: -1,
            known_next_read: None,
            records: records
        }
    }
}

impl<'a> Iterator for ReferenceSortedBamFilter<'a> {
    type Item = &'a bam::Record;

    fn next(&mut self) -> Option<Self::Item> {
        match self.known_next_read {
            Some(ref read) => {
                // Previously found a read pair, return the second of that pair.
                let read = Some(read);
                self.known_next_read = None;
                return read
            },
            None => {
                loop {
                    // read one of pair
                    match self.records.next() {
                        None => return None,
                        Some(record_result) => {
                            let record = record_result.expect("BAM read error");

                            // TODO: make usage ensure flag_filtering when mapping
                            if record.is_secondary() ||
                                record.is_supplementary() ||
                                !record.is_proper_pair() {
                                    continue;
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
                                    self.first_set.insert(Rc::new(qname), record.clone());
                                    return Some(&record)
                                }
                                // pairs from different contigs are ignored.
                            }
                            else { // Second read in insert
                                let qname = String::from(str::from_utf8(record.qname())
                                                         .expect("UTF8 error in conversion of read name"));
                                match self.first_set.get(&qname) {
                                    Some(record1) => {
                                        // if passes %ID and length thresholds, TODO
                                        self.first_set.remove(&qname);
                                        self.known_next_read = Some(record);
                                        return Some(record1)
                                    },
                                    // if pair is not in first read set, ignore it
                                    None => {}
                                }

                            }
                        }
                    }
                }
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hello_world(){
        let reader = bam::Reader::from_path(
            &"test/data/7seqs.reads_for_seq1_and_seq2.bam").unwrap().records();
        let sorted = ReferenceSortedBamFilter::new(
            &mut reader);
        println!("{:?}",sorted.next());
    }
}
