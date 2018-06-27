use std;
use rust_htslib::bam;

struct ReferenceSortedBamFilter {
    bam_reader: bam::Reader,
    first_set: BTreeSet<String, bam::Record>,
    current_reference: u32,
    known_next_read: Option<bam::Record>,
}

impl Iterator for ReferenceSortedBamFilter {
    type Item = bam::Record;

    fn next(&mut self) -> Option<Self::Item> {
        match self.known_next_read {
            Some(read) => {
                // Previously found a read pair, return the second of that pair.
                let read = Some(read);
                self.known_next_read = None;
            },
            None => {
                // read one of pair
                let record
                let read = self.next
                // if a new reference ID is encountered, instantiate a new first read set
                // if this is a first read
                // if tlen is +ve and < threshold
                // add to first read set
                // else
                // if pair is not in first read set, ignore it
                // if passes %ID and length thresholds, send it via the channel
            }
        }
    }
}

