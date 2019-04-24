
use tempfile::NamedTempFile;

#[derive(Clone)]
pub enum ReadFormat {
    Coupled,
    Interleaved,
    Single,
}

pub struct MappingParameters<'a> {
    references: Vec<&'a str>,
    threads: u16,
    read1: Vec<&'a str>,
    read2: Vec<&'a str>,
    interleaved: Vec<&'a str>,
    unpaired: Vec<&'a str>,
    iter_reference_index: usize,
    bwa_options: Option<&'a str>
}

impl<'a> MappingParameters<'a> {
    pub fn generate_from_clap(
        m: &'a clap::ArgMatches,
        reference_tempfile: &'a Option<NamedTempFile>)
        -> MappingParameters<'a> {

        let mut read1: Vec<&str> = vec!();
        let mut read2: Vec<&str> = vec!();
        let mut interleaved: Vec<&str> = vec!();
        let mut unpaired: Vec<&str> = vec!();

        if m.is_present("read1") {
            read1 = m.values_of("read1").unwrap().collect();
            read2 = m.values_of("read2").unwrap().collect();
            if read1.len() != read2.len() {
                panic!("When specifying paired reads with the -1 and -2 flags, there must be equal numbers specified. Instead found {} and {} respectively", read1.len(), read2.len())
            }
        }

        // Parse --coupled
        if m.is_present("coupled") {
            let coupled: Vec<&str> = m.values_of("coupled").unwrap().collect();
            if coupled.len() % 2 != 0 {
                panic!(
                    "The --coupled flag must be set with pairs of read sets, but an odd number ({}) was specified",
                    coupled.len()
                )
            }
            let mut i = 0;
            while i < coupled.len() {
                read1.push(coupled[i]);
                read2.push(coupled[i+1]);
                i += 2;
            }
        }

        if m.is_present("interleaved") {
            interleaved = m.values_of("interleaved").unwrap().collect();
        }
        if m.is_present("single") {
            unpaired = m.values_of("single").unwrap().collect();
        }

        let bwa_options = match m.is_present("bwa-params") {
            true => {
                let params = m.value_of("bwa-params");
                params
            }
            false => None
        };
        debug!("Setting BWA options as '{:?}'", bwa_options);

        return MappingParameters {
            references: match reference_tempfile {
                Some(r) => vec!(r.path().to_str().unwrap()),
                None => m.values_of("reference").unwrap().collect()
            },
            threads: m.value_of("threads").unwrap().parse::<u16>()
                .expect("Failed to convert threads argument into integer"),
            read1: read1,
            read2: read2,
            interleaved: interleaved,
            unpaired: unpaired,
            iter_reference_index: 0,
            bwa_options: bwa_options,
        }
    }
}

pub struct SingleReferenceMappingParameters<'a> {
    pub reference: &'a str,
    threads: u16,
    read1: Vec<&'a str>,
    read2: Vec<&'a str>,
    interleaved: Vec<&'a str>,
    unpaired: Vec<&'a str>,
    bwa_options: Option<&'a str>,

    iter_read_pair_index: usize,
    iter_interleaved_index: usize,
    iter_unpaired_index: usize,
}

impl<'a> Iterator for MappingParameters<'a> {
    type Item = SingleReferenceMappingParameters<'a>;

    fn next(&mut self) -> Option<SingleReferenceMappingParameters<'a>> {
        if self.iter_reference_index < self.references.len() {
            let i = self.iter_reference_index;
            self.iter_reference_index += 1;
            return Some(SingleReferenceMappingParameters {
                reference: self.references[i],
                threads: self.threads,
                read1: self.read1.clone(),
                read2: self.read2.clone(),
                interleaved: self.interleaved.clone(),
                unpaired: self.unpaired.clone(),
                bwa_options: self.bwa_options,
                iter_read_pair_index: 0,
                iter_interleaved_index: 0,
                iter_unpaired_index: 0,
            })
        } else {
            return None
        }
    }
}

impl<'a> Iterator for SingleReferenceMappingParameters<'a> {
    type Item = OneSampleMappingParameters<'a>;

    fn next(&mut self) -> Option<OneSampleMappingParameters<'a>> {
        if self.iter_read_pair_index < self.read1.len() {
            let i = self.iter_read_pair_index;
            self.iter_read_pair_index += 1;
            return Some(OneSampleMappingParameters {
                reference: self.reference,
                read_format: ReadFormat::Coupled,
                read1: self.read1[i],
                read2: Some(self.read2[i]),
                threads: self.threads,
                bwa_options: self.bwa_options,
            })
        } else if self.iter_interleaved_index < self.interleaved.len() {
            let i = self.iter_interleaved_index;
            self.iter_interleaved_index += 1;
            return Some(OneSampleMappingParameters {
                reference: self.reference,
                read_format: ReadFormat::Interleaved,
                read1: self.interleaved[i],
                read2: None,
                threads: self.threads,
                bwa_options: self.bwa_options,
            })
        } else if self.iter_unpaired_index < self.unpaired.len() {
            let i = self.iter_unpaired_index;
            self.iter_unpaired_index += 1;
            return Some(OneSampleMappingParameters {
                reference: self.reference,
                read_format: ReadFormat::Single,
                read1: self.unpaired[i],
                read2: None,
                threads: self.threads,
                bwa_options: self.bwa_options,
            })
        } else {
            return None
        }
    }
}

pub struct OneSampleMappingParameters<'a> {
    pub reference: &'a str,
    pub read_format: ReadFormat,
    pub read1: &'a str,
    pub read2: Option<&'a str>,
    pub threads: u16,
    pub bwa_options: Option<&'a str>
}
