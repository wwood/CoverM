use std::process;

use tempfile::NamedTempFile;

use bam_generator::MappingProgram;
use mapping_index_maintenance::check_reference_existence;

#[derive(Clone, Debug)]
pub enum ReadFormat {
    Coupled,
    Interleaved,
    Single,
}

#[derive(Debug)]
pub struct MappingParameters<'a> {
    references: Vec<&'a str>,
    threads: u16,
    read1: Vec<&'a str>,
    read2: Vec<&'a str>,
    interleaved: Vec<&'a str>,
    unpaired: Vec<&'a str>,
    iter_reference_index: usize,
    mapping_options: Option<&'a str>,
}

impl<'a> MappingParameters<'a> {
    pub fn generate_from_clap(
        m: &'a clap::ArgMatches,
        mapping_program: MappingProgram,
        reference_tempfile: &'a Option<NamedTempFile>,
    ) -> MappingParameters<'a> {
        let mut read1: Vec<_> = vec![];
        let mut read2: Vec<_> = vec![];
        let mut interleaved: Vec<_> = vec![];
        let mut unpaired: Vec<_> = vec![];

        if m.contains_id("read1") {
            read1 = m
                .get_many::<String>("read1")
                .unwrap()
                .map(|s| s.as_str())
                .collect();
            read2 = m
                .get_many::<String>("read2")
                .unwrap()
                .map(|s| s.as_str())
                .collect();
            if read1.len() != read2.len() {
                error!(
                    "When specifying paired reads with the -1 and -2 flags, \
                        there must be equal numbers specified. Instead found \
                        {} and {} respectively",
                    read1.len(),
                    read2.len()
                );
                process::exit(1);
            }
        }

        // Parse --coupled
        if m.contains_id("coupled") {
            let coupled: Vec<_> = m
                .get_many::<String>("coupled")
                .unwrap()
                .map(|s| s.as_str())
                .collect();
            if coupled.len() % 2 != 0 {
                error!(
                    "The --coupled flag must be set with pairs of read \
                     sets, but an odd number ({}) was specified",
                    coupled.len()
                );
                process::exit(1);
            }
            let mut i = 0;
            while i < coupled.len() {
                read1.push(coupled[i]);
                read2.push(coupled[i + 1]);
                i += 2;
            }
        }

        if m.contains_id("interleaved") {
            interleaved = m
                .get_many::<String>("interleaved")
                .unwrap()
                .map(|s| s.as_str())
                .collect();
        }
        if m.contains_id("single") {
            unpaired = m
                .get_many::<String>("single")
                .unwrap()
                .map(|s| s.as_str())
                .collect();
        }

        match mapping_program {
            MappingProgram::MINIMAP2_ONT
            | MappingProgram::MINIMAP2_PB
            | MappingProgram::MINIMAP2_HIFI => {
                if !read1.is_empty() || !interleaved.is_empty() {
                    error!(
                        "Paired-end read input specified to be mapped \
                        with minimap2-ont, minimap2-pb, or minimap2-hifi which is presumably \
                        incorrect. Mapping paired reads can be run via \
                        minimap2-no-params if -ont or -pb mapping \
                        is desired."
                    );
                    process::exit(1);
                }
            }
            MappingProgram::X_MAPPER => {
                if !interleaved.is_empty() {
                    error!(
                        "Interleaved read input specified to be mapped with x-mapper which is not supported."
                    );
                    process::exit(1);
                }
            }
            _ => {}
        }

        let mapping_parameters_arg = match mapping_program {
            MappingProgram::BWA_MEM | MappingProgram::BWA_MEM2 => "bwa-params",
            MappingProgram::MINIMAP2_SR
            | MappingProgram::MINIMAP2_ONT
            | MappingProgram::MINIMAP2_HIFI
            | MappingProgram::MINIMAP2_PB
            | MappingProgram::MINIMAP2_NO_PRESET => "minimap2-params",
            MappingProgram::STROBEALIGN => "strobealign-params",
            MappingProgram::X_MAPPER => "x-mapper-params",
        };
        let mapping_options = match m.contains_id(mapping_parameters_arg) {
            true => {
                let params = m.get_one::<String>(mapping_parameters_arg);
                params
            }
            false => None,
        };
        debug!("Setting mapper {mapping_program:?} options as '{mapping_options:?}'");

        MappingParameters {
            references: match reference_tempfile {
                Some(r) => vec![r.path().to_str().unwrap()],
                None => match m.get_many::<String>("reference") {
                    Some(refs) => refs
                        .collect::<Vec<_>>()
                        .into_iter()
                        .map(|r| {
                            check_reference_existence(r, &mapping_program);
                            r.as_str()
                        })
                        .collect(),
                    None => vec![],
                },
            },
            threads: *m.get_one::<u16>("threads").unwrap(),
            read1,
            read2,
            interleaved,
            unpaired,
            iter_reference_index: 0,
            mapping_options: mapping_options.map(|x| &**x),
        }
    }

    // Return a Vec of str + Option<str> where each entry is a read pair or
    // single with None.
    pub fn readsets(&self) -> Vec<(&str, Option<&str>)> {
        let mut to_return: Vec<(&str, Option<&str>)> = vec![];

        for (r1, r2) in self.read1.iter().zip(self.read2.iter()) {
            to_return.push((r1, Some(r2)))
        }
        for s in self.unpaired.iter() {
            to_return.push((s, None))
        }
        to_return
    }
}

pub struct SingleReferenceMappingParameters<'a> {
    pub reference: &'a str,
    threads: u16,
    read1: Vec<&'a str>,
    read2: Vec<&'a str>,
    interleaved: Vec<&'a str>,
    unpaired: Vec<&'a str>,
    mapping_options: Option<&'a str>,

    iter_read_pair_index: usize,
    iter_interleaved_index: usize,
    iter_unpaired_index: usize,
}

impl<'a> SingleReferenceMappingParameters<'a> {
    pub fn len(&self) -> usize {
        self.read1.len() + self.interleaved.len() + self.unpaired.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<'a> Iterator for MappingParameters<'a> {
    type Item = SingleReferenceMappingParameters<'a>;

    fn next(&mut self) -> Option<SingleReferenceMappingParameters<'a>> {
        if self.iter_reference_index < self.references.len() {
            let i = self.iter_reference_index;
            self.iter_reference_index += 1;
            Some(SingleReferenceMappingParameters {
                reference: self.references[i],
                threads: self.threads,
                read1: self.read1.clone(),
                read2: self.read2.clone(),
                interleaved: self.interleaved.clone(),
                unpaired: self.unpaired.clone(),
                mapping_options: self.mapping_options,
                iter_read_pair_index: 0,
                iter_interleaved_index: 0,
                iter_unpaired_index: 0,
            })
        } else {
            None
        }
    }
}

impl<'a> Iterator for SingleReferenceMappingParameters<'a> {
    type Item = OneSampleMappingParameters<'a>;

    fn next(&mut self) -> Option<OneSampleMappingParameters<'a>> {
        if self.iter_read_pair_index < self.read1.len() {
            let i = self.iter_read_pair_index;
            self.iter_read_pair_index += 1;
            Some(OneSampleMappingParameters {
                reference: self.reference,
                read_format: ReadFormat::Coupled,
                read1: self.read1[i],
                read2: Some(self.read2[i]),
                threads: self.threads,
                mapping_options: self.mapping_options,
            })
        } else if self.iter_interleaved_index < self.interleaved.len() {
            let i = self.iter_interleaved_index;
            self.iter_interleaved_index += 1;
            Some(OneSampleMappingParameters {
                reference: self.reference,
                read_format: ReadFormat::Interleaved,
                read1: self.interleaved[i],
                read2: None,
                threads: self.threads,
                mapping_options: self.mapping_options,
            })
        } else if self.iter_unpaired_index < self.unpaired.len() {
            let i = self.iter_unpaired_index;
            self.iter_unpaired_index += 1;
            Some(OneSampleMappingParameters {
                reference: self.reference,
                read_format: ReadFormat::Single,
                read1: self.unpaired[i],
                read2: None,
                threads: self.threads,
                mapping_options: self.mapping_options,
            })
        } else {
            None
        }
    }
}

#[derive(Debug)]
pub struct OneSampleMappingParameters<'a> {
    pub reference: &'a str,
    pub read_format: ReadFormat,
    pub read1: &'a str,
    pub read2: Option<&'a str>,
    pub threads: u16,
    pub mapping_options: Option<&'a str>,
}
