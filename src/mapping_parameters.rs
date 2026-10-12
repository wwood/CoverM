use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

use tempfile::NamedTempFile;

use bam_generator::is_long_read_only_preset;
use bam_generator::MappingProgram;
use mapping_index_maintenance::check_reference_existence;

#[derive(Clone, Debug)]
pub enum ReadFormat {
    Coupled,
    Interleaved,
    Single,
}

/// Return the read file paths specified directly with `arg_id` (e.g.
/// "single"), followed by those listed in files given with `<arg_id>-list`,
/// one per line. For "coupled", each line of a list file must contain a
/// forward and reverse path separated by a tab, and both are returned in
/// order, as if they were specified with --coupled.
pub fn read_paths_from_clap<'a>(m: &'a clap::ArgMatches, arg_id: &str) -> Vec<&'a str> {
    let mut paths: Vec<&'a str> = match m.get_many::<String>(arg_id) {
        Some(values) => values.map(|s| s.as_str()).collect(),
        None => vec![],
    };
    let list_arg_id = format!("{arg_id}-list");
    if let Some(list_files) = m.get_many::<String>(&list_arg_id) {
        for list_file in list_files {
            let file = File::open(list_file).unwrap_or_else(|e| {
                error!("Failed to open --{list_arg_id} file '{list_file}': {e}");
                process::exit(1);
            });
            for (i, line) in BufReader::new(file).lines().enumerate() {
                let line = line.unwrap_or_else(|e| {
                    error!("Failed to read --{list_arg_id} file '{list_file}': {e}");
                    process::exit(1);
                });
                let line = line.trim();
                if line.is_empty() {
                    continue;
                }
                let fields: Vec<&str> = line.split('\t').map(|f| f.trim()).collect();
                let expected_fields = if arg_id == "coupled" { 2 } else { 1 };
                if fields.len() != expected_fields || fields.iter().any(|f| f.is_empty()) {
                    error!(
                        "Line {} of --{list_arg_id} file '{list_file}' should contain {}, \
                        but found '{line}'",
                        i + 1,
                        if expected_fields == 2 {
                            "a forward and reverse read path separated by a tab"
                        } else {
                            "a single read path"
                        }
                    );
                    process::exit(1);
                }
                // Leak so the paths live as long as those from clap. These are
                // read once at startup, so the memory cost is negligible.
                for field in fields {
                    paths.push(field.to_string().leak());
                }
            }
        }
    }
    paths
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

        if m.contains_id("read1") || m.contains_id("read1-list") {
            read1 = read_paths_from_clap(m, "read1");
            read2 = read_paths_from_clap(m, "read2");
            if read1.len() != read2.len() {
                error!(
                    "When specifying paired reads with the -1 and -2 flags \
                        (or --read1-list and --read2-list), there must be equal \
                        numbers specified. Instead found {} and {} respectively",
                    read1.len(),
                    read2.len()
                );
                process::exit(1);
            }
        }

        // Parse --coupled
        if m.contains_id("coupled") || m.contains_id("coupled-list") {
            let coupled = read_paths_from_clap(m, "coupled");
            if !coupled.len().is_multiple_of(2) {
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

        let interleaved = read_paths_from_clap(m, "interleaved");
        let unpaired = read_paths_from_clap(m, "single");

        if is_long_read_only_preset(mapping_program) {
            if !read1.is_empty() || !interleaved.is_empty() {
                error!(
                    "Paired-end read input specified to be mapped \
                    with a long-read preset (e.g. minimap2-ont, minimap2-lr-hq, \
                    rammap-ont, rammap-lr-hq) which is presumably \
                    incorrect. Mapping paired reads can be run via \
                    minimap2-no-preset if -ont or -pb mapping \
                    is desired."
                );
                process::exit(1);
            }
        } else if matches!(mapping_program, MappingProgram::MINIBWA) && !interleaved.is_empty() {
            error!(
                "Interleaved read input was specified to be mapped with minibwa, \
                but minibwa has no interleaved-pairing option. Supply the reads as \
                coupled pairs (with -1/-2 or --coupled) instead."
            );
            process::exit(1);
        }

        let mapping_parameters_arg = match mapping_program {
            MappingProgram::BWA_MEM | MappingProgram::BWA_MEM2 => "bwa-params",
            MappingProgram::MINIMAP2_SR
            | MappingProgram::MINIMAP2_ONT
            | MappingProgram::MINIMAP2_HIFI
            | MappingProgram::MINIMAP2_PB
            | MappingProgram::MINIMAP2_LR_HQ
            | MappingProgram::MINIMAP2_NO_PRESET => "minimap2-params",
            MappingProgram::STROBEALIGN => "strobealign-params",
            MappingProgram::MINIBWA => "minibwa-params",
            MappingProgram::RAMMAP_SR
            | MappingProgram::RAMMAP_ONT
            | MappingProgram::RAMMAP_PB
            | MappingProgram::RAMMAP_HIFI
            | MappingProgram::RAMMAP_LR_HQ
            | MappingProgram::RAMMAP_NO_PRESET => "rammap-params",
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
