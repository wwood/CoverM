use needletail::parse_fastx_file;
use std;
use std::collections::HashSet;
use std::io::Read;
use std::path::Path;
use std::process;

use bam_generator::MappingProgram;
use CONCATENATED_FASTA_FILE_SEPARATOR;

use tempdir::TempDir;
use tempfile::{Builder, NamedTempFile};

pub trait MappingIndex {
    fn index_path(&self) -> &String;
    fn command_prefix(&self) -> &str {
        ""
    }
}

pub struct VanillaIndexStruct {
    index_path_internal: String,
}
impl VanillaIndexStruct {
    pub fn new(reference_path: &str) -> VanillaIndexStruct {
        VanillaIndexStruct {
            index_path_internal: reference_path.to_string(),
        }
    }
}
impl MappingIndex for VanillaIndexStruct {
    fn index_path(&self) -> &String {
        &self.index_path_internal
    }
}

pub struct TemporaryIndexStruct {
    #[allow(dead_code)] // field is never used, it just needs to be kept in scope.
    tempdir: TempDir,
    index_path_internal: String,
}

impl TemporaryIndexStruct {
    pub fn new(
        mapping_program: MappingProgram,
        reference_path: &str,
        num_threads: Option<u16>,
        index_creation_options: Option<&str>,
    ) -> TemporaryIndexStruct {
        // Generate a BWA/minimap index in a temporary directory, where the
        // temporary directory does not go out of scope until the struct does.
        let td =
            TempDir::new("coverm-mapping-index").expect("Unable to create temporary directory");
        let index_path = std::path::Path::new(td.path()).join(
            std::path::Path::new(reference_path)
                .file_name()
                .expect("Failed to glean file stem from reference DB. Strange."),
        );

        info!(
            "Generating {:?} index for {} ..",
            mapping_program, reference_path
        );
        let mut cmd = match mapping_program {
            MappingProgram::BWA_MEM => std::process::Command::new("bwa"),
            MappingProgram::BWA_MEM2 => std::process::Command::new("bwa-mem2"),
            MappingProgram::MINIMAP2_SR
            | MappingProgram::MINIMAP2_ONT
            | MappingProgram::MINIMAP2_PB
            | MappingProgram::MINIMAP2_HIFI
            | MappingProgram::MINIMAP2_NO_PRESET => std::process::Command::new("minimap2"),
            MappingProgram::STROBEALIGN => std::process::Command::new("strobealign"),
        };
        match &mapping_program {
            MappingProgram::BWA_MEM | MappingProgram::BWA_MEM2 => {
                cmd.arg("index")
                    .arg("-p")
                    .arg(&index_path)
                    .arg(reference_path);
            }
            MappingProgram::MINIMAP2_SR
            | MappingProgram::MINIMAP2_ONT
            | MappingProgram::MINIMAP2_HIFI
            | MappingProgram::MINIMAP2_PB
            | MappingProgram::MINIMAP2_NO_PRESET => {
                match &mapping_program {
                    MappingProgram::MINIMAP2_SR => {
                        cmd.arg("-x").arg("sr");
                    }
                    MappingProgram::MINIMAP2_ONT => {
                        cmd.arg("-x").arg("map-ont");
                    }
                    MappingProgram::MINIMAP2_HIFI => {
                        cmd.arg("-x").arg("map-hifi");
                    }
                    MappingProgram::MINIMAP2_PB => {
                        cmd.arg("-x").arg("map-pb");
                    }
                    MappingProgram::MINIMAP2_NO_PRESET
                    | MappingProgram::BWA_MEM
                    | MappingProgram::BWA_MEM2
                    | MappingProgram::STROBEALIGN => {}
                };
                if let Some(t) = num_threads {
                    cmd.arg("-t").arg(format!("{}", t));
                }
                cmd.arg("-d").arg(&index_path).arg(reference_path);
            }
            MappingProgram::STROBEALIGN => {
                warn!("STROBEALIGN pre-indexing is not supported currently, so skipping index generation.");
            }
        };
        if let Some(params) = index_creation_options {
            for s in params.split_whitespace() {
                cmd.arg(s);
            }
        };
        // Some BWA versions output log info to stdout. Ignore this.
        cmd.stdout(std::process::Stdio::piped());
        cmd.stderr(std::process::Stdio::piped());
        debug!("Running DB indexing command: {:?}", cmd);

        let mut process = cmd
            .spawn()
            .unwrap_or_else(|_| panic!("Failed to start {:?} index process", mapping_program));
        let es = process.wait().unwrap_or_else(|_| {
            panic!(
                "Failed to glean exitstatus from failing {:?} index process",
                mapping_program
            )
        });
        if !es.success() {
            error!("Error when running {:?} index process.", mapping_program);
            let mut err = String::new();
            process
                .stderr
                .unwrap_or_else(|| {
                    panic!(
                        "Failed to grab stderr from failed {:?} index process",
                        mapping_program
                    )
                })
                .read_to_string(&mut err)
                .expect("Failed to read stderr into string");
            error!("The STDERR was: {:?}", err);
            error!("Cannot continue after {:?} index failed.", mapping_program);
            process::exit(1);
        }
        info!("Finished generating {:?} index.", mapping_program);
        return TemporaryIndexStruct {
            index_path_internal: index_path.to_string_lossy().to_string(),
            tempdir: td,
        };
    }
}
impl MappingIndex for TemporaryIndexStruct {
    fn index_path(&self) -> &String {
        &self.index_path_internal
    }
}
impl Drop for TemporaryIndexStruct {
    fn drop(&mut self) {
        debug!(
            "Dropping index tempdir {}",
            self.tempdir.path().to_string_lossy()
        )
    }
}

fn check_for_bwa_index_existence(reference_path: &str, mapping_program: &MappingProgram) -> bool {
    let bwa_extensions = match mapping_program {
        MappingProgram::BWA_MEM => vec!["amb", "ann", "bwt", "pac", "sa"],
        MappingProgram::BWA_MEM2 => vec!["0123", "amb", "ann", "bwt.2bit.64", "pac"],
        _ => unreachable!(),
    };
    let num_extensions = bwa_extensions.len();
    let mut num_existing: usize = 0;
    for extension in bwa_extensions {
        if std::path::Path::new(&format!("{}.{}", reference_path, extension)).exists() {
            num_existing += 1;
        }
    }
    if num_existing == 0 {
        false
    } else if num_existing == num_extensions {
        return true;
    } else {
        error!("BWA index appears to be incomplete, cannot continue.");
        process::exit(1);
    }
}

/// Check that a reference exists, or that a corresponding index exists.
pub fn check_reference_existence(reference_path: &str, mapping_program: &MappingProgram) {
    let ref_path = std::path::Path::new(reference_path);
    match mapping_program {
        MappingProgram::BWA_MEM | MappingProgram::BWA_MEM2 => {
            if check_for_bwa_index_existence(reference_path, mapping_program) {
                return;
            }
        }
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_HIFI
        | MappingProgram::MINIMAP2_PB
        | MappingProgram::MINIMAP2_NO_PRESET
        | MappingProgram::STROBEALIGN => {}
    };

    if !ref_path.exists() {
        panic!(
            "The reference specified '{}' does not appear to exist",
            &reference_path
        );
    } else if !ref_path.is_file() {
        panic!(
            "The reference specified '{}' should be a file, not e.g. a directory",
            &reference_path
        );
    }
}

pub fn generate_bwa_index(
    reference_path: &str,
    index_creation_parameters: Option<&str>,
    mapping_program: MappingProgram,
) -> Box<dyn MappingIndex> {
    if check_for_bwa_index_existence(reference_path, &mapping_program) {
        info!("BWA index appears to be complete, so going ahead and using it.");
        Box::new(VanillaIndexStruct::new(reference_path))
    } else {
        Box::new(TemporaryIndexStruct::new(
            mapping_program,
            reference_path,
            None,
            index_creation_parameters,
        ))
    }
}

pub fn generate_minimap2_index(
    reference_path: &str,
    num_threads: Option<u16>,
    index_creation_parameters: Option<&str>,
    mapping_program: MappingProgram,
) -> Box<dyn MappingIndex> {
    Box::new(TemporaryIndexStruct::new(
        mapping_program,
        reference_path,
        num_threads,
        index_creation_parameters,
    ))
}

pub fn generate_concatenated_fasta_file(fasta_file_paths: &Vec<String>) -> NamedTempFile {
    let tmpfile: NamedTempFile = Builder::new()
        .prefix("coverm-concatenated-fasta")
        .tempfile()
        .unwrap();
    let mut something_written_at_all = false;
    {
        // scope so writer, which borrows tmpfile, goes out of scope.
        let mut writer = bio::io::fasta::Writer::new(&tmpfile);
        let mut genome_names: HashSet<String> = HashSet::new();

        // NOTE: A lot of this code is shared with genome_parsing#read_genome_fasta_files
        for file in fasta_file_paths {
            let mut something_written = false;
            let path = std::path::Path::new(file);
            let mut reader = parse_fastx_file(path)
                .unwrap_or_else(|_| panic!("Unable to read fasta file {}", file));

            // Remove .gz .bz .xz from file names if present
            let mut genome_name1 =
                String::from(path.to_str().expect("File name string conversion problem"));
            if let Some(i) = genome_name1.rfind(".gz") {
                genome_name1.truncate(i);
            } else if let Some(i) = genome_name1.rfind(".bz") {
                genome_name1.truncate(i);
            } else if let Some(i) = genome_name1.rfind(".xz") {
                genome_name1.truncate(i);
            }
            let path1 = Path::new(&genome_name1);

            let genome_name = String::from(
                path1
                    .file_stem()
                    .expect("Problem while determining file stem")
                    .to_str()
                    .expect("File name string conversion problem"),
            );
            if genome_names.contains(&genome_name) {
                error!("The genome name {} was derived from >1 file", genome_name);
                process::exit(1);
            }
            while let Some(record) = reader.next() {
                let record_expected = record
                    .unwrap_or_else(|_| panic!("Failed to parse record in fasta file {:?}", path));

                if record_expected.format() != needletail::parser::Format::Fasta {
                    panic!(
                        "File {:?} is not a fasta file, but a {:?}",
                        path,
                        record_expected.format()
                    );
                }

                something_written = true;
                something_written_at_all = true;
                let contig_name = String::from(
                    std::str::from_utf8(record_expected.id())
                        .expect("UTF-8 conversion problem in contig name"),
                );
                writer
                    .write(
                        &format!(
                            "{}{}{}",
                            genome_name,
                            CONCATENATED_FASTA_FILE_SEPARATOR,
                            match contig_name.split_once(' ') {
                                Some((contig, _)) => contig.to_string(),
                                None => contig_name,
                            }
                        ),
                        None,
                        &record_expected.seq(),
                    )
                    .unwrap()
            }
            genome_names.insert(genome_name);
            if !something_written {
                error!(
                    "FASTA file {} appears to be empty as no sequences were contained in it",
                    file
                );
                process::exit(1);
            }
        }
    }
    if !something_written_at_all {
        error!("Concatenated FASTA file to use as a reference is empty");
        process::exit(1);
    }
    tmpfile
}

pub struct PregeneratedStrobealignIndexStruct {
    index_path_internal: String,
}
impl PregeneratedStrobealignIndexStruct {
    pub fn new(reference_path: &str) -> PregeneratedStrobealignIndexStruct {
        PregeneratedStrobealignIndexStruct {
            index_path_internal: reference_path.to_string(),
        }
    }
}
impl MappingIndex for PregeneratedStrobealignIndexStruct {
    fn index_path(&self) -> &String {
        &self.index_path_internal
    }

    fn command_prefix(&self) -> &str {
        "--use-index"
    }
}
