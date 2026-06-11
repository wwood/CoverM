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

        run_index_command(
            mapping_program,
            reference_path,
            &index_path,
            num_threads,
            index_creation_options,
        );

        TemporaryIndexStruct {
            index_path_internal: index_path.to_string_lossy().to_string(),
            tempdir: td,
        }
    }
}

/// Build the index-generation command for the given mapping program, writing
/// the resulting index to `index_path`. Returns None for STROBEALIGN, which
/// does not support standalone pre-indexing in this manner.
fn build_index_command(
    mapping_program: MappingProgram,
    reference_path: &str,
    index_path: &Path,
    num_threads: Option<u16>,
    index_creation_options: Option<&str>,
) -> Option<std::process::Command> {
    let mut cmd = match mapping_program {
        MappingProgram::BWA_MEM => std::process::Command::new("bwa"),
        MappingProgram::BWA_MEM2 => std::process::Command::new("bwa-mem2"),
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_PB
        | MappingProgram::MINIMAP2_HIFI
        | MappingProgram::MINIMAP2_LR_HQ
        | MappingProgram::MINIMAP2_NO_PRESET => std::process::Command::new("minimap2"),
        MappingProgram::STROBEALIGN => {
            warn!("STROBEALIGN pre-indexing is not supported currently, so skipping index generation.");
            return None;
        }
    };
    match &mapping_program {
        MappingProgram::BWA_MEM | MappingProgram::BWA_MEM2 => {
            cmd.arg("index")
                .arg("-p")
                .arg(index_path)
                .arg(reference_path);
        }
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_HIFI
        | MappingProgram::MINIMAP2_PB
        | MappingProgram::MINIMAP2_LR_HQ
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
                MappingProgram::MINIMAP2_LR_HQ => {
                    cmd.arg("-x").arg("lr:hq");
                }
                MappingProgram::MINIMAP2_NO_PRESET
                | MappingProgram::BWA_MEM
                | MappingProgram::BWA_MEM2
                | MappingProgram::STROBEALIGN => {}
            };
            if let Some(t) = num_threads {
                cmd.arg("-t").arg(format!("{t}"));
            }
            cmd.arg("-d").arg(index_path).arg(reference_path);
        }
        MappingProgram::STROBEALIGN => unreachable!(),
    };
    if let Some(params) = index_creation_options {
        for s in params.split_whitespace() {
            cmd.arg(s);
        }
    };
    Some(cmd)
}

/// Run the index-generation command for the given mapping program, writing the
/// resulting index to `index_path`. Exits the process on failure.
fn run_index_command(
    mapping_program: MappingProgram,
    reference_path: &str,
    index_path: &Path,
    num_threads: Option<u16>,
    index_creation_options: Option<&str>,
) {
    info!("Generating {mapping_program:?} index for {reference_path} ..");
    let cmd = match build_index_command(
        mapping_program,
        reference_path,
        index_path,
        num_threads,
        index_creation_options,
    ) {
        Some(cmd) => cmd,
        None => return,
    };
    execute_index_command(cmd, mapping_program);
}

/// Spawn and wait on a pre-built index-generation command, exiting the process
/// on failure.
fn execute_index_command(mut cmd: std::process::Command, mapping_program: MappingProgram) {
    // Some BWA versions output log info to stdout. Ignore this.
    cmd.stdout(std::process::Stdio::piped());
    cmd.stderr(std::process::Stdio::piped());
    debug!("Running DB indexing command: {cmd:?}");

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
        error!("Error when running {mapping_program:?} index process.");
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
        error!("The STDERR was: {err:?}");
        error!("Cannot continue after {mapping_program:?} index failed.");
        process::exit(1);
    }
    info!("Finished generating {mapping_program:?} index.");
}

/// Generate a strobealign index for `reference_path` inside `output_directory`.
///
/// Unlike minimap2/BWA, `strobealign --create-index` writes its `.sti` index
/// alongside the reference FASTA (and the FASTA is still required at mapping
/// time, since strobealign reads the sequences from it). To keep the generated
/// database self-contained within `output_directory`, the reference is first
/// copied there and the index is created next to the copy. Returns the path to
/// the copied reference, which is what should be passed to
/// `--reference .. --strobealign-use-index`.
///
/// Strobealign indexes are read-length specific. The canonical read length can
/// be set by passing `-r <length>` via `index_creation_options`, or estimated
/// by passing a reads file there (e.g. `reads.fq`), matching the positional
/// argument that `strobealign --create-index` accepts.
fn generate_strobealign_persistent_index(
    reference_path: &str,
    output_directory: &str,
    num_threads: Option<u16>,
    index_creation_options: Option<&str>,
) -> String {
    let reference_file_name = std::path::Path::new(reference_path)
        .file_name()
        .expect("Failed to glean file name from reference path");
    let copied_reference = std::path::Path::new(output_directory).join(reference_file_name);

    if copied_reference == std::path::Path::new(reference_path) {
        // The reference already lives in the output directory (e.g. a
        // concatenated_genomes.fasta generated by `coverm makedb` from a genome
        // definition), so the database is already self-contained and the index
        // can be created in place.
        info!(
            "Reference {reference_path} is already inside the output directory, \
            so indexing it in place."
        );
    } else {
        info!(
            "Copying reference {} into {} so the strobealign database is self-contained ..",
            reference_path, output_directory
        );
        std::fs::copy(reference_path, &copied_reference).unwrap_or_else(|e| {
            panic!(
                "Failed to copy reference {} into output directory: {}",
                reference_path, e
            )
        });
    }

    info!("Generating STROBEALIGN index for {reference_path} ..");
    let mut cmd = std::process::Command::new("strobealign");
    if let Some(t) = num_threads {
        cmd.arg("-t").arg(format!("{t}"));
    }
    cmd.arg("--create-index").arg(&copied_reference);
    if let Some(params) = index_creation_options {
        for s in params.split_whitespace() {
            cmd.arg(s);
        }
    };
    execute_index_command(cmd, MappingProgram::STROBEALIGN);

    copied_reference.to_string_lossy().to_string()
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
        if std::path::Path::new(&format!("{reference_path}.{extension}")).exists() {
            num_existing += 1;
        }
    }
    if num_existing == 0 {
        false
    } else if num_existing == num_extensions {
        true
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
        | MappingProgram::MINIMAP2_LR_HQ
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

/// Short name used to label the kind of database generated for a given mapping
/// program, e.g. as a filename suffix in `coverm makedb`.
pub fn mapping_program_db_name(mapping_program: MappingProgram) -> &'static str {
    match mapping_program {
        MappingProgram::BWA_MEM => "bwa-mem",
        MappingProgram::BWA_MEM2 => "bwa-mem2",
        MappingProgram::MINIMAP2_SR => "minimap2-sr",
        MappingProgram::MINIMAP2_ONT => "minimap2-ont",
        MappingProgram::MINIMAP2_PB => "minimap2-pb",
        MappingProgram::MINIMAP2_HIFI => "minimap2-hifi",
        MappingProgram::MINIMAP2_LR_HQ => "minimap2-lr-hq",
        MappingProgram::MINIMAP2_NO_PRESET => "minimap2-no-preset",
        MappingProgram::STROBEALIGN => "strobealign",
    }
}

/// Generate a mapping index for `reference_path` in `output_directory` that is
/// persisted on disk (unlike [`TemporaryIndexStruct`]), so it can later be fed
/// back into `coverm contig`/`coverm genome`. Returns the path to the generated
/// database. Used by `coverm makedb`.
pub fn generate_persistent_index(
    mapping_program: MappingProgram,
    reference_path: &str,
    output_directory: &str,
    num_threads: Option<u16>,
    index_creation_options: Option<&str>,
) -> String {
    // Strobealign writes its index next to the reference (and needs the
    // reference at mapping time), so it is handled separately.
    if let MappingProgram::STROBEALIGN = mapping_program {
        return generate_strobealign_persistent_index(
            reference_path,
            output_directory,
            num_threads,
            index_creation_options,
        );
    }

    let reference_stem = std::path::Path::new(reference_path)
        .file_name()
        .expect("Failed to glean file name from reference path")
        .to_string_lossy();
    let db_program_name = mapping_program_db_name(mapping_program);

    // For minimap2 the index is a single .mmi file, while for BWA the supplied
    // path is used as a prefix for the several files BWA generates.
    let index_path = match mapping_program {
        MappingProgram::MINIMAP2_SR
        | MappingProgram::MINIMAP2_ONT
        | MappingProgram::MINIMAP2_PB
        | MappingProgram::MINIMAP2_HIFI
        | MappingProgram::MINIMAP2_LR_HQ
        | MappingProgram::MINIMAP2_NO_PRESET => std::path::Path::new(output_directory)
            .join(format!("{reference_stem}.{db_program_name}.mmi")),
        MappingProgram::BWA_MEM | MappingProgram::BWA_MEM2 => {
            std::path::Path::new(output_directory)
                .join(format!("{reference_stem}.{db_program_name}"))
        }
        MappingProgram::STROBEALIGN => unreachable!(),
    };

    run_index_command(
        mapping_program,
        reference_path,
        &index_path,
        num_threads,
        index_creation_options,
    );

    index_path.to_string_lossy().to_string()
}

/// Write a single concatenated FASTA file from the given per-genome FASTA
/// files to `writer`, prefixing every contig name with its genome name (the
/// file stem) and [`CONCATENATED_FASTA_FILE_SEPARATOR`] so the genome of each
/// contig can later be recovered. Returns whether any sequence was written.
fn write_concatenated_fasta<W: std::io::Write>(writer: W, fasta_file_paths: &[String]) -> bool {
    let mut writer = bio::io::fasta::Writer::new(writer);
    let mut genome_names: HashSet<String> = HashSet::new();
    let mut something_written_at_all = false;

    // NOTE: A lot of this code is shared with genome_parsing#read_genome_fasta_files
    for file in fasta_file_paths {
        let mut something_written = false;
        let path = std::path::Path::new(file);
        let mut reader =
            parse_fastx_file(path).unwrap_or_else(|_| panic!("Unable to read fasta file {}", file));

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
            error!("The genome name {genome_name} was derived from >1 file");
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
            error!("FASTA file {file} appears to be empty as no sequences were contained in it");
            process::exit(1);
        }
    }
    something_written_at_all
}

pub fn generate_concatenated_fasta_file(fasta_file_paths: &[String]) -> NamedTempFile {
    let tmpfile: NamedTempFile = Builder::new()
        .prefix("coverm-concatenated-fasta")
        .tempfile()
        .unwrap();
    // scope so writer, which borrows tmpfile, goes out of scope.
    let something_written_at_all = write_concatenated_fasta(&tmpfile, fasta_file_paths);
    if !something_written_at_all {
        error!("Concatenated FASTA file to use as a reference is empty");
        process::exit(1);
    }
    tmpfile
}

/// Like [`generate_concatenated_fasta_file`], but writes the concatenated FASTA
/// to a persistent file at `output_path` (rather than a temporary file). Used
/// by `coverm makedb` so a database built from a genome definition keeps the
/// reference FASTA around to be re-used as a `coverm genome` reference.
pub fn write_concatenated_fasta_file(fasta_file_paths: &[String], output_path: &Path) {
    let file = std::fs::File::create(output_path).unwrap_or_else(|e| {
        panic!(
            "Failed to create concatenated FASTA file {:?}: {}",
            output_path, e
        )
    });
    let something_written_at_all = write_concatenated_fasta(file, fasta_file_paths);
    if !something_written_at_all {
        error!("Concatenated FASTA file to use as a reference is empty");
        process::exit(1);
    }
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
