use std;
use std::io::Read;
use std::collections::HashSet;
use std::process;

use CONCATENATED_FASTA_FILE_SEPARATOR;
use bam_generator::MappingProgram;

use tempdir::TempDir;
use tempfile::NamedTempFile;


/// Actually a trait for all kinds of mapping indices, just too lazy to change
/// the name.
pub trait MappingIndex {
    fn index_path(&self) -> &String;
}

pub struct VanillaBwaIndexStuct {
    index_path_internal: String
}
impl VanillaBwaIndexStuct {
    pub fn new(reference_path: &str) -> VanillaBwaIndexStuct {
        return VanillaBwaIndexStuct {
            index_path_internal: reference_path.to_string()
        }
    }
}
impl MappingIndex for VanillaBwaIndexStuct {
    fn index_path(&self) -> &String {
        return &self.index_path_internal
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
        index_creation_options: Option<&str>)
        -> TemporaryIndexStruct {

        // Generate a BWA/minimap index in a temporary directory, where the
        // temporary directory does not go out of scope until the struct does.
        let td = TempDir::new("coverm-mapping-index")
            .expect("Unable to create temporary directory");
        let index_path = std::path::Path::new(td.path())
            .join(std::path::Path::new(reference_path).file_name()
                  .expect("Failed to glean file stem from reference DB. Strange."));

        info!("Generating {:?} index for {} ..", mapping_program, reference_path);
        let mut cmd = match mapping_program {
            MappingProgram::BWA_MEM => {
                std::process::Command::new("bwa")
            },
            MappingProgram::MINIMAP2_SR |
            MappingProgram::MINIMAP2_ONT |
            MappingProgram::MINIMAP2_PB |
            MappingProgram::MINIMAP2_NO_PRESET => {
                std::process::Command::new("minimap2")
            }
        };
        match &mapping_program {
            MappingProgram::BWA_MEM => {
                cmd
                    .arg("index")
                    .arg("-p")
                    .arg(&index_path)
                    .arg(&reference_path);
            },
            MappingProgram::MINIMAP2_SR |
            MappingProgram::MINIMAP2_ONT |
            MappingProgram::MINIMAP2_PB |
            MappingProgram::MINIMAP2_NO_PRESET => {
                match &mapping_program {
                    MappingProgram::MINIMAP2_SR =>  { cmd.arg("-x").arg("sr"); }
                    MappingProgram::MINIMAP2_ONT => { cmd.arg("-x").arg("map-ont"); }
                    MappingProgram::MINIMAP2_PB =>  { cmd.arg("-x").arg("map-pb"); }
                    MappingProgram::MINIMAP2_NO_PRESET |
                    MappingProgram::BWA_MEM => { }
                };
                cmd
                    .arg("-d")
                    .arg(&index_path)
                    .arg(&reference_path);
            }
        };
        match index_creation_options {
            Some(params) => {
                for s in params.split_whitespace() {
                    cmd.arg(s);
                }
            },
            None => {}
        };
        // Some BWA versions output log info to stdout. Ignore this.
        cmd.stdout(std::process::Stdio::piped());
        cmd.stderr(std::process::Stdio::piped());
        debug!("Running DB indexing command: {:?}", cmd);

        let mut process = cmd.spawn().expect(
            &format!("Failed to start {:?} index process", mapping_program));
        let es = process.wait().expect(
            &format!("Failed to glean exitstatus from failing {:?} index process", mapping_program));
        if !es.success() {
            error!("Error when running {:?} index process.", mapping_program);
            let mut err = String::new();
            process.stderr.expect(&format!(
                "Failed to grab stderr from failed {:?} index process", 
                mapping_program))
                .read_to_string(&mut err).expect("Failed to read stderr into string");
            error!("The STDERR was: {:?}", err);
            error!("Cannot continue after {:?} index failed.", mapping_program);
            process::exit(1);
        }
        info!("Finished generating {:?} index.", mapping_program);
        return TemporaryIndexStruct {
            index_path_internal: index_path.to_string_lossy().to_string(),
            tempdir: td
        }
    }
}
impl MappingIndex for TemporaryIndexStruct {
    fn index_path(&self) -> &String {
        return &self.index_path_internal
    }
}
impl Drop for TemporaryIndexStruct {
    fn drop(&mut self) {
        debug!("Dropping index tempdir ..")
    }
}

pub fn generate_bwa_index(
    reference_path: &str,
    index_creation_parameters: Option<&str>) -> Box<dyn MappingIndex> {
    let bwa_extensions = vec!("amb","ann","bwt","pac","sa");
    let num_extensions = bwa_extensions.len();
    let mut num_existing: u8 = 0;
    for extension in bwa_extensions {
        if std::path::Path::new(&format!("{}.{}", reference_path, extension)).exists() {
            num_existing += 1;
        }
    }
    if num_existing == 0 {
        return Box::new(TemporaryIndexStruct::new(
            MappingProgram::BWA_MEM, reference_path, index_creation_parameters));
    } else if num_existing as usize != num_extensions {
        error!("BWA index appears to be incomplete, cannot continue.");
        process::exit(1);
    } else {
        info!("BWA index appears to be complete, so going ahead and using it.");
        return Box::new(VanillaBwaIndexStuct::new(
            reference_path));
    }
}

pub fn generate_minimap2_index(
    reference_path: &str,
    index_creation_parameters: Option<&str>,
    mapping_program: MappingProgram)
    -> Box<dyn MappingIndex> {

    return Box::new(TemporaryIndexStruct::new(
            mapping_program,
            reference_path,
            index_creation_parameters));
}

pub fn generate_concatenated_fasta_file(
    fasta_file_paths: &Vec<String>)
    -> NamedTempFile {

    let tmpfile: NamedTempFile = NamedTempFile::new().unwrap();
    let mut something_written_at_all = false;
    { // scope so writer, which borrows tmpfile, goes out of scope.
        let mut writer = bio::io::fasta::Writer::new(&tmpfile);
        let mut genome_names: HashSet<String> = HashSet::new();

        for file in fasta_file_paths {
            let mut something_written = false;
            let path = std::path::Path::new(file);
            let reader = bio::io::fasta::Reader::from_file(path)
                .expect(&format!("Unable to read fasta file {}", file));

            let genome_name = String::from(
                path.file_stem()
                    .expect("Problem while determining file stem")
                    .to_str()
                    .expect("File name string conversion problem"));
            if genome_names.contains(&genome_name) {
                error!("The genome name {} was derived from >1 file", genome_name);
                process::exit(1);
            }
            for record in reader.records() {
                something_written = true;
                something_written_at_all = true;
                let r = record.unwrap();
                writer.write(
                    &format!("{}{}{}", genome_name, CONCATENATED_FASTA_FILE_SEPARATOR, r.id()),
                    r.desc(),
                    r.seq()
                ).unwrap()
            }
            genome_names.insert(genome_name);
            if !something_written {
                error!(
                    "FASTA file {} appears to be empty as no sequences were contained in it",
                    file);
                process::exit(1);
            }
        }
    }
    if !something_written_at_all {
        error!("Concatenated FASTA file to use as a reference is empty");
        process::exit(1);
    }
    return tmpfile;
}
