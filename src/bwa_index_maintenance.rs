use std;
use std::io::Read;

use tempdir::TempDir;


pub trait BwaIndexStruct {
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
impl BwaIndexStruct for VanillaBwaIndexStuct {
    fn index_path(&self) -> &String {
        return &self.index_path_internal
    }
}

pub struct TemporaryBwaIndexStruct {
    #[allow(dead_code)] // field is never used, it just needs to be kept in scope.
    tempdir: TempDir,
    index_path_internal: String
}

impl TemporaryBwaIndexStruct {
    pub fn new(reference_path: &str) -> TemporaryBwaIndexStruct {
        // Generate a BWA index in a temporary directory, where the temporary
        // directory does not go out of scope until the struct does.
        let td = TempDir::new("coverm-bwa-index")
            .expect("Unable to create temporary directory");
        let index_path = std::path::Path::new(td.path())
            .join(std::path::Path::new(reference_path).file_name()
                  .expect("Failed to glean file stem from reference DB. Strange."));

        info!("Generating BWA index for {} ..", reference_path);
        let mut cmd = std::process::Command::new("bwa");
        cmd
            .arg("index")
            .arg("-p")
            .arg(&index_path)
            .arg(&reference_path)
            .stderr(std::process::Stdio::piped());
        let mut process = cmd.spawn().expect("Failed to start BWA index process");
        let es = process.wait().expect("Failed to glean exitstatus from failing BWA index process");
        if !es.success() {
            error!("Error when running bwa index process.");
            let mut err = String::new();
            process.stderr.expect("Failed to grab stderr from failed BWA index process")
                .read_to_string(&mut err).expect("Failed to read stderr into string");
            error!("The STDERR was: {:?}", err);
            panic!("Cannot continue after BWA index failed.");
        }
        info!("Finished generating BWA index.");
        return TemporaryBwaIndexStruct {
            index_path_internal: index_path.to_string_lossy().to_string(),
            tempdir: td
        }
    }
}
impl BwaIndexStruct for TemporaryBwaIndexStruct {
    fn index_path(&self) -> &String {
        return &self.index_path_internal
    }
}
impl Drop for TemporaryBwaIndexStruct {
    fn drop(&mut self) {
        debug!("Dropping index tempdir ..")
    }
}

pub fn generate_bwa_index(reference_path: &str) -> Box<dyn BwaIndexStruct> {
    let bwa_extensions = vec!("amb","ann","bwt","pac","sa");
    let num_extensions = bwa_extensions.len();
    let mut num_existing: u8 = 0;
    for extension in bwa_extensions {
        if std::path::Path::new(&format!("{}.{}", reference_path, extension)).exists() {
            num_existing += 1;
        }
    }
    if num_existing == 0 {
        return Box::new(TemporaryBwaIndexStruct::new(reference_path));
    } else if num_existing as usize != num_extensions {
        panic!("BWA index appears to be incomplete, cannot continue.");
    } else {
        info!("BWA index appears to be complete, so going ahead and using it.");
        return Box::new(VanillaBwaIndexStuct::new(reference_path));
    }
}
