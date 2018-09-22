use std;
use std::io::Read;

pub trait BwaIndexStruct<'a> {
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
impl<'a> BwaIndexStruct<'a> for VanillaBwaIndexStuct {
    fn index_path(&self) -> &String {
        return &self.index_path_internal
    }
}

pub fn generate_bwa_index(reference_path: &str) -> VanillaBwaIndexStuct {
    let bwa_extensions = vec!("amb","ann","bwt","pac","sa");
    let num_extensions = bwa_extensions.len();
    let mut num_existing: u8 = 0;
    for extension in bwa_extensions {
        if std::path::Path::new(&format!("{}.{}", reference_path, extension)).exists() {
            num_existing += 1;
        }
    }
    if num_existing == 0 {
        info!("Generating BWA index for {} ..", reference_path);
        let mut cmd = std::process::Command::new("bwa");
        cmd
            .arg("index")
            .arg(reference_path)
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
        return VanillaBwaIndexStuct::new(reference_path);
    } else if num_existing as usize != num_extensions {
        panic!("BWA index appears to be incomplete, cannot continue.");
    } else {
        info!("BWA index appears to be complete, so going ahead and using it.");
        return VanillaBwaIndexStuct::new(reference_path);
    }
}
