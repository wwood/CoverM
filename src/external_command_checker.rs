use std;
use std::io::Read;
use std::process;

pub fn check_for_bwa() {
    self::check_for_external_command_presence("BWA", "which bwa");
}

pub fn check_for_samtools() {
    self::check_for_external_command_presence("samtools", "which samtools");
}

fn check_for_external_command_presence(
    executable_name: &str, testing_cmd: &str) {
    debug!("Checking for {} ..", executable_name);
    let mut cmd = std::process::Command::new("bash");
    cmd
        .arg("-c")
        .arg(testing_cmd)
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    let mut process = cmd.spawn().expect("Unable to execute bash");
    let es = process.wait()
        .expect(&format!(
            "Failed to glean exitstatus while checking for presence of {}",
            executable_name));
    if !es.success() {
        error!("Could not find an available {} executable.", executable_name);
        let mut err = String::new();
        process.stderr.expect("Failed to grab stderr from failed executable finding process")
            .read_to_string(&mut err).expect("Failed to read stderr into string");
        error!("The STDERR was: {:?}", err);
        error!("Cannot continue without {}. Testing for presence with `{}` failed",
               executable_name, testing_cmd);
        process::exit(1);
    }
}
