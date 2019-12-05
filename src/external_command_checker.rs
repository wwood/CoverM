use std;
use std::io::Read;
use std::process;
use version_compare::Version;

pub fn check_for_bwa() {
    self::check_for_external_command_presence("BWA", "which bwa");
}

pub fn check_for_minimap2() {
    self::check_for_external_command_presence("minimap2", "which minimap2");
    self::default_version_check("minimap2","2.17-r941");
}

pub fn check_for_samtools() {
    self::check_for_external_command_presence("samtools", "which samtools");
    self::default_version_check("samtools","1.9");
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

fn default_version_check(executable_name: &str, min_version: &str) {
    let version_command = format!("{} --version", executable_name);
    let mut cmd = std::process::Command::new("bash");
    cmd
        .arg("-c")
        .arg(&version_command)
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
        error!("Cannot continue without {}. Finding version of `{}` failed",
               executable_name, &version_command);
        process::exit(1);
    }
    let mut version = String::new();
    process.stdout.expect("Failed to grab stdout from failed command version finding process")
        .read_to_string(&mut version).expect("Failed to read stdout into string");
    version = version.trim().to_string();

    let expected_version = Version::from(min_version)
        .expect("Programming error: failed to parse code-specified version");
    let found_version = Version::from(
        version
            .lines()
            .next()
            .expect(&format!("Unable to parse version for {} (error 1)", &executable_name))
            .trim()
            .rsplit(' ')
            .next()
            .expect(&format!("Unable to parse version for {} (error 2)", &executable_name)))

        .expect(&format!("Unable to parse version number '{}' from executable {}",
            version, executable_name));

    info!("Found {} version {} ", executable_name, found_version);
    if found_version < expected_version {
        error!("It appears the available version of {} is too old \
            (found version {}, required is {})",
            executable_name, found_version, expected_version);
        process::exit(11);
    }
}
