use std;

use csv;

use mapping_parameters::*;

pub fn screen_and_print_matching_genomes(
    mash_db_path_str: &str,
    read_inputs: &MappingParameters,
    num_threads: usize,
    min_identity: f32) {

    info!("Running mash screen against db {} ..", mash_db_path_str);

    for readset in read_inputs.readsets() {
        let mut cmd = std::process::Command::new("mash");
        cmd
            .arg("screen")
            .arg("-p")
            .arg(format!("{}", num_threads))
            .arg("-i")
            .arg(&format!("{}",min_identity))
            .arg(mash_db_path_str)
            .arg(readset.0) // Ignore r2 for now. TODO: Fix?
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped());
        info!("Running mash screen process on {}", readset.0);
        let process = cmd.spawn().expect("Failed to start mash screen process");
        // process.wait().expect("Failed to finish process");
        // let mut s = String::new();
        // process.stdout.expect("f").read_to_string(&mut s).expect("g");
        // debug!("stdout: {}", s);


        let mut child = process.stdout.unwrap();
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(&mut child);
        for result in rdr.records() {
            match result {
                Err(e) => {
                    error!("Failed to read CSV output from mash screen output. The error was {}", e);
                    std::process::exit(1);
                }
                Ok(res) => {
                    // The output fields are
                    //     [identity, shared-hashes, median-multiplicity, p-value, query-ID,
                    //      query-comment], where median-multiplicity is computed for shared hashes, based
                    //     on the number of observations of those hashes within the pool.
                    debug!("Got output from mash: {:?}", res);
                    if res.len() != 6 {
                        error!("Unexpected mash screen output line encountered: {:?}", res);
                        std::process::exit(1);
                    }
                    println!("{}\t{}", &res[4], &res[0]);
                }
            }
        }

        // TODO: Parse exitstatus etc. Having borrow checker problems.
        // let es = process.wait().expect("Failed to glean exitstatus from failing mash screen process");
        // if !es.success() {
        //     error!("Error when running mash screen process.");
        //     let mut err = String::new();
        //     process.stderr.expect("Failed to grab stderr from failed mash screen process")
        //         .read_to_string(&mut err).expect("Failed to read stderr into string");
        //     error!("The STDERR was: {:?}", err);
        //     error!("Cannot continue after mash screen failed.");
        //     process::exit(1);
        // }
        info!("Finished mash screen subcommand.");
    }
}
