use std;
use std::io::Read;
use std::path::Path;

use core_genome::CoreGenomicRegion;
use csv;
use tempdir::TempDir;
use bio::io::fasta;

use run_command_safely;

pub fn parsnp_core_genomes_from_genome_fasta_files(
    genome_fasta_paths: &[&str],
    clade_id: u32
) -> Vec<Vec<CoreGenomicRegion>> {
    // Get absolute paths of each genome fasta file
    let absolute_fasta_paths: Vec<std::path::PathBuf> = genome_fasta_paths.iter().map(|f| {
        std::path::Path::new(f).canonicalize()
            .expect(&format!("Failed to canonicalised input fasta path {}", f))
    }).collect();

    // Make a tempdir
    let tmp_dir = TempDir::new("coverm-parsnp")
        .expect("Failed to make temporary directory for parsnp");

    // In the tempdir, make a genomes directory
    let genomes_path = tmp_dir.path().join("genomes");
    std::fs::create_dir(&genomes_path)
        .expect("Failed to create genomes folder in parsnp tempdir");

    // Change working directory to genome
    let original_working_directory = std::env::current_dir()
        .expect("Failed to get working directory");
    std::env::set_current_dir(&genomes_path)
        .expect("Failed to run set working directory");

    // Write each genome as a new file
    debug!("Writing temporary FASTA files for parsnp ..");
    for (i, abs_path) in absolute_fasta_paths.iter().enumerate() {
        let output_path = genomes_path.join(std::path::Path::new(
            &format!("genome{}.fasta", i+1)));
        let reader = fasta::Reader::from_file(abs_path)
            .expect("Failed to open read fasta file for parsnp");
        let mut writer = fasta::Writer::to_file(output_path)
            .expect("Failed to open write fasta file for parsnp");
        for (contig_idx, record_res) in reader.records().enumerate() {
            let record = record_res.expect("Failed to parse FASTA");

            writer.write(
                &format!("{} contig{}",
                         i, contig_idx),
                None,
                &record.seq())
                .expect("Failed to write FASTA for parsnp input");
        }
        writer.flush().expect("Failed to flush parsnp temp genome writer");
    }
    debug!("Finished writing temporary FASTA files for parsnp");

    // Change directory to genome/..
    std::env::set_current_dir(tmp_dir.path())
        .expect("Failed to run set working directory");

    // Run parsnp using the first genome as the reference
    let mut cmd = std::process::Command::new("parsnp");
    cmd
        .arg("-r")
        // Use the second genome as the reference one, otherwise it seems not
        // detect them. I dunno.
        .arg(tmp_dir.path().join(Path::new("genomes")).join(Path::new(&format!("genome2.fasta"))))
        .arg("-c")
        .arg("genomes")
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    debug!("Running parsnp cmd: {:?}", &cmd);
    // Parsnp can fail at steps after MUMi part, and that's OK.
    let mut process = cmd.spawn().expect(&format!("Failed to spawn {}", "parsnp"));
    process.wait()
        .expect(&format!("Failed to glean exitstatus from failing {} process", "parsnp"));
    // if the stdout has certain words, we are in business.
    let mut stdout = String::new();
    process.stdout.unwrap().read_to_string(&mut stdout).expect("Failed to read parsnp stdout");
    if stdout.find("Running PhiPack").is_none() {
        error!("Parsnp failed.");
        let mut err = String::new();
        process.stderr.expect(&format!("Failed to grab stderr from failed {} process", "parsnp"))
            .read_to_string(&mut err).expect("Failed to read stderr into string");
        error!("The STDERR was: {:?}", err);
        error!("The STDOUT was: {:?}", stdout);
    }

    // TODO: Comment out this debug code.
    let mut cp = std::process::Command::new("cp");
    cp.arg("-r")
        .arg(tmp_dir.path())
        .arg("/tmp/tester").spawn().unwrap().wait().unwrap();

    // Find the P_ directory in the current directory. parsnp -o does not appear
    // to work as intended.
    let mut p_dir = None;
    for entry in std::fs::read_dir(tmp_dir.path())
        .expect("Failed to iterate parsnp tmpdir") {
            let e2 = entry
                .expect("Failed to readdir parsnp tmpdir entry");
            let e3 = e2.file_name();
            let name = e3
                .to_str()
                .expect("Failed to convert parsnp tmpdir entry to str");
            if name.get(0..1) == Some(&"P") {
                debug!("Found P directory {:?}", name);
                p_dir = match p_dir {
                    None => Some(name.to_string()),
                    Some(x) => panic!("Found 2 or more P dirs: {} and {}", x, name)
                }
            }
        }

    // Run harvesttools to get the backbone intervals
    let mut harvest = std::process::Command::new("harvesttools");
    harvest
        .arg("-i")
        .arg(
            Path::new(&p_dir.expect("Found no P directory output from parsnp"))
                .join("parsnp.ggr"))
        .arg("-B")
        .arg("-")
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    let harvest_child = run_command_safely(&mut harvest, "harvest");

    let cores = harvest_child.stdout.unwrap();

    std::env::set_current_dir(original_working_directory)
        .expect("Failed to run set working directory");

    // Parse backbone file
    let named_regionset = core_genomes_from_backbone_intervals(
        cores, clade_id);

    // Reorder the set as the parsnp output order I don't think is
    // deterministic.
    let indices: Vec<usize> = named_regionset.names.iter().map(|name| {
        name.split(" ").next().unwrap().parse::<usize>().unwrap()
    }).collect();
    debug!("New indices for core-genomes: {:?}", indices);

    let mut new_vec = vec![vec![]; indices.len()];
    for (prev_i, new_i) in indices.iter().enumerate() {
        // TODO: Probably better ways than clone here
        new_vec[*new_i] = named_regionset.regions[prev_i].clone();
    }

    return new_vec;
}

#[derive(Debug)]
struct NamedCoreGenomicRegionSet {
    names: Vec<String>,
    regions: Vec<Vec<CoreGenomicRegion>>,
}

fn core_genomes_from_backbone_intervals<R: Read>(
    reader: R,
    clade_id: u32)
    -> NamedCoreGenomicRegionSet {

    // example files:
    // >genome1 random_start	>genome1 random_end	>genome2 one A to T SNP at position 100_start	>genome2 one A to T SNP at position 100_end	>genome1_near_duplicate genome1_plus_A_at_start_start	>genome1_near_duplicate genome1_plus_A_at_start_end
    // 1	200	1	200	2	201
    //
    // NOTE: Parsnp treats the genome fasta files as one continuous sequence,
    // ignoring lines starting with '>' except for the first line of the file.
    // This means we have to map the entire alignment positions to contig IDs.
    // Alternately, as a shortcut we can just do the same here, treating each
    // genome as a single contig.
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);

    // Read one line, which will be a single LCB
    let headers = rdr.headers()
        .expect("Failed to get header line from backbone interval stream");
    let num_headers = headers.len();
    assert!(num_headers > 0);
    assert!(num_headers % 2 == 0);
    let num_genomes = num_headers / 2;

    // Extract the names of each of the genomes
    let mut names = vec![];
    let mut header_i = 0;
    while header_i < num_headers {
        // Name is the header minus the _start which is 6 chars long
        let mut name = headers[header_i].to_string();
        name.replace_range(name.len()-6.., "");
        name.replace_range(0..1,""); // Remove the '>'.
        names.push(name);
        header_i += 2;
    }
    debug!("Got names from backbone: {:?}", names);

    let mut core_regions: Vec<Vec<CoreGenomicRegion>> = vec![vec![]; num_genomes];

    for record_res in rdr.records() {
        let mut i = 0;
        match record_res {
            Ok(record) => {
                debug!("Parsing backbone interval line: {:?}", record);
                while i/2 < num_genomes {
                    let genome_idx = i/2;
                    core_regions[genome_idx].push(
                        CoreGenomicRegion {
                            clade_id: clade_id,
                            contig_id: 0,
                            start: record[i].parse()
                                .expect("Failed to parse LCB start number from backbone intervals"),
                            stop: record[i+1].parse()
                                .expect("Failed to parse LCB stop number from backbone intervals"),
                        }
                    );
                    i += 2;
                }
            },
            Err(e) => {
                error!("Error parsing backbone intervals: {}", e);
                std::process::exit(1);
            }
        }
    }
    debug!("Got original regions from backbone: {:?}", core_regions);

    return NamedCoreGenomicRegionSet {
        names: names,
        regions: core_regions
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn test_backbone_parser_hello_world() {
        init();
        let backbones =
            ">genome1_near_duplicate genome1_plus_A_at_start_start	>genome1_near_duplicate genome1_plus_A_at_start_end	>genome1_2_start	>genome1_2_end	>genome2 one A to T SNP at position 100_start	>genome2 one A to T SNP at position 100_end\n\
             2	88	114	200	1	87\n\
             89	201	1	113	88	200\n".as_bytes();
        assert_eq!(
            vec![
                vec![
                    CoreGenomicRegion {
                        clade_id: 7,
                        contig_id: 0,
                        start: 2,
                        stop: 88,
                    },
                    CoreGenomicRegion {
                        clade_id: 7,
                        contig_id: 0,
                        start: 89,
                        stop: 201,
                    },
                ],
                vec![
                    CoreGenomicRegion {
                        clade_id: 7,
                        contig_id: 0,
                        start: 114,
                        stop: 200,
                    },
                    CoreGenomicRegion {
                        clade_id: 7,
                        contig_id: 0,
                        start: 1,
                        stop: 113,
                    },
                ],
                vec![
                    CoreGenomicRegion {
                        clade_id: 7,
                        contig_id: 0,
                        start: 1,
                        stop: 87,
                    },
                    CoreGenomicRegion {
                        clade_id: 7,
                        contig_id: 0,
                        start: 88,
                        stop: 200,
                    },
                ]
            ],
            core_genomes_from_backbone_intervals(backbones, 7).regions);
    }

    #[test]
    fn test_parsnp_hello_world() {
        init();
        let cores = parsnp_core_genomes_from_genome_fasta_files(
            &vec![
                "tests/data/2_single_species_dummy_dataset/2genomes/parsnp/g2.fasta",
                "tests/data/2_single_species_dummy_dataset/2genomes/parsnp/g1_split.fasta",
                "tests/data/2_single_species_dummy_dataset/2genomes/parsnp/g1_duplicate.fasta",
            ],
            7
        );
        assert_eq!(
            vec![

                vec![
                    CoreGenomicRegion {
                        clade_id: 7,
                        contig_id: 0,
                        start: 88,
                        stop: 200,
                    },
                    CoreGenomicRegion {
                        clade_id: 7,
                        contig_id: 0,
                        start: 1,
                        stop: 87,
                    },
                ],

                vec![
                    CoreGenomicRegion {
                        clade_id: 7,
                        contig_id: 0,
                        start: 1,
                        stop: 113,
                    },
                    CoreGenomicRegion {
                        clade_id: 7,
                        contig_id: 0,
                        start: 114,
                        stop: 200,
                    },
                ],

                vec![
                    CoreGenomicRegion {
                        clade_id: 7,
                        contig_id: 0,
                        start: 89,
                        stop: 201,
                    },
                    CoreGenomicRegion {
                        clade_id: 7,
                        contig_id: 0,
                        start: 2,
                        stop: 88,
                    },
                ],

            ],
            cores);
    }
}
