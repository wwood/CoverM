// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Utility methods.
use std::collections::HashMap;
use std::fmt::Debug;
use std::fs::{File};
use std::io::{self, Write, BufRead, BufReader, BufWriter};
use std::path::Path;
use std::sync::{Arc, Mutex};

use bincode::{self, deserialize_from, serialize_into};
use failure::Error;
use flate2::read::MultiGzDecoder;
use serde::{Serialize, de::DeserializeOwned};

use bio::io::{fasta, fastq};
use debruijn::dna_string::DnaString;

use pseudoaligner::config;

const MAPPABILITY_HEADER_STRING: &'static str = "tx_name\tgene_name\ttx_kmer_count\ttx_fraction_unique\tgene_fraction_unique\n";

pub fn write_obj<T: Serialize, P: AsRef<Path> + Debug>(
    g: &T,
    filename: P,
) -> Result<(), bincode::Error> {
    let f = match File::create(&filename) {
        Err(err) => panic!("couldn't create file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut writer = BufWriter::new(f);
    serialize_into(&mut writer, &g)
}

pub fn read_obj<T: DeserializeOwned, P: AsRef<Path> + Debug>(
    filename: P,
) -> Result<T, bincode::Error> {
    let f = match File::open(&filename) {
        Err(err) => panic!("couldn't open file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut reader = BufReader::new(f);
    deserialize_from(&mut reader)
}

/// Open a (possibly gzipped) file into a BufReader.
fn _open_with_gz<P: AsRef<Path>>(p: P) -> Result<Box<BufRead>, Error> {
    let r = File::open(p.as_ref())?;

    if p.as_ref().extension().unwrap() == "gz" {
        let gz = MultiGzDecoder::new(r);
        let buf_reader = BufReader::with_capacity(32 * 1024, gz);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::with_capacity(32 * 1024, r);
        Ok(Box::new(buf_reader))
    }
}

pub fn read_transcripts(
    reader: fasta::Reader<File>,
) -> Result<(Vec<DnaString>, Vec<String>, HashMap<String, String>), Error> {
    let mut seqs = Vec::new();
    let mut transcript_counter = 0;
    let mut tx_ids = Vec::new();
    let mut tx_to_gene_map = HashMap::new();

    let mut fasta_format: Option<u8> = None;

    info!("Starting reading the Fasta file\n");
    for result in reader.records() {
        // obtain record or fail with error
        let record = result?;

        // Sequence
        let dna_string = DnaString::from_acgt_bytes_hashn(record.seq(), record.id().as_bytes());
        seqs.push(dna_string);

        if let None = fasta_format {
            fasta_format = detect_fasta_format(&record);
        }

        let (tx_id, gene_id) = extract_tx_gene_id(&record, fasta_format)?;

        tx_ids.push(tx_id.clone());
        tx_to_gene_map.insert(tx_id, gene_id);

        transcript_counter += 1;
        if transcript_counter % 100 == 0 {
            print!("\r Done reading {} sequences", transcript_counter);
            io::stdout().flush().expect("Could not flush stdout");
        }
    }

    println!();
    info!(
        "Done reading the Fasta file; Found {} sequences",
        transcript_counter
    );

    Ok((seqs, tx_ids, tx_to_gene_map))
}

pub fn detect_fasta_format(record: &fasta::Record) -> Option<u8> {
    let id_tokens: Vec<&str> = record.id().split('|').collect();
    if id_tokens.len() == 9 {
        return Some(config::FASTA_FORMAT_GENCODE)
    }
    match record.desc() {
        Some(desc) => {
            let desc_tokens: Vec<&str> = desc.split(' ').collect();
            if desc_tokens.len() == 5 {
                Some(config::FASTA_FORMAT_ENSEMBL)
            } else {
                None
            }
        },
        None => None
    }
}

pub fn extract_tx_gene_id(record: &fasta::Record, fasta_format: Option<u8>) -> Result<(String, String), Error>{
    match fasta_format {
        Some(config::FASTA_FORMAT_GENCODE) => {
            let id_tokens: Vec<&str> = record.id().split('|').collect();
            let tx_id = id_tokens[0].to_string();
            let gene_id = id_tokens[1].to_string();
            Ok((tx_id, gene_id))
        },
        Some(config::FASTA_FORMAT_ENSEMBL) => {
            let tx_id = record.id().to_string();
            let desc_tokens: Vec<&str> = record.desc().unwrap().split(' ').collect();
            let gene_tmp: Vec<&str> = desc_tokens[2].split(':').collect();
            let gene_id = gene_tmp[1].to_string();
            Ok((tx_id, gene_id))
        },
        _ => {
            Ok((record.id().to_string(), record.id().to_string()))
        }
    }
}

pub fn get_next_record<R: io::Read>(
    reader: &Arc<Mutex<fastq::Records<R>>>,
) -> Option<Result<fastq::Record, io::Error>> {
    let mut lock = reader.lock().unwrap();
    lock.next()
}

pub fn open_file<P: AsRef<Path>>(
    filename: &str, outdir: P
) -> Result<File, Error> {
    let out_fn = outdir.as_ref().join(filename);
    let outfile = File::create(&out_fn)?;
    Ok(outfile)
}
