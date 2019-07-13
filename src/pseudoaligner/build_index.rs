// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

use std::sync::Arc;
use std::collections::HashMap;

use boomphf::hashmap::{BoomHashMap2, NoKeyBoomHashMap};
use pseudoaligner::config::{KmerType, MEM_SIZE, REPORT_ALL_KMER, STRANDED};
use debruijn;
use debruijn::compression::*;
use debruijn::dna_string::{DnaString, DnaStringSlice};
use debruijn::filter::*;
use debruijn::graph::*;
use debruijn::*;

use boomphf;
use failure::Error;
use pseudoaligner::config::{MAX_WORKER, MIN_KMERS, U32_MAX};
use pseudoaligner::pseudoaligner::Pseudoaligner;
use rayon;
use rayon::prelude::*;

const MIN_SHARD_SEQUENCES: usize = 2000;

pub fn build_index<K: Kmer + Sync + Send>(
    seqs: &[DnaString],
    tx_names: &Vec<String>,
    tx_gene_map: &HashMap<String, String>
) -> Result<Pseudoaligner<K>, Error> {
    // Thread pool Configuration for calling BOOMphf
    rayon::ThreadPoolBuilder::new()
        .num_threads(MAX_WORKER)
        .build_global()?;

    if seqs.len() >= U32_MAX {
        panic!("Too many ({}) sequences to handle.", seqs.len());
    }

    println!("Sharding sequences...");

    let mut buckets: Vec<_> = seqs
        .into_par_iter()
        .enumerate()
        .flat_map(|(id, seq)| partition_contigs::<KmerType>(seq, id as u32))
        .collect();

    buckets.par_sort_unstable_by_key(|x| x.0);
    println!("Got {} sequence chunks", buckets.len());

    let summarizer = Arc::new(debruijn::filter::CountFilterEqClass::new(MIN_KMERS));
    let sequence_shards = group_by_slices(&buckets, |x| x.0, MIN_SHARD_SEQUENCES);

    let mut shard_dbgs = Vec::with_capacity(sequence_shards.len());

    println!("Assembling {} shards...", sequence_shards.len());

    sequence_shards
        .into_par_iter()
        .map_with(summarizer.clone(), |s, strings| {
            assemble_shard::<K>(strings, s)
        }).collect_into_vec(&mut shard_dbgs);

    println!();
    println!("Done separate de Bruijn graph construction");
    println!("Starting merging disjoint graphs");

    //println!("{:?}", summarizer);
    let dbg = merge_shard_dbgs(shard_dbgs);
    println!("Merger of graphs complete");

    // TODO update rust-debruijn version and fix this
    let eq_classes = summarizer.get_eq_classes();

    println!("Indexing de Bruijn graph");
    let dbg_index = make_dbg_index(&dbg);
    Ok(Pseudoaligner::new(
        dbg, eq_classes, dbg_index, tx_names.clone(), tx_gene_map.clone()
    ))
}

type PmerType = debruijn::kmer::Kmer6;

lazy_static! {
    static ref PERM: Vec<usize> = {
        let maxp = 1 << (2 * PmerType::k());
        let mut permutation = Vec::with_capacity(maxp);
        for i in 0..maxp {
            permutation.push(i);
        }
        permutation
    };
}

fn partition_contigs<'a, K: Kmer>(
    contig: &'a DnaString,
    contig_id: u32,
) -> Vec<(u16, u32, DnaStringSlice<'a>, Exts)> {
    // One FASTA entry possibly broken into multiple contigs
    // based on the location of `N` int he sequence.

    let mut bucket_slices = Vec::new();

    if contig.len() >= K::k() {
        let msps = debruijn::msp::simple_scan::<_, PmerType>(K::k(), contig, &PERM, false);
        for msp in msps {
            let bucket_id = msp.bucket();
            let slice = contig.slice(msp.start(), msp.end());
            let exts = Exts::from_dna_string(contig, msp.start(), msp.len());
            bucket_slices.push((bucket_id, contig_id, slice, exts));
        }
    }

    bucket_slices
}

fn assemble_shard<K: Kmer>(
    shard_data: &[(u16, u32, DnaStringSlice, Exts)],
    summarizer: &Arc<CountFilterEqClass<u32>>,
) -> BaseGraph<K, EqClassIdType> {
    let filter_input: Vec<_> = shard_data
        .into_iter()
        .cloned()
        .map(|(_, seqid, string, exts)| (string, exts, seqid))
        .collect();

    let (phf, _): (BoomHashMap2<K, Exts, EqClassIdType>, _) = filter_kmers(
        &filter_input,
        summarizer,
        STRANDED,
        REPORT_ALL_KMER,
        MEM_SIZE,
    );

    compress_kmers_with_hash(STRANDED, &ScmapCompress::new(), &phf)
}

fn merge_shard_dbgs<K: Kmer + Sync + Send>(
    uncompressed_dbgs: Vec<BaseGraph<K, EqClassIdType>>,
) -> DebruijnGraph<K, EqClassIdType> {
    let combined_graph = BaseGraph::combine(uncompressed_dbgs.into_iter()).finish();
    compress_graph(STRANDED, &ScmapCompress::new(), combined_graph, None)
}

#[inline(never)]
fn make_dbg_index<K: Kmer + Sync + Send>(
    dbg: &DebruijnGraph<K, EqClassIdType>,
) -> NoKeyBoomHashMap<K, (u32, u32)> {
    let mut total_kmers = 0;
    let kmer_length = K::k();
    for node in dbg.iter_nodes() {
        total_kmers += node.len() - kmer_length + 1;
    }

    println!("Total {:?} kmers to process in dbg", total_kmers);
    println!("Making mphf of kmers");
    let mphf = boomphf::Mphf::new_parallel_with_keys(1.7, dbg, None, total_kmers, MAX_WORKER);

    println!("Assigning offsets to kmers");
    let mut node_and_offsets = Vec::with_capacity(total_kmers);
    node_and_offsets.resize(total_kmers, (U32_MAX as u32, U32_MAX as u32));

    for node in dbg {
        let node_id = node.node_id;

        for (offset, kmer) in node.into_iter().enumerate() {
            let index = match mphf.try_hash(&kmer) {
                None => panic!("can't find kmer"),
                Some(index) => index,
            };

            node_and_offsets[index as usize] = (node_id as u32, offset as u32);
        }
    }

    boomphf::hashmap::NoKeyBoomHashMap::new_with_mphf(mphf, node_and_offsets)
}

fn group_by_slices<T, K: PartialEq, F: Fn(&T) -> K>(
    data: &[T],
    f: F,
    min_size: usize,
) -> Vec<&[T]> {
    let mut slice_start = 0;
    let mut result = Vec::new();
    for i in 1..data.len() {
        if !(f(&data[i - 1]) == f(&data[i])) && (i - slice_start) > min_size {
            result.push(&data[slice_start..i]);
            slice_start = i;
        }
    }
    if slice_start > 0 {
        result.push(&data[slice_start..]);
    }
    result
}
