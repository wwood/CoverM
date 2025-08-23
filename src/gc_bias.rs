use bio::io::fasta;
use plotters::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use splines::{Interpolation, Key, Spline};
use std::collections::HashMap;
use std::error::Error;

/// Fit a spline to GC fraction and relative coverage pairs.
/// Returns the spline and the binned points used for fitting.
pub(crate) fn fit_gc_bias_spline(
    gc: &[f64],
    rel_cov: &[f64],
) -> (Spline<f64, f64>, Vec<(f64, f64)>) {
    assert_eq!(gc.len(), rel_cov.len());
    let bins = 20usize;
    let mut bin_sums = vec![0.0f64; bins];
    let mut bin_counts = vec![0usize; bins];
    for (&g, &r) in gc.iter().zip(rel_cov.iter()) {
        let g_clamped = g.clamp(0.0, 1.0);
        let idx = ((g_clamped * (bins as f64 - 1.0)).round()) as usize;
        bin_sums[idx] += r;
        bin_counts[idx] += 1;
    }
    let mut points: Vec<(f64, f64)> = Vec::new();
    let mut keys: Vec<Key<f64, f64>> = Vec::new();
    for i in 0..bins {
        if bin_counts[i] > 0 {
            let x = (i as f64 + 0.5) / bins as f64;
            let y = bin_sums[i] / bin_counts[i] as f64;
            points.push((x, y));
            keys.push(Key::new(x, y, Interpolation::Linear));
        }
    }
    (Spline::from_vec(keys.clone()), points)
}

fn plot_spline(
    points: &[(f64, f64)],
    spline: &Spline<f64, f64>,
    path: &str,
) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new(path, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let y_min = points.iter().map(|(_, y)| *y).fold(f64::INFINITY, f64::min);
    let y_max = points
        .iter()
        .map(|(_, y)| *y)
        .fold(f64::NEG_INFINITY, f64::max);
    let mut chart = ChartBuilder::on(&root)
        .caption("GC bias spline", ("sans-serif", 20))
        .margin(20)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(0f64..1f64, y_min..y_max)?;
    chart
        .configure_mesh()
        .x_desc("GC fraction")
        .y_desc("Relative coverage")
        .draw()?;
    chart.draw_series(
        points
            .iter()
            .map(|(x, y)| Circle::new((*x, *y), 3, BLUE.filled())),
    )?;
    chart.draw_series(LineSeries::new(
        (0..100).map(|i| {
            let gc = i as f64 / 99.0;
            let y = spline.clamped_sample(gc).unwrap_or(1.0);
            (gc, y)
        }),
        &RED,
    ))?;
    root.present()?;
    Ok(())
}

/// Calculate GC bias adjusted coverage for each contig given a reference fasta
/// and BAM files. Returns a vector of (contig_name, adjusted_coverage, original_coverage, gc_percent).
/// Only contigs with mean coverage >= `min_contig_cov` are used to fit the spline.
#[allow(clippy::type_complexity)]
pub fn gc_bias_correct(
    reference: &str,
    bam_files: &[&str],
    window_size: usize,
    min_contig_cov: f64,
    plot_path: Option<&str>,
) -> Result<Vec<(String, f64, f64, f64)>, Box<dyn Error>> {
    // Read reference sequences
    let fasta_reader = fasta::Reader::from_file(reference)?;
    let mut sequences: Vec<(String, Vec<u8>)> = Vec::new();
    for r in fasta_reader.records() {
        let rec = r?;
        sequences.push((rec.id().to_string(), rec.seq().to_vec()));
    }
    let mut name_to_index = HashMap::new();
    for (i, (name, _)) in sequences.iter().enumerate() {
        name_to_index.insert(name.clone(), i);
    }

    // Prepare coverage vectors
    let mut coverage: Vec<Vec<u32>> = sequences
        .iter()
        .map(|(_, seq)| vec![0u32; seq.len()])
        .collect();

    // Iterate through BAM files and accumulate coverage using pileup
    for bam_path in bam_files {
        let mut reader = bam::Reader::from_path(bam_path)?;
        let header = reader.header().to_owned();
        let mut tid_to_index: HashMap<i32, usize> = HashMap::new();
        for (tid, name_bytes) in header.target_names().iter().enumerate() {
            if let Ok(name) = std::str::from_utf8(name_bytes) {
                if let Some(idx) = name_to_index.get(name) {
                    tid_to_index.insert(tid as i32, *idx);
                }
            }
        }
        for p in reader.pileup() {
            let pile = p?;
            let tid = pile.tid();
            if let Some(&idx) = tid_to_index.get(&(tid as i32)) {
                let pos = pile.pos() as usize;
                if pos < coverage[idx].len() {
                    coverage[idx][pos] += pile.depth();
                }
            }
        }
    }

    struct ContigData {
        gc: Vec<f64>,
        cov: Vec<f64>,
    }
    let mut contig_data: Vec<ContigData> = Vec::new();
    let mut contig_stats: Vec<(f64, f64)> = Vec::new(); // (mean_cov, gc_frac)
    let mut fit_gc = Vec::new();
    let mut fit_rel_cov = Vec::new();

    for (idx, (_name, seq)) in sequences.iter().enumerate() {
        let cov_vec = &coverage[idx];
        let len = seq.len();
        let sum_cov: u64 = cov_vec.iter().map(|v| *v as u64).sum();
        let contig_mean = if len > 0 {
            sum_cov as f64 / len as f64
        } else {
            0.0
        };
        let gc_total = seq
            .iter()
            .filter(|b| matches!(**b, b'G' | b'g' | b'C' | b'c'))
            .count();
        let contig_gc_frac = if len > 0 {
            gc_total as f64 / len as f64
        } else {
            0.0
        };
        contig_stats.push((contig_mean, contig_gc_frac));
        let mut data = ContigData {
            gc: Vec::new(),
            cov: Vec::new(),
        };
        let mut pos = 0usize;
        while pos < len {
            let end = std::cmp::min(pos + window_size, len);
            let window_len = end - pos;
            let gc_count = seq[pos..end]
                .iter()
                .filter(|b| matches!(**b, b'G' | b'g' | b'C' | b'c'))
                .count();
            let gc_frac = gc_count as f64 / window_len as f64;
            let cov_sum: u64 = cov_vec[pos..end].iter().map(|v| *v as u64).sum();
            let mean_cov = cov_sum as f64 / window_len as f64;
            data.gc.push(gc_frac);
            data.cov.push(mean_cov);
            if contig_mean >= min_contig_cov && contig_mean > 0.0 {
                fit_gc.push(gc_frac);
                fit_rel_cov.push(mean_cov / contig_mean);
            }
            pos += window_size;
        }
        contig_data.push(data);
    }

    if fit_gc.is_empty() {
        return Err("No contigs with sufficient coverage for GC bias modelling".into());
    }
    let (spline, points) = fit_gc_bias_spline(&fit_gc, &fit_rel_cov);
    if let Some(p) = plot_path {
        plot_spline(&points, &spline, p)?;
    }

    let mut results = Vec::new();
    for (((name, _), data), (orig_cov, gc_frac)) in sequences
        .iter()
        .zip(contig_data.iter())
        .zip(contig_stats.iter())
    {
        if data.gc.is_empty() {
            results.push((name.clone(), 0.0, *orig_cov, gc_frac * 100.0));
            continue;
        }
        let mut sum_adj = 0.0;
        let mut count = 0usize;
        for (g, c) in data.gc.iter().zip(data.cov.iter()) {
            if let Some(pred) = spline.clamped_sample(*g) {
                if pred > 0.0 {
                    sum_adj += c / pred;
                    count += 1;
                }
            }
        }
        let adj_mean = if count > 0 {
            sum_adj / count as f64
        } else {
            0.0
        };
        results.push((name.clone(), adj_mean, *orig_cov, gc_frac * 100.0));
    }
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fit_gc_bias_spline_quadratic() {
        let gc = vec![0.0, 0.5, 1.0];
        let rel: Vec<f64> = gc.iter().map(|g| 1.0 + (g - 0.5) * (g - 0.5)).collect();
        let (spline, _pts) = fit_gc_bias_spline(&gc, &rel);
        for g in &gc {
            let pred = spline.clamped_sample(*g).unwrap();
            let true_val = 1.0 + (g - 0.5) * (g - 0.5);
            assert!((pred - true_val).abs() < 0.1);
        }
    }

    #[test]
    fn test_adjustment_restores_mean() {
        let gc = vec![0.0, 0.5, 1.0];
        let cov: Vec<f64> = gc
            .iter()
            .map(|g| 10.0 * (1.0 + (g - 0.5) * (g - 0.5)))
            .collect();
        let contig_mean = 10.0;
        let rel: Vec<f64> = cov.iter().map(|c| c / contig_mean).collect();
        let (spline, _pts) = fit_gc_bias_spline(&gc, &rel);
        let mut adj = Vec::new();
        for (g, c) in gc.iter().zip(cov.iter()) {
            let pred = spline.clamped_sample(*g).unwrap();
            adj.push(c / pred);
        }
        let mean_adj: f64 = adj.iter().sum::<f64>() / adj.len() as f64;
        assert!((mean_adj - 10.0).abs() < 0.1);
    }
}
