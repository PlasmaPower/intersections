use std::f64;
use std::ops::Range;
use std::sync::atomic::{AtomicU32, AtomicU64};
use std::sync::atomic::Ordering::Relaxed;

use find_genomes::Gene;

#[derive(Default)]
pub struct AtomicGeneData {
    overlap: AtomicU64,
    genome_count: AtomicU32,
    data_count: AtomicU32,
    start_sum: AtomicU32,
    start_squared_sum: AtomicU64,
    end_sum: AtomicU32,
    end_squared_sum: AtomicU64,
    length_sum: AtomicU32,
}

#[derive(Serialize)]
pub struct FinalGeneData<'a> {
    pub gene: Gene<'a>,
    pub total_overlap: u64,
    pub genome_count: u32,
    pub start_avg: f64,
    pub start_stdev: f64,
    pub end_avg: f64,
    pub end_stdev: f64,
    pub length_avg: f64,
}

fn stdev(n: u32, avg: f64, squared_sum: u64) -> f64 {
    if n <= 1 {
        return f64::NAN;
    }
    ((((squared_sum as f64) - (n as f64)*avg*avg) as f64) / ((n - 1) as f64)).sqrt()
}

impl AtomicGeneData {
    pub fn first_overlap(&self, gene_range: Range<usize>, first_in_genome: bool) {
        if first_in_genome {
            self.genome_count.fetch_add(1, Relaxed);
        }
        self.data_count.fetch_add(1, Relaxed);
        self.start_sum.fetch_add(gene_range.start as u32, Relaxed);
        self.start_squared_sum.fetch_add((gene_range.start as u64) * (gene_range.start as u64), Relaxed);
        self.end_sum.fetch_add(gene_range.end as u32, Relaxed);
        self.end_squared_sum.fetch_add((gene_range.end as u64) * (gene_range.end as u64), Relaxed);
        self.length_sum.fetch_add((gene_range.end - gene_range.start) as u32, Relaxed);
    }

    pub fn total_overlap(&self, overlap: u64) {
        self.overlap.fetch_add(overlap, Relaxed);
    }

    pub fn finalize(self, gene: Gene) -> FinalGeneData {
        let n = self.data_count.into_inner();
        let n_float = n as f64;
        let start_avg = (self.start_sum.into_inner() as f64) / n_float;
        let end_avg = (self.end_sum.into_inner() as f64) / n_float;
        FinalGeneData {
            gene: gene,
            total_overlap: self.overlap.into_inner(),
            genome_count: self.genome_count.into_inner(),
            start_avg: start_avg,
            end_avg: end_avg,
            start_stdev: stdev(n, start_avg, self.start_squared_sum.into_inner()),
            end_stdev: stdev(n, end_avg, self.end_squared_sum.into_inner()),
            length_avg: (self.length_sum.into_inner() as f64) / n_float,
        }
    }
}
