#![feature(integer_atomics)]

use std::io::{self, Write};
use std::sync::Arc;
use std::sync::atomic::{self, AtomicU16};
use std::collections::HashMap;

extern crate csv;

extern crate jobsteal;
use jobsteal::make_pool;

#[macro_use]
extern crate log;
extern crate env_logger;

#[macro_use]
extern crate serde_derive;

#[macro_use]
extern crate clap;

mod find_genomes;
use find_genomes::{find_genomes, Gene};

fn print_counts<I: IntoIterator<Item=(Gene, u16)>, W: Write>(writer: &mut csv::Writer<W>, gene_counts: I) {
    for (gene, count) in gene_counts {
        let tmp;
        let name = match gene.name {
            Some(name) => {
                tmp = name;
                tmp.as_str()
            }
            None => "unknown",
        };
        writer.write_record(&[name, gene.product.as_str(), count.to_string().as_str()]).expect("IO Error writing to CSV");
    }
}

fn main() {
    let args = load_yaml!("../cli.yml");
    let args = clap::App::from_yaml(args).get_matches();
    env_logger::init().unwrap();
    let thread_count: usize = args.value_of("threads").unwrap().parse().expect("Failed to parse thread count");
    let mut pool = make_pool(thread_count - 1).unwrap();
    let genomes = find_genomes(args.value_of("directory").unwrap()).expect("Failed to find genomes");
    let mut gene_counts: HashMap<Gene, Arc<AtomicU16>> = HashMap::new();
    let genomes = genomes.map(|genome| {
        let genome = genome.expect("Failed to list and open genomes");
        info!("Found genome {}", genome.name);
        let gff = genome.gff_iter
            .map(|item| item.expect("Error reading from GFF file"))
            .map(|(gene, range)| {
                let count_ref = gene_counts.entry(gene).or_insert_with(|| Arc::new(AtomicU16::new(0)));
                (count_ref.clone(), range)
            })
            .collect::<Vec<_>>();
        (genome.blast_iter, gff)
    }).collect::<Vec<_>>();
    if genomes.is_empty() {
        warn!("Failed to find genomes to process");
    }
    pool.scope(|scope| {
        for genome in genomes {
            scope.submit(move || {
                let mut sequence_count: Vec<u16> = Vec::new();
                for item in genome.0 {
                    let range = item.expect("IO Error reading from BLAST file");
                    if range.end > sequence_count.len() {
                        sequence_count.resize(range.end, 0);
                    }
                    for index in range {
                        sequence_count[index] += 1;
                    }
                }
                for (count, range) in genome.1 {
                    let mut tmp_count = 0;
                    for index in range {
                        if let Some(x) = sequence_count.get(index) {
                            tmp_count += *x;
                        }
                    }
                    count.fetch_add(tmp_count, atomic::Ordering::Relaxed);
                }
            });
        }
    });
    let stdout = io::stdout();
    let mut writer = csv::Writer::from_writer(stdout.lock());
    writer.write_record(&["name", "product", "count"]).expect("IO Error writing to CSV");
    let gene_counts = gene_counts
        .into_iter()
        .map(|(gene, count)| (gene, Arc::try_unwrap(count).unwrap().into_inner()));
    if args.is_present("sort") {
        let mut gene_counts = gene_counts.collect::<Vec<_>>();
        // descending sort by count
        gene_counts.sort_by(|ref a, ref b| b.1.cmp(&a.1));
        print_counts(&mut writer, gene_counts);
    } else {
        print_counts(&mut writer, gene_counts);
    }
    writer.flush().expect("IO Error flushing CSV");
}
