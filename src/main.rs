#![feature(integer_atomics)]

use std::thread;
use std::io::{self, Write};
use std::sync::Arc;
use std::sync::atomic::{self, AtomicU32};
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

fn print_counts<I: IntoIterator<Item=(Gene, u32)>, W: Write>(writer: &mut csv::Writer<W>, gene_counts: I) {
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
    let mut gene_counts: HashMap<Gene, Arc<AtomicU32>> = HashMap::new();
    let genomes = genomes.map(|genome| {
        let genome = genome.expect("Failed to list and open genomes");
        info!("Found genome {}", genome.name);
        let gff = genome.gff_iter
            .map(|item| item.expect("Error reading from GFF file"))
            .map(|(gene, range)| {
                let count_ref = gene_counts.entry(gene).or_insert_with(|| Arc::new(AtomicU32::new(0)));
                (count_ref.clone(), range)
            })
            .collect::<Vec<_>>();
        (genome.name, genome.blast_iter, gff)
    }).collect::<Vec<_>>();
    if genomes.is_empty() {
        warn!("Failed to find genomes to process");
    }
    info!("Created gene_counts HashMap");
    pool.scope(|scope| {
        for genome in genomes {
            scope.submit(move || {
                let our_thread = thread::current();
                let thread_name = our_thread.name().unwrap_or("[unknown]");
                info!("{}: now working on: {}", thread_name, genome.0);
                // TODO this can probably be changed to Vec<u8>
                // However, we want to avoid silent overflows
                let mut sequence_count: Vec<u16> = Vec::new();
                for item in genome.1 {
                    let range = item.expect("IO Error reading from BLAST file");
                    if range.end > sequence_count.len() {
                        sequence_count.resize(range.end, 0);
                    }
                    for index in range {
                        sequence_count[index] += 1;
                    }
                }
                for (count, range) in genome.2 {
                    let mut tmp_count: u32 = 0;
                    for index in range {
                        if let Some(x) = sequence_count.get(index) {
                            tmp_count += (*x).into();
                        }
                    }
                    count.fetch_add(tmp_count, atomic::Ordering::Relaxed);
                }
            });
        }
    });
    info!("Finished counting");
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
