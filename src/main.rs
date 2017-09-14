#![feature(integer_atomics)]

use std::thread;
use std::io;
use std::sync::Arc;
use std::ops::Deref;
use std::collections::{HashSet, HashMap};
use std::hash::BuildHasherDefault;

extern crate twox_hash;
use twox_hash::XxHash;

extern crate fnv;
use fnv::FnvHashMap;

extern crate csv;

extern crate jobsteal;
use jobsteal::make_pool;

#[macro_use]
extern crate log;
extern crate env_logger;

extern crate serde;
#[macro_use]
extern crate serde_derive;

#[macro_use]
extern crate clap;

extern crate typed_arena;

mod find_genomes;
use find_genomes::{find_genomes, Gene, create_arenas};

mod gene_data;
use gene_data::AtomicGeneData;

fn main() {
    let arenas = create_arenas();
    let args = load_yaml!("../cli.yml");
    let args = clap::App::from_yaml(args).get_matches();
    let delimiter = args.value_of("delimiter").unwrap().as_bytes();
    if delimiter.len() != 1 {
        panic!("Expected 1 byte delimiter, but found multiple bytes: {:?}", delimiter);
    }
    let delimiter = delimiter[0];
    env_logger::init().unwrap();
    let thread_count: usize = args.value_of("threads").unwrap().parse().expect("Failed to parse thread count");
    let mut pool = make_pool(thread_count - 1).unwrap();
    let mut gene_counts: HashMap<Gene, Arc<AtomicGeneData>, BuildHasherDefault<XxHash>> = HashMap::default();
    let genomes = find_genomes(&arenas, args.value_of("directory").unwrap())
        .expect("Failed to find genomes")
        .map(|genome| {
            let genome = genome.expect("Failed to list and open genomes");
            info!("Found genome {}", genome.name);
            let gff = genome.gff_iter
                .map(|item| item.expect("Error reading from GFF file"))
                .map(|(accn, gene, range)| {
                    let counts_ref = gene_counts.entry(gene)
                        .or_insert_with(Default::default);
                    (accn, counts_ref.clone(), range)
                })
                .collect::<Vec<_>>();
            (genome.name, genome.blast_iter, gff)
        })
        .collect::<Vec<_>>();
    if genomes.is_empty() {
        warn!("Failed to find genomes to process");
    }
    info!("Created gene_counts HashMap");
    pool.scope(|scope| {
        for genome in genomes {
            scope.submit(move || {
                let our_thread = thread::current();
                let thread_name = our_thread.name().unwrap_or("[unknown]");
                info!("Thread {}: now working on: {}", thread_name, genome.0);
                let mut sequence_counts: FnvHashMap<u64, Vec<u8>> = HashMap::default();
                for item in genome.1 { // blast
                    let (accn, range) = item.expect("IO Error reading from BLAST file");
                    let sequence_count = sequence_counts.entry(accn).or_insert_with(Vec::new);
                    if range.end > sequence_count.len() {
                        sequence_count.resize(range.end, 0);
                    }
                    for index in range {
                        sequence_count[index] += 1;
                    }
                }
                let mut genes_encountered: HashSet<_, BuildHasherDefault<XxHash>> = HashSet::default();
                for (accn, data, range) in genome.2 { // gff (preprocessed)
                    let mut tmp_count: u64 = 0;
                    let mut index = range.start;
                    let sequence_count = match sequence_counts.get(&accn) {
                        Some(x) => x,
                        None => continue,
                    };
                    while index < range.end {
                        if let Some(&x) = sequence_count.get(index) {
                            if x > 0 {
                                let range_end = range.end;
                                // A bit hacky, but this should always work
                                let first_in_genome = genes_encountered.insert(data.deref() as *const _);
                                data.first_overlap(range, first_in_genome);
                                let x: u64 = x.into();
                                tmp_count += x;
                                // Loop unrolling might allow this code to be simpler
                                // However, assuming it gets unrolled might be a bad idea
                                // This way it'll always be fast
                                for i in (index + 1)..range_end {
                                    if let Some(&x) = sequence_count.get(i) {
                                        let x: u64 = x.into();
                                        tmp_count += x;
                                    }
                                }
                                break;
                            }
                        }
                        index += 1;
                    }
                    data.total_overlap(tmp_count);
                }
            });
        }
    });
    info!("Finished counting");
    let stdout = io::stdout();
    let mut writer = csv::WriterBuilder::new()
        .delimiter(delimiter)
        .has_headers(false) // we'll write them manually
        .from_writer(stdout.lock());
    writer.write_record(&["name", "product", "total_overlap",
                          "genome_count", "start_avg", "start_stdev",
                          "end_avg", "end_stdev", "length_avg"])
        .expect("Failed to write headers");
    let gene_counts = gene_counts
        .into_iter()
        .map(|(gene, data)| {
            match Arc::try_unwrap(data) {
                Ok(data) => data.finalize(gene),
                Err(_) => panic!("Internal error: duplicate count for gene {:?}", gene),
            }
        })
        .filter(|data| data.total_overlap != 0);
    if args.is_present("sort") {
        let mut gene_counts = gene_counts.collect::<Vec<_>>();
        // descending sort by count
        gene_counts.sort_by(|ref a, ref b| b.total_overlap.cmp(&a.total_overlap));
        for count in gene_counts {
            writer.serialize(count).expect("Failed to write gene counts");
        }
    } else {
        for count in gene_counts {
            writer.serialize(count).expect("Failed to write gene counts");
        }
    }
    writer.flush().expect("IO Error flushing CSV");
}
