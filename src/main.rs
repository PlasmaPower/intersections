use std::io;
use std::sync::Mutex;
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
use find_genomes::find_genomes;

fn main() {
    let args = load_yaml!("../cli.yml");
    let args = clap::App::from_yaml(args).get_matches();
    env_logger::init().unwrap();
    let mut pool = make_pool(0).unwrap();
    let genomes = find_genomes(args.value_of("directory").unwrap()).expect("Failed to find genomes");
    let gene_counts = Mutex::new(HashMap::new());
    pool.scope(|scope| {
        let gene_counts = &gene_counts;
        let mut found_genome = false;
        for genome in genomes {
            found_genome = true;
            let genome = genome.expect("Filed to open genomes");
            info!("Found genome {}", genome.name);
            scope.submit(move || {
                let mut sequence_count: Vec<u16> = Vec::new();
                for item in genome.blast_iter {
                    let range = item.expect("IO Error reading from BLAST file");
                    if range.end < sequence_count.len() {
                        sequence_count.resize(range.end, 0);
                    }
                    for index in range {
                        sequence_count[index] += 1;
                    }
                }
                let gff_iter = genome.gff_iter;
                let res = (|| -> csv::Result<()> {
                    let mut gene_counts = gene_counts.lock().unwrap();
                    for item in gff_iter {
                        let (gene, range) = item?;
                        let mut count = 0;
                        for index in range {
                            if let Some(x) = sequence_count.get(index) {
                                count += *x;
                            }
                        }
                        *gene_counts.entry(gene).or_insert(0) += count;
                    }
                    Ok(())
                })();
                // as to not break the mutex
                res.expect("IO Error reading from GFF file");
            });
        }
        if !found_genome {
            warn!("Failed to find genomes to process");
        }
    });
    let stdout = io::stdout();
    let mut writer = csv::Writer::from_writer(stdout.lock());
    writer.write_record(&["name", "product", "count"]).expect("IO Error writing to CSV");
    let gene_counts = gene_counts.into_inner().unwrap();
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
    writer.flush().expect("IO Error flushing CSV");
}
