use std::str;
use std::ops::Range;
use std::path::Path;
use std::collections::HashMap;
use std::fs::{self, File, read_dir, ReadDir};
use std::io::{self, BufRead, BufReader};

use csv;

#[derive(PartialEq, Eq, Hash)]
pub struct Gene {
    pub name: Option<String>,
    pub product: String,
}

pub struct GffIter {
    reader: BufReader<File>,
    buffer: Vec<u8>,
}

impl Iterator for GffIter {
    type Item = io::Result<(Gene, Range<usize>)>;

    fn next(&mut self) -> Option<io::Result<(Gene, Range<usize>)>> {
        loop {
            self.buffer.clear();
            match self.reader.read_until(b'\n', &mut self.buffer) {
                Ok(x) => x,
                Err(e) => return Some(Err(e)),
            };
            let line = &mut self.buffer;
            if line.get(0) == Some(&b'#') {
                continue;
            }
            let mut line = line.splitn(8, |x| *x == b'\t');
            let start: usize = match line.nth(3) {
                Some(x) => match str::from_utf8(x).ok().and_then(|x| x.parse().ok()) {
                    Some(x) => x,
                    None => {
                        warn!("Failed to parse start: {:?}", x);
                        continue;
                    }
                },
                None => continue,
            };
            let end: usize = match line.next() {
                Some(x) => match str::from_utf8(x).ok().and_then(|x| x.parse().ok()) {
                    Some(x) => x,
                    None => {
                        warn!("Failed to parse end: {:?}", x);
                        continue;
                    }
                },
                None => continue,
            };
            // This iterator has already been paritally consumed
            // 7th index - nth(3) - next()
            // 7 - (3 + 1) - 1
            // 3
            let descriptors = match line.nth(3) {
                Some(x) => x.split(|x| *x == b';'),
                None => continue,
            };
            let mut name = None;
            let mut product = None;
            for descriptor in descriptors {
                let mut parts = descriptor.splitn(2, |x| *x == b'=');
                let key = match parts.next() {
                    Some(x) => x,
                    None => continue,
                };
                let value = match parts.next() {
                    Some(x) => x,
                    None => {
                        warn!("Descriptor had a key but no value");
                        continue;
                    }
                };
                if key == b"Name" {
                    name = Some(String::from_utf8_lossy(value).into_owned());
                } else if key == b"product" {
                    product = Some(String::from_utf8_lossy(value).into_owned());
                }
            }
            let product = match product {
                Some(x) => x,
                None => {
                    warn!("Encountered a gene with no product, name: {:?}", name);
                    continue;
                }
            };
            return Some(Ok((Gene {
                name: name,
                product: product,
            }, start..end)));
        }
    }
}

// a lot of these fields are numbers,
// but we only parse what we need
#[allow(dead_code)]
#[derive(Deserialize)]
struct BlastLine {
    query: String,
    subject: String,
    percent_id: String,
    align_len: String,
    mismatches: String,
    gap_openings: String,
    q_start: String,
    q_end: String,
    s_start: usize,
    s_end: usize,
    e_value: String,
    bit_score: String,
}

pub struct BlastIter(csv::DeserializeRecordsIntoIter<File, BlastLine>);

impl Iterator for BlastIter {
    type Item = csv::Result<Range<usize>>;

    fn next(&mut self) -> Option<csv::Result<Range<usize>>> {
        self.0.next().map(|line| {
            line.map(|line| {
                (line.s_start..line.s_end)
            })
        })
    }
}

pub struct Genome {
    pub name: String,
    pub gff_iter: GffIter,
    pub blast_iter: BlastIter,
}

pub struct GenomeIter {
    dir: ReadDir,
    found: HashMap<String, (Option<String>, Option<String>)>,
}

const GFF_POSTFIX: &'static str = ".gff";
const BLAST_POSTFIX: &'static str = ".bla";

const OPTIONAL_PREFIX: &'static str = "R";

impl Iterator for GenomeIter {
    type Item = io::Result<Genome>;

    fn next(&mut self) -> Option<io::Result<Genome>> {
        loop {
            let entry = match self.dir.next() {
                Some(Ok(x)) => x,
                Some(Err(e)) => return Some(Err(e)),
                None => return None,
            };
            let metadata = match fs::metadata(entry.path()) {
                Ok(x) => x,
                Err(e) => return Some(Err(e)),
            };;
            if !metadata.is_file() {
                continue;
            }
            let file_name = match entry.file_name().into_string() {
                Ok(x) => x,
                Err(_) => {
                    warn!("Found non-UTF8 file");
                    continue;
                }
            };
            let name = file_name.clone();
            let mut name = name.as_str();
            let is_blast = if file_name.ends_with(GFF_POSTFIX) {
                name = &name[..(name.len() - GFF_POSTFIX.len())];
                false
            } else if file_name.ends_with(BLAST_POSTFIX) {
                name = &name[..(name.len() - BLAST_POSTFIX.len())];
                true
            } else {
                continue;
            };
            if name.starts_with(OPTIONAL_PREFIX) {
                name = &name[OPTIONAL_PREFIX.len()..];
            }
            let genome_files = self.found.entry(name.to_string()).or_insert((None, None));
            {
                let ours = if is_blast {
                    &mut genome_files.0
                } else {
                    &mut genome_files.1
                };
                if ours.is_some() {
                    warn!("Encountered duplicate file: {}", file_name);
                }
                *ours = Some(file_name);
            }
            if let Some(ref blast_file) = genome_files.0 {
                if let Some(ref gff_file) = genome_files.1 {
                    let gff_file = match File::open(gff_file) {
                        Ok(x) => x,
                        Err(e) => return Some(Err(e)),
                    };
                    let blast_file = match File::open(blast_file) {
                        Ok(x) => x,
                        Err(e) => return Some(Err(e)),
                    };
                    let gff_iter = GffIter {
                        reader: BufReader::new(gff_file),
                        buffer: Vec::new(),
                    };
                    let blast_reader = csv::ReaderBuilder::new()
                        .delimiter(b'\t')
                        .has_headers(false)
                        .from_reader(blast_file);
                    let blast_iter = BlastIter(blast_reader.into_deserialize());
                    return Some(Ok(Genome {
                        name: name.to_string(),
                        blast_iter: blast_iter,
                        gff_iter: gff_iter,
                    }));
                }
            }
        }
    }
}

pub fn find_genomes<P: AsRef<Path>>(path: P) -> io::Result<GenomeIter> {
    Ok(GenomeIter {
        dir: read_dir(path)?,
        found: HashMap::new(),
    })
}
