use std::str;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::collections::{HashMap, VecDeque};
use std::fs::{self, File, read_dir, ReadDir};
use std::io::{self, BufRead, BufReader};
use std::marker::PhantomData;
use std::fmt;

use csv;
use serde::{Serialize, Serializer};
use serde::ser::SerializeTuple;
use typed_arena::Arena;

pub struct Arenas {
    gene_base: Arena<GeneBase>,
}

pub struct LabelString<D: fmt::Display + ?Sized, T>(Option<&'static str>, T, PhantomData<D>);

impl<'a, D: fmt::Display + ?Sized, T> LabelString<D, T>
    where Self: fmt::Display
{
    fn new(label: Option<&'static str>, contents: T) -> LabelString<D, T> {
        LabelString(label, contents, PhantomData {})
    }
}

impl From<&'static str> for LabelString<str, [&'static str; 1]> {
    fn from(message: &'static str) -> Self {
        LabelString(None, [message], PhantomData {})
    }
}

impl<'a, D: fmt::Display + 'a + ?Sized, T> fmt::Display for LabelString<D, T>
    where T: AsRef<[&'a D]>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if let Some(label) = self.0 {
            f.write_str(label)?;
            f.write_str("(")?;
        }
        for (i, ref item) in self.1.as_ref().iter().enumerate() {
            if i != 0 {
                f.write_str(", ")?;
            }
            item.fmt(f)?;
        }
        if self.0.is_some() {
            f.write_str(")")?;
        }
        Ok(())
    }
}

impl<'a, D: fmt::Display + ?Sized, T> Serialize for LabelString<D, T> where Self: fmt::Display {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.collect_str(&self)
    }
}

// We don't implement `Clone` so we don't accidentally clone it
#[derive(PartialEq, Eq, Hash, Debug)]
pub enum GeneBase {
    Normal {
        name: Option<String>,
        product: String,
    },
    Hypothetical,
}

impl GeneBase {
    fn get_name(&self) -> &str {
        match *self {
            GeneBase::Normal { name: Some(ref name), .. } => &name,
            _ => "UNKNOWN",
        }
    }

    fn get_product(&self) -> &str {
        match *self {
            GeneBase::Normal { ref product, .. } => &product,
            GeneBase::Hypothetical => "hypothetical protein",
        }
    }

    fn is_normal(&self) -> bool {
        match *self {
            GeneBase::Normal { .. } => true,
            GeneBase::Hypothetical => false,
        }
    }
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum RelativePosition {
    Before,
    After,
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct GeneInfo<'a> {
    hypothetical: Option<RelativePosition>,
    base: &'a GeneBase,
}

impl<'a> From<&'a GeneBase> for GeneInfo<'a> {
    fn from(base: &'a GeneBase) -> Self {
        GeneInfo {
            hypothetical: None,
            base: base,
        }
    }
}

impl<'a> GeneInfo<'a> {
    fn get_label(&self) -> Option<&'static str> {
        match self.hypothetical {
            Some(RelativePosition::Before) => Some("HypotheticalBefore"),
            Some(RelativePosition::After) => Some("HypotheticalAfter"),
            None => None,
        }
    }

    fn display_name(&self) -> LabelString<str, [&str; 1]> {
        LabelString::new(self.get_label(), [self.base.get_name()])
    }

    fn display_product(&self) -> LabelString<str, [&str; 1]> {
        LabelString::new(self.get_label(), [self.base.get_product()])
    }
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum GeneDescriptor<'a> {
    Normal(GeneInfo<'a>),
    Start,
    End,
}

impl<'a> From<&'a GeneBase> for GeneDescriptor<'a> {
    fn from(base: &'a GeneBase) -> Self {
        GeneDescriptor::Normal(base.into())
    }
}

impl<'a> From<GeneInfo<'a>> for GeneDescriptor<'a> {
    fn from(info: GeneInfo<'a>) -> Self {
        GeneDescriptor::Normal(info)
    }
}

impl<'a> GeneDescriptor<'a> {
    fn display_name(&self) -> LabelString<str, [&str; 1]> {
        match *self {
            GeneDescriptor::Normal(ref info) => info.display_name(),
            GeneDescriptor::Start => LabelString::new(None, ["START"]),
            GeneDescriptor::End => LabelString::new(None, ["END"]),
        }
    }

    fn display_product(&self) -> LabelString<str, [&str; 1]> {
        match *self {
            GeneDescriptor::Normal(ref info) => info.display_product(),
            GeneDescriptor::Start => LabelString::new(None, ["START"]),
            GeneDescriptor::End => LabelString::new(None, ["END"]),
        }
    }
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum Gene<'a> {
    Normal(GeneDescriptor<'a>),
    Space(GeneDescriptor<'a>, GeneDescriptor<'a>),
}

impl<'a> Serialize for Gene<'a> {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let mut tuple = s.serialize_tuple(2)?;
        match *self {
            Gene::Normal(ref base) => {
                tuple.serialize_element(&base.display_name())?;
                tuple.serialize_element(&base.display_product())?;
            }
            Gene::Space(ref before, ref after) => {
                tuple.serialize_element(&LabelString::<LabelString<_, _>, _>::new(Some("Between"), &[&before.display_name(), &after.display_name()]))?;
                tuple.serialize_element(&LabelString::<LabelString<_, _>, _>::new(Some("Between"), &[&before.display_product(), &after.display_product()]))?;
            }
        }
        tuple.end()
    }
}

fn get_accn(bytes: &[u8]) -> u64 {
    let mut last_index = bytes.len();
    for i in bytes.len()..0 {
        match bytes[i] {
            b'0' => {}, // slight optimization: strip leading 0s
            b'1'...b'9' => last_index = i,
            _ => break,
        }
    }
    if last_index == bytes.len() {
        warn!("Found ACCN without trailing numbers: {}", String::from_utf8_lossy(bytes));
        return 0;
    }
    str::from_utf8(&bytes[last_index..]).unwrap().parse().expect("Found ACCN larger than 2^64 - 1")
}

struct GffBaseIter<'a> {
    arenas: &'a Arenas,
    reader: BufReader<File>,
    buffer: Vec<u8>,
}

const HYPOTHETICAL: &'static GeneBase = &GeneBase::Hypothetical;

impl<'a> Iterator for GffBaseIter<'a> {
    type Item = io::Result<(u64, &'a GeneBase, Range<usize>)>;

    fn next(&mut self) -> Option<io::Result<(u64, &'a GeneBase, Range<usize>)>> {
        loop {
            self.buffer.clear();
            match self.reader.read_until(b'\n', &mut self.buffer) {
                Ok(x) => x,
                Err(e) => return Some(Err(e)),
            };
            let line = &mut self.buffer;
            if line.is_empty() {
                return None;
            }
            if line[0] == b'#' {
                continue;
            }
            let mut line = line.splitn(9, |x| *x == b'\t');
            let accn = line.next().unwrap();
            let start: usize = match line.nth(2) {
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
            let range = start..(end + 1);
            // This iterator has already been paritally consumed
            // 7th index - next() - nth(2) - next()
            // 7 - 1 - (2 + 1) - 1
            // 3
            let descriptors = match line.nth(3) {
                Some(x) => x.split(|x| *x == b';'),
                None => continue,
            };
            let accn = get_accn(accn);
            let mut name = None;
            let mut product = None;
            for descriptor in descriptors {
                let mut parts = descriptor.splitn(2, |x| *x == b'=');
                let key = match parts.next() {
                    Some(x) => x,
                    None => continue,
                };
                let mut value = match parts.next() {
                    Some(x) => x,
                    None => {
                        warn!("Descriptor had a key but no value");
                        continue;
                    }
                };
                if value.last() == Some(&b'\n') {
                    value = &value[0..(value.len() - 1)];
                }
                let to_set = if key == b"Name" {
                    Some(&mut name)
                } else if key == b"product" {
                    if value == b"hypothetical protein" {
                        return Some(Ok((accn, HYPOTHETICAL, range)));
                    }
                    Some(&mut product)
                } else {
                    None
                };
                if let Some(var) = to_set {
                    let value = String::from_utf8_lossy(value).into();
                    *var = Some(value);
                }
            }
            let product = match product {
                Some(x) => x,
                None => {
                    warn!("Encountered a gene with no product, name: {:?}", name);
                    continue;
                }
            };
            let gene_base = self.arenas.gene_base.alloc(GeneBase::Normal {
                name: name,
                product: product,
            });
            return Some(Ok((accn, gene_base, range)));
        }
    }
}

pub struct GffIter<'a> {
    base: GffBaseIter<'a>,
    gene_buffer: VecDeque<(u64, &'a GeneBase, Range<usize>)>,
    next_index: Option<(u64, usize)>,
    // ACCN, base
    prev_useful_gene: Option<(u64, &'a GeneBase)>,
    next_useful_gene: Option<(u64, &'a GeneBase)>,
}

impl<'a> GffIter<'a> {
    fn next_useful_gene(&mut self) -> io::Result<Option<(u64, &'a GeneBase)>> {
        if let Some(gene) = self.next_useful_gene {
            return Ok(Some(gene));
        }
        let mut buf_iter = self.gene_buffer.iter().cloned().map(|x| (true, Ok(x))).collect::<Vec<_>>().into_iter();
        let mut first_accn = None;
        let mut next_normal = None;
        let mut next_hypothetical = None;
        while let Some((from_buf, base)) = buf_iter.next().or_else(|| self.base.next().map(|x| (false, x))) {
            let base = base?;
            if !from_buf {
                self.gene_buffer.push_back(base.clone());
            }
            let base = (base.0, base.1);
            match first_accn {
                None => first_accn = Some(base.0),
                Some(x) if base.0 != x => break,
                _ => {},
            }
            if base.1.is_normal() {
                next_normal = Some(base);
                break;
            } else {
                next_hypothetical = Some(next_hypothetical.unwrap_or(base));
            }
        }
        let next_useful = next_normal.or(next_hypothetical);
        self.next_useful_gene = next_useful;
        Ok(next_useful)
    }
}

impl<'a> Iterator for GffIter<'a> {
    type Item = io::Result<(u64, Gene<'a>, Range<usize>)>;

    fn next(&mut self) -> Option<io::Result<(u64, Gene<'a>, Range<usize>)>> {
        if let Some(base) = self.gene_buffer.pop_front().map(Ok).or_else(|| self.base.next()) {
            let base = match base {
                Ok(x) => x,
                Err(e) => return Some(Err(e)),
            };
            if let Some((accn, _)) = self.prev_useful_gene {
                if base.0 != accn {
                    self.prev_useful_gene = None;
                }
            }
            let ret = (|| {
                if let Some(next_index) = self.next_index {
                    if base.0 == next_index.0 && base.2.start > next_index.1 {
                        self.gene_buffer.push_front(base.clone());
                        let before = self.prev_useful_gene
                            .map(|(_, x)| x)
                            .map(Into::into)
                            .unwrap_or(GeneDescriptor::Start);
                        let after = if base.1.is_normal() {
                            base.1.into()
                        } else {
                            match self.next_useful_gene() {
                                Ok(Some(x)) if x.0 == next_index.0 => x.1.into(),
                                Ok(_) => GeneDescriptor::End,
                                Err(e) => return Err(e),
                            }
                        };
                        return Ok((base.0, Gene::Space(before, after), next_index.1..base.2.start));
                    }
                }
                if base.1.is_normal() {
                    return Ok((base.0, Gene::Normal(base.1.into()), base.2.clone()));
                }
                if let Some((_, before)) = self.prev_useful_gene {
                    if before.is_normal() {
                        return Ok((base.0, Gene::Normal(GeneInfo {
                            hypothetical: Some(RelativePosition::After),
                            base: before,
                        }.into()), base.2.clone()));
                    }
                }
                if let Some((accn, next)) = self.next_useful_gene()? {
                    if accn == base.0 {
                        return Ok((base.0, Gene::Normal(GeneInfo {
                            hypothetical: Some(RelativePosition::Before),
                            base: next,
                        }.into()), base.2.clone()));
                    }
                }
                error!("Encountered ACCN with only one gene, a hypothetical protein");
                Ok((base.0, Gene::Normal(GeneDescriptor::Start), base.2.clone()))
            })();
            let prev_normal = self.prev_useful_gene.as_ref().map(|x| x.1.is_normal()).unwrap_or(false);
            if base.1.is_normal() || !prev_normal {
                self.prev_useful_gene = Some((base.0, base.1));
            }
            self.next_index = Some((base.0, base.2.end));
            Some(ret)
        } else {
            None
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
    type Item = csv::Result<(u64, Range<usize>)>;

    fn next(&mut self) -> Option<csv::Result<(u64, Range<usize>)>> {
        self.0.next().map(|line| {
            line.map(|line| {
                let accn = get_accn(line.subject.as_bytes());
                (accn, line.s_start..(line.s_end + 1))
            })
        })
    }
}

pub struct Genome<'a> {
    pub name: String,
    pub gff_iter: GffIter<'a>,
    pub blast_iter: BlastIter,
}

pub struct GenomeIter<'a> {
    arenas: &'a Arenas,
    dir: ReadDir,
    found: HashMap<String, (Option<PathBuf>, Option<PathBuf>)>,
}

const GFF_POSTFIX: &'static str = ".gff";
const BLAST_POSTFIX: &'static str = ".bla";

const OPTIONAL_PREFIX: &'static str = "R_";

impl<'a> Iterator for GenomeIter<'a> {
    type Item = io::Result<Genome<'a>>;

    fn next(&mut self) -> Option<io::Result<Genome<'a>>> {
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
            let file_path = entry.path();
            let name = match entry.file_name().into_string() {
                Ok(x) => x,
                Err(_) => {
                    warn!("Found non-UTF8 file");
                    continue;
                }
            };
            let mut name = name.as_str();
            let is_blast = if name.ends_with(GFF_POSTFIX) {
                name = &name[..(name.len() - GFF_POSTFIX.len())];
                false
            } else if name.ends_with(BLAST_POSTFIX) {
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
                    warn!("Encountered duplicate file: {:?}", file_path);
                }
                *ours = Some(file_path);
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
                        base: GffBaseIter {
                            arenas: self.arenas,
                            reader: BufReader::new(gff_file),
                            buffer: Vec::new(),
                        },
                        gene_buffer: VecDeque::new(),
                        next_index: None,
                        prev_useful_gene: None,
                        next_useful_gene: None,
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

pub fn find_genomes<'a, P: AsRef<Path>>(arenas: &'a Arenas, path: P) -> io::Result<GenomeIter<'a>> {
    Ok(GenomeIter {
        arenas: arenas,
        dir: read_dir(path)?,
        found: HashMap::new(),
    })
}

pub fn create_arenas() -> Arenas {
    Arenas {
        gene_base: Arena::new(),
    }
}
