# gene-seq-intersections

Finds intersections between genes and sequences.
More documentation coming soon™. Run with `--help` for CLI options.

# Intersections
This program finds the overlap of sequences and genes using format 6 blastn output (http://www.metagenomics.wiki/tools/blast/blastn-output-format-6)
  
  ```
  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
  Query_1	accn|JISN01000002	100.000	28	0	0	29	56	37930	37957	1.32e-08	52.8
  ```

and gff3 output (from prokka)
  
  ```
  ##gff-version 3
  ##sequence-region accn_JISN01000001 1 334949
  ...
  accn_JISN01000001	Prodigal:2.6	CDS	240	2849	.	+	0	ID=NKHGEDLF_00001;Name=clpB;gene=clpB;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:Q7A6G6;locus_tag=NKHGEDLF_00001;product=Chaperone protein ClpB
  ...
  >accn_JISN01000001
  AATTAATTATCGACCAAGAAAGTGTTTAAATTGGAAGTTTCCTTATGAAGTTTTAT
  ...
  ```

Lines 9 and 10 of the blastn output are compared to lines 4 and 5 of the gff3 file (section type 2) for overlap. Any number of bla files can be intersected with an equal number of MATCHING gff files.

## Prerequisites
Folder of .bla files and .gff files MATCHED by NAME (I.E. genome1.bla genome1.gff genome2.bla genome2.gff). Bla files are files created in blastn format 6 by the blasting of one or more sequences against the respective genome. Gff3 files are created (for example) by prokka v1.12 (http://www.vicbioinformatics.com/software.prokka.shtml) for a respective genome.

## Installing
First download rust (instructions from https://rustup.rs/)

```
curl https://sh.rustup.rs -sSf | sh
```

Then download the crate for intersections

```
cargo +nightly install sequence-intersections
```

Intersections can then be found in ~/.cargo/bin/
If a previous version of intersections already exists in the directory use

```
cargo +nightly install -f sequence-intersections
```

## Output and Options

| Column | Description |
| --- | --- |
| name | Name of gene according to gff file. Regions between two genes are denoted Between(GeneNameBefore, GeneNameAfter). Hypothetical proteins are denoted HypotheticalAfter(GeneName) or HypotheticalBefore(GeneName) |
| product | Product of gene according to gff file. Same style as name. |
| total_overlap | Amount of sequence which intersected at this gene. If a sequence of 31 in the blast in put file completely overlapped with this gene (IE blast was in ID_1 and spanned 1000-1031 and the gene was in ID_1 and spanned 1000-1500) then the total_overlap for this gene would add +31. |
| genome_count | The number of genomes which had at least one sequence overlap this gene with at least 1 total_overlap. |
| start_avg | The average start for this gene according to the gff file. | 
| start_stdev | The standard deviation of the start of this gene. |
| end_avg | The average end for this gene according to the gff file. | 
| end_stdev | The standard deviation of the end of this gene. |
| length_avg | The average span of each gene (# of nucleotides long). Is not related to start or end location but only length of the gene. | 

## Example
Example blast and gff intersections at: https://github.com/dUmich/intersections-example

## Errors
Run with this command preceding to get warnings

```
RUST_LOG=warn 
```

## Built with

## Versioning 

## Authors
* **Lee Bousfield** - *Free-lance code wizard* - [PlasmaPower](https://github.com/PlasmaPower)
* **Daniel Harris** - *Researcher, Snitkin Lab, University of Michigan* - [dUmich](https://github.com/dUmich)
