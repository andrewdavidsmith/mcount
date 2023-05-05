/* MIT License
 *
 * Copyright (c) 2023 Andrew Smith
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

use std::error::Error;
use std::io::{prelude::*,BufWriter};

use std::collections::HashSet;
use std::collections::HashMap;

use rust_htslib::{bam,bam::Read,bam::Reader};

use msite::MSite;

pub fn read_fasta(
    filename: &String
) -> Result<(Vec<String>, Vec<Vec<u8>>), std::io::Error> {

    let f = std::fs::File::open(&filename).expect("bad fasta file");
    let reader = std::io::BufReader::new(f);

    let mut names = Vec::new();
    let mut seqs = Vec::new();

    for line in reader.lines() {
        let line = line.expect("bad fasta line");
        let line = line.as_bytes();

        if line[0] == ('>' as u8) {
            let x = match line.iter().position(|&x| x.is_ascii_whitespace()) {
                Some(x) => x,
                None => line.len(),
            };
            names.push(String::from_utf8_lossy(&line[1..x]).into_owned());
            seqs.push(Vec::new());
        }
        else {
            let l = seqs.len() - 1;
            seqs[l].extend_from_slice(line);
        }
    }
    Ok((names, seqs))
}

pub fn is_cytosine(base: u8) -> bool {
    base == b'C' || base == b'c'
}

pub fn is_guanine(base: u8) -> bool {
    base == b'G' || base == b'g'
}

pub fn is_ddg(s: &Vec<u8>, i: usize) -> bool {
    i < (s.len() - 2) &&
        !is_cytosine(s[i]) &&
        !is_cytosine(s[i + 1]) &&
        is_guanine(s[i + 2])
}

pub fn is_c_at_g(s: &Vec<u8>, i: usize) -> bool {
    i < (s.len() - 2) &&
        is_cytosine(s[i]) &&
        !is_cytosine(s[i + 1]) &&
        !is_guanine(s[i + 1]) &&
        is_guanine(s[i + 2])
}

pub fn is_chh(s: &Vec<u8>, i: usize) -> bool {
    i < (s.len() - 2) &&
        is_cytosine(s[i]) &&
        !is_guanine(s[i + 1]) &&
        !is_guanine(s[i + 2])
}

pub fn is_cpg(s: &Vec<u8>, i: usize) -> bool {
    i < (s.len() - 1) && is_cytosine(s[i]) && is_guanine(s[i+1])
}


pub fn get_chrom_sizes(
    infile_name: &String
) -> Result<Vec<(String, u64)>, Box<dyn Error>> {

    let bam = Reader::from_path(&infile_name)?;
    let header = bam::Header::from_template(bam.header());

    let mut chrom_sizes: Vec<(String, u64)> = Vec::new();
    for (key, records) in header.to_hashmap() {
        if key == "SQ" {
            for record in records {
                let chrom_name = match record.get("SN") {
                    Some(cn) => cn.clone(),
                    None => return Err("failed to parse header".into()),
                };
                let chrom_size = match record.get("LN") {
                    Some(cs) => cs.parse::<u64>().unwrap(),
                    None => return Err("failed to parse header".into()),
                };
                chrom_sizes.push((chrom_name, chrom_size));
            }
        }
    }
    Ok(chrom_sizes)
}


fn get_context_tag(s: &Vec<u8>, p: usize) -> Vec<u8> {
    if is_cytosine(s[p]) {
        if is_cpg(s, p) {
            return b"CpG".to_vec();
        }
        else if is_chh(s, p) {
            return b"CHH".to_vec();
        }
        else if is_c_at_g(s, p) {
            return b"CXG".to_vec();
        }
        else {
            return b"CCG".to_vec();
        }
    }
    if is_guanine(s[p]) {
        if is_cpg(s, p - 1) {
            return b"CpG".to_vec();
        }
        else if is_ddg(s, p - 2) {
            return b"CHH".to_vec();
        }
        else if is_c_at_g(s, p - 2) {
            return b"CXG".to_vec();
        }
        else {
            return b"CCG".to_vec();
        }
    }
    b"N".to_vec()
}


#[derive(Clone,Default)]
#[allow(non_snake_case)]
struct CountSet {
    pA: u32,
    pC: u32,
    pG: u32,
    pT: u32,

    nA: u32,
    nC: u32,
    nG: u32,
    nT: u32,

    N: u32,
}

impl CountSet {
    pub fn count_pos(&mut self, x: u8) {
        match x {
            b'A' => self.pA += 1,
            b'C' => self.pC += 1,
            b'G' => self.pG += 1,
            b'T' => self.pT += 1,
            _ => self.N += 1
        }
    }
    pub fn count_neg(&mut self, x: u8) {
        match x {
            b'A' => self.nA += 1,
            b'C' => self.nC += 1,
            b'G' => self.nG += 1,
            b'T' => self.nT += 1,
            _ => self.N += 1
        }
    }
    pub fn pos_total(&self) -> u32 {
        self.pA + self.pC + self.pG + self.pT
    }
    pub fn neg_total(&self) -> u32 {
        self.nA + self.nC + self.nG + self.nT
    }
    pub fn unconverted_cytosine(&self) -> u32 { self.pC }
    pub fn converted_cytosine(&self) -> u32 { self.pT }
    pub fn unconverted_guanine(&self) -> u32 { self.nC }
    pub fn converted_guanine(&self) -> u32 { self.nT }

    /* "has_mutated" looks on the opposite strand to see if the apparent
     * conversion from C->T was actually already in the DNA because of a
     * mutation or SNP.
     */
    pub fn has_mutated(&self, base: u8) -> bool {
        const MUT_FRACTION: f64 = 0.5;
        let x = is_cytosine(base);
        (x  && (self.nG as f64) < MUT_FRACTION*(self.neg_total() as f64)) ||
        (!x && (self.pG as f64) < MUT_FRACTION*(self.pos_total() as f64))
    }
}

impl std::fmt::Display for CountSet {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
               self.pA, self.pC, self.pG, self.pT,
               self.nA, self.nC, self.nG, self.nT, self.N)
    }
}


fn write_output<W: Write>(
    out: &mut BufWriter<W>,
    chrom_name: &String,
    chrom: &Vec<u8>,
    counts: &Vec<CountSet>,
) {

    let mut the_site = MSite::new();
    the_site.chrom = chrom_name.as_bytes().to_vec();

    for i in 0..counts.len() {
        let base = chrom[i];
        if is_cytosine(base) || is_guanine(base) {
            the_site.pos = i as u64;
            let is_c = is_cytosine(base);
            the_site.strand = match is_c { true => '+', _ => '-' };
            let unconverted = match is_c {
                true => counts[i].unconverted_cytosine(),
                _ => counts[i].unconverted_guanine()
            } as f64;
            let converted = match is_c {
                true => counts[i].converted_cytosine(),
                _ => counts[i].converted_guanine()
            } as f64;
            the_site.n_reads = (unconverted + converted) as u64;

            the_site.meth = match the_site.n_reads {
                0 => 0.0,
                _ => unconverted/(converted + unconverted),
            };
            the_site.context = get_context_tag(chrom, i);
            if counts[i].has_mutated(base) {
                the_site.context.push(b'x');
            }

            writeln!(out, "{}", the_site).unwrap_or_else(|err| {
                eprintln!("failed writing site: {err}");
                std::process::exit(1);
            });
        }
    }
}


pub fn process_reads(
    verbose: bool,
    bam_filename: String,
    out_filename: String,
    names_in: Vec<String>,
    chroms: Vec<Vec<u8>>,
) {

    if verbose {
        eprintln!("opening BAM file: {}", bam_filename);
    }

    let chrom_info = get_chrom_sizes(&bam_filename)
        .unwrap_or_else(|err| {
            eprintln!("{err}");
            std::process::exit(1);
        });

    let (names, sizes): (Vec<String>, Vec<u64>) =
        chrom_info.into_iter().unzip();

    let mut chrom_lookup: HashMap<String, usize> = HashMap::new();
    for i in 0..names_in.len() {
        chrom_lookup.insert(names_in[i].clone(), i);
    }

    for i in &names {
        if !chrom_lookup.contains_key(i) {
            eprintln!("chrom not found: {}", i);
            std::process::exit(1);
        }
    }

    let valid_names: HashSet<String> = names.iter().cloned().collect();

    let mut bam = Reader::from_path(&bam_filename).unwrap_or_else(|err| {
        eprintln!("{err}");
        std::process::exit(1);
    });

    use std::fs::File;
    let out = File::create(out_filename).unwrap_or_else(|err| {
        eprintln!("{err}");
        std::process::exit(1);
    });
    let mut out = BufWriter::new(out);

    let mut counts: Vec<CountSet> = Vec::new();
    let mut chroms_seen: HashSet<i32> = HashSet::new();

    let mut tid: i32 = -1; // i32 is used in htslib for tid
    let mut rname = String::new();

    // ADS: AbstractInterval gives us aln.contig()
    use bio_types::genome::AbstractInterval;

    for r in bam.records() {
        let aln = r.ok().expect("error reading BAM record");

        // if chrom changes, output previous results, get new one
        if aln.tid() != tid {
            if aln.tid() != tid {
                if !valid_names.contains(aln.contig()) {
                    eprintln!("chrom seq not found: {}", aln.contig());
                    std::process::exit(1);
                }
                if chroms_seen.contains(&aln.tid()) {
                    eprintln!("reads from same chrom not consecutive");
                    std::process::exit(1);
                }
                chroms_seen.insert(aln.tid());
            }
            if !counts.is_empty() {
                let idx = chrom_lookup[&rname] as usize;
                write_output(&mut out, &rname, &chroms[idx], &counts);
            }

            // update current target id and reference name
            tid = aln.tid();
            rname = aln.contig().to_string();
            counts.clear();
            counts.resize_with(sizes[tid as usize] as usize, Default::default);

            assert!(rname == names[tid as usize]);

            if verbose {
                eprintln!("PROCESSING:\t{}", names[tid as usize]);
            }
        }

        // process current read, depending on strand
        use rust_htslib::bam::ext::BamRecordExtensions;
        if aln.is_reverse() {
            let l = aln.seq_len() as i64 - 1;
            for i in aln.aligned_pairs() {
                counts[i[1] as usize].count_neg(aln.seq()[(l - i[0]) as usize]);
            }
        }
        else {
            for i in aln.aligned_pairs() {
                counts[i[1] as usize].count_pos(aln.seq()[i[0] as usize]);
            }
        }
    }
    if !counts.is_empty() {
        let idx = chrom_lookup.get(&rname).unwrap();
        write_output(&mut out, &rname, &chroms[*idx], &counts);
    }
}
