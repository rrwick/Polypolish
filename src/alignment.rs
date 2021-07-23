// Copyright 2021 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Polypolish

//This file is part of Polypolish. Polypolish is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Polypolish
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Polypolish. If not, see <http://www.gnu.org/licenses/>.

use crate::pileup::Pileup;
use crate::misc::{quit_with_error, reverse_complement};

use std::path::PathBuf;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{prelude::*, BufReader};
use std::result::Result;

use regex::Regex;
use lazy_static::lazy_static;


lazy_static! {
    static ref RE: Regex = Regex::new(r"\d+[MIDNSHP=X]").unwrap();
}


#[derive(Debug)]
pub struct Alignment {
    read_name: String,
    ref_name: String,
    sam_flags: u32,
    pub ref_start: usize,
    cigar: String,
    expanded_cigar: String,
    pub read_seq: String,
    mismatches: u32,
}

impl Alignment {
    fn new(sam_line: &str) -> Result<Alignment, &str> {
        let parts = sam_line.split('\t').collect::<Vec<&str>>();
        if parts.len() < 11 {
            return Err("too few columns");
        }

        let read_name = parts[0];
        let sam_flags = parts[1].parse::<u32>().unwrap();
        let ref_name = parts[2];
        let ref_start = parts[3].parse::<usize>().unwrap() - 1;
        let cigar = parts[5];
        let read_seq = parts[9];

        let mut mismatches = u32::MAX;
        for p in &parts[11..] {
            if p.starts_with("NM:i:") {
                let nm = p[5..].to_string();
                mismatches = nm.parse::<u32>().unwrap();
            }
        }
        if mismatches == u32::MAX && sam_flags & 4 == 0 {
            return Err("missing NM tag");
        }

        Ok(Alignment {
            read_name: read_name.to_string(),
            ref_name: ref_name.to_string(),
            sam_flags: sam_flags,
            ref_start: ref_start,
            cigar: cigar.to_string(),
            expanded_cigar: get_expanded_cigar(&cigar, read_seq.len()),
            read_seq: read_seq.to_ascii_uppercase(),
            mismatches: mismatches,
        })
    }

    fn is_aligned(&self) -> bool {
        (self.sam_flags & 4) == 0
    }

    fn get_strand(&self) -> i8 {
        if self.is_on_forward_strand() {
            1
        } else {
            -1
        }
    }

    fn is_on_forward_strand(&self) -> bool {
        (self.sam_flags & 16) == 0
    }

    fn starts_and_ends_with_match(&self) -> bool {
        self.expanded_cigar.chars().next().unwrap() == 'M' &&
            self.expanded_cigar.chars().last().unwrap() == 'M'
    }

    fn add_read_seq(&mut self, read_seq: &str, strand: i8) {
        if self.get_strand() == strand {
            self.read_seq = read_seq.to_string();
        } else {
            self.read_seq = reverse_complement(read_seq);
        }
    }

    /// This function returns a vector giving the read base(s) for each position of the target
    /// sequence. Instead of returning these as strings (which would involve a lot of allocation
    /// of new strings to memory which is slow), it returns them as start/end indices of the read
    /// sequence. Most values will have an end one more than the start (e.g. 5,6) indicating a
    /// single base. However, insertions can lead to bigger ranges (e.g. 5,7) and deletions to
    /// zero-length ranges (e.g. 5,5).
    pub fn get_read_bases_for_each_target_base(&self) -> Vec<(usize, usize)> {
        let mut i = 0;
        let mut read_bases = Vec::with_capacity(self.expanded_cigar.len());
        for c in self.expanded_cigar.chars() {
            if c == 'M' {
                read_bases.push((i, i+1));
                i += 1;
            } else if c == 'I' {
                read_bases.last_mut().unwrap().1 = i+1;
                i += 1;
            } else if c == 'D' {
                read_bases.push((i, i));
            } else {
                quit_with_error(&format!("unexpected character in CIGAR string for read {}: {:?}",
                                         self.read_name, self.cigar));
            }
        }
        if i != self.read_seq.len() {
            quit_with_error(&format!("CIGAR string for read {} does not match read sequence",
                                     self.read_name));
        }
        trim_bases_for_homopolymers(&mut read_bases, &self.read_seq);
        read_bases
    }
}


pub fn process_sam(filename: &PathBuf, pileups: &mut HashMap<String, Pileup>,
                   max_errors: u32) -> (usize, usize, usize) {
    let result = add_to_pileup(filename, pileups, max_errors);
    match result {
        Ok((_,_,_)) => (),
        Err(_)      => quit_with_error(&format!("unable to load alignments from {:?}", filename)),
    }
    result.unwrap()
}


pub fn add_to_pileup(filename: &PathBuf, pileups: &mut HashMap<String, Pileup>,
                     max_errors: u32) -> io::Result<(usize, usize, usize)> {
    let file = File::open(&filename)?;
    let reader = BufReader::new(file);

    let mut current_read_name = String::new();
    let mut current_read_alignments = Vec::new();

    let mut line_count: usize = 0;
    let mut alignment_count: usize = 0;
    let mut used_count: usize = 0;
    let mut read_count: usize = 0;

    for line in reader.lines() {
        line_count += 1;
        let sam_line = line?;
        if sam_line.len() == 0 {continue;}
        if sam_line.starts_with('@') {continue;}

        let alignment_result = Alignment::new(&sam_line);
        match alignment_result {
            Ok(_)  => (),
            Err(e) => quit_with_error(&format!("{} in {:?} (line {})", e, filename, line_count)),
        }
        let alignment = alignment_result.unwrap();
        if !alignment.is_aligned() {continue;}

        alignment_count += 1;
        let read_name = alignment.read_name.clone();

        if current_read_name.is_empty() || current_read_name == alignment.read_name {
            current_read_alignments.push(alignment);
        } else {
            used_count += process_one_read(current_read_alignments, pileups, max_errors);
            read_count += 1;
            current_read_alignments = Vec::new();
            current_read_alignments.push(alignment);
        }
        current_read_name = read_name;
    }
    used_count += process_one_read(current_read_alignments, pileups, max_errors);
    read_count += 1;

    if alignment_count == 0 {
        quit_with_error(&format!("no alignments in {:?}", filename))
    }
    Ok((alignment_count, used_count, read_count))
}


fn process_one_read(alignments: Vec<Alignment>, pileups: &mut HashMap<String, Pileup>,
                    max_errors: u32) -> usize {
    let (read_seq, strand) = get_read_seq_from_alignments(&alignments);

    let mut good_alignments = Vec::new();
    for a in alignments {
        if a.starts_and_ends_with_match() && a.mismatches <= max_errors {
            good_alignments.push(a);
        }
    }
    let depth_contribution = 1.0 / good_alignments.len() as f64;

    for a in &mut good_alignments {
        let needs_length = a.read_seq == "*";
        if needs_length {
            a.add_read_seq(&read_seq, strand);
        }
    }

    for a in &good_alignments {
        if !pileups.contains_key(&a.ref_name) {
            quit_with_error(&format!("query name {} in SAM but not in assembly", a.ref_name))
        }
        let pileup = pileups.get_mut(&a.ref_name).unwrap();
        pileup.add_alignment(a, depth_contribution);
    }
    good_alignments.len()
}


/// This function takes a vector of all the alignments for one read. At least one of these
/// alignments should have the read seq included (i.e. not just "*"). This function will return
/// that sequence and its strand.
fn get_read_seq_from_alignments(alignments: &Vec<Alignment>) -> (String, i8) {
    for a in alignments {
        if a.read_seq == "*" {
            continue;
        } else {
            return (a.read_seq.clone(), a.get_strand());
        }
    }
    let read_name = &alignments.first().unwrap().read_name;
    quit_with_error(&format!("no alignments for read {} contain sequence", read_name));
    ("".to_string(), 0)  // never reached
}


fn get_expanded_cigar(cigar: &str, read_seq_len: usize) -> String {
    let mut expanded_cigar = String::with_capacity(read_seq_len);
    for m in RE.find_iter(cigar) {
        let num: u32 = cigar[m.start()..m.end()-1].parse().unwrap();
        let letter = &cigar[m.end()-1..m.end()];
        for _ in 0..num {
            expanded_cigar.push_str(letter);
        }
    }
    expanded_cigar
}


/// Alignments that end in a homopolymer can cause trouble, as they can align cleanly
/// (without an indel) even when an indel is needed.
///
/// For example, an alignment should look like this:
///   read: ... T G A G T A C AG G G G G A A G T
///   ref:  ... T G A G T A C A  G G G G A A G T C C A G T ...
///
/// But if the read ends in the homopolymer, it could look like this:
///   read: ... T G A G T A C A G G
///   ref:  ... T G A G T A C A G G G G A A G T C C A G T ...
///
/// Which results in a clean alignment on the 'A' that should be 'AG'. To avoid this, we
/// trim off the last couple unique bases of the alignment, so the example becomes:
///   read: ... T G A G T A C
///   ref:  ... T G A G T A C A G G G G A A G T C C A G T ...
fn trim_bases_for_homopolymers(read_bases: &mut Vec<(usize, usize)>, read_seq: &str) {
    let (last_start, last_end) = *read_bases.last().unwrap();
    let last_base = &read_seq[last_start..last_end];
    while read_bases.len() > 0 {
        let (current_last_start, current_last_end) = *read_bases.last().unwrap();
        let current_last_base = &read_seq[current_last_start..current_last_end];
        if current_last_base != last_base {
            break;
        }
        read_bases.pop();
    }
    if read_bases.len() > 0 {
        read_bases.pop();
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_expanded_cigar() {
        assert_eq!(get_expanded_cigar("10M", 10), "MMMMMMMMMM");
        assert_eq!(get_expanded_cigar("3M1I7M", 11), "MMMIMMMMMMM");
        assert_eq!(get_expanded_cigar("5M2D4M", 9), "MMMMMDDMMMM");
    }
}
