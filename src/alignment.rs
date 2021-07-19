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
use std::path::PathBuf;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{prelude::*, BufReader};
use crate::load_assembly;


#[derive(Debug)]
pub struct Alignment {
    read_name: String,
    ref_name: String,
    sam_flags: i32,
    ref_start: i64,
    cigar: String,
    read_seq: String,
}

struct AlignmentError;

impl Alignment {
    fn new(sam_line: &str) -> Alignment {
        let parts = sam_line.split('\t').collect::<Vec<&str>>();
        if parts.len() < 11 {
            panic!();
        }

        let read_name = parts[0];
        let sam_flags = parts[1].parse::<i32>().unwrap();
        let ref_name = parts[2];
        let ref_start = parts[3].parse::<i64>().unwrap() - 1;
        let cigar = parts[5];
        let read_seq = parts[9];

        Alignment {
            read_name: read_name.to_string(),
            ref_name: ref_name.to_string(),
            sam_flags: sam_flags,
            ref_start: ref_start,
            cigar: cigar.to_string(),
            read_seq: read_seq.to_string(),
        }
    }
}


// TODO: define reverse complement function



pub fn process_sam(filename: &PathBuf, pileups: &HashMap<String, Pileup>) -> io::Result<()> {
    let file = File::open(&filename)?;
    let reader = BufReader::new(file);

    let mut current_read_name = String::new();
    let mut current_read_alignments = Vec::new();

    for line in reader.lines() {
        let sam_line = line?;
        if sam_line.len() == 0 {continue;}
        if sam_line.starts_with('@') {continue;}
        let alignment = Alignment::new(&sam_line);
        let read_name = alignment.read_name.clone();

        if current_read_name.is_empty() || current_read_name == alignment.read_name {
            current_read_alignments.push(alignment);
        } else {
            process_one_read(current_read_alignments, &pileups);
            current_read_alignments = Vec::new();
            current_read_alignments.push(alignment);
        }
        current_read_name = read_name;
    }
    process_one_read(current_read_alignments, &pileups);

    Ok(())
}


fn process_one_read(alignments: Vec<Alignment>, pileups: &HashMap<String, Pileup>) {
    for a in alignments {
        eprintln!("{:?}", a);
    }
    eprintln!();
}
