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
use crate::misc::quit_with_error;
use std::path::PathBuf;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{prelude::*, BufReader};
use std::result::Result;


#[derive(Debug)]
pub struct Alignment {
    read_name: String,
    ref_name: String,
    sam_flags: i32,
    ref_start: i64,
    cigar: String,
    read_seq: String,
}


impl Alignment {
    fn new(sam_line: &str) -> Result<Alignment, &str> {
        let parts = sam_line.split('\t').collect::<Vec<&str>>();
        if parts.len() < 11 {
            return Err("bad SAM format");
        }

        let read_name = parts[0];
        let sam_flags = parts[1].parse::<i32>().unwrap();
        let ref_name = parts[2];
        let ref_start = parts[3].parse::<i64>().unwrap() - 1;
        let cigar = parts[5];
        let read_seq = parts[9];

        Ok(Alignment {
            read_name: read_name.to_string(),
            ref_name: ref_name.to_string(),
            sam_flags: sam_flags,
            ref_start: ref_start,
            cigar: cigar.to_string(),
            read_seq: read_seq.to_string(),
        })
    }
}


pub fn process_sam(filename: &PathBuf, pileups: &mut HashMap<String, Pileup>) {
    let result = add_to_pileup(filename, pileups);
    match result {
        Ok((_, _)) => ( ),
        Err(_) => quit_with_error(&format!("unable to load alignments from {:?}", filename)),
    }
    let (line_count, used_count) = result.unwrap();
    eprintln!("{} alignments loaded from {:?} ({} used)", line_count, filename, used_count)
}


pub fn add_to_pileup(filename: &PathBuf,
                     pileups: &mut HashMap<String, Pileup>) -> io::Result<(usize, usize)> {
    let file = File::open(&filename)?;
    let reader = BufReader::new(file);

    let mut current_read_name = String::new();
    let mut current_read_alignments = Vec::new();

    let mut line_count: usize = 0;
    let mut used_count: usize = 0;

    for line in reader.lines() {
        let sam_line = line?;
        if sam_line.len() == 0 {continue;}
        if sam_line.starts_with('@') {continue;}

        let alignment_result = Alignment::new(&sam_line);
        match alignment_result {
            Ok(_) => ( ),
            Err(_) => quit_with_error(&format!("{:?} is not correctly formatted", filename)),
        }
        let alignment = alignment_result.unwrap();

        line_count += 1;
        let read_name = alignment.read_name.clone();

        if current_read_name.is_empty() || current_read_name == alignment.read_name {
            current_read_alignments.push(alignment);
        } else {
            used_count += process_one_read(current_read_alignments, pileups);
            current_read_alignments = Vec::new();
            current_read_alignments.push(alignment);
        }
        current_read_name = read_name;
    }
    used_count += process_one_read(current_read_alignments, pileups);

    if line_count == 0 {
        quit_with_error(&format!("no alignments in {:?}", filename))
    }
    Ok((line_count, used_count))
}


fn process_one_read(alignments: Vec<Alignment>, pileups: &mut HashMap<String, Pileup>) -> usize {
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    // TODO
    for a in alignments {
        eprintln!("{:?}", a);
    }
    eprintln!();

    1  //TEMP
}



// TODO: define reverse complement function
