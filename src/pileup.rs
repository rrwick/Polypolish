// Copyright 2021 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Polypolish

//This file is part of Polypolish. Polypolish is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Polypolish
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Polypolish. If not, see <http://www.gnu.org/licenses/>.

use std::collections::HashMap;
use std::cmp;
use math::round::half_to_even;
use crate::alignment::Alignment;


#[derive(Debug)]
pub struct PileupBase {
    original: char,
    depth: f64,
    counts: HashMap<String, usize>,
}

impl PileupBase {
    fn new(original: char) -> PileupBase {
        PileupBase {
            original: original.to_ascii_uppercase(),
            depth: 0.0,
            counts: HashMap::new(),
        }
    }
    pub fn add_seq(&mut self, seq: &str, depth_contribution: f64) {
        *self.counts.entry(seq.to_string()).or_insert(0) += 1;
        self.depth += depth_contribution;
    }
    // TODO: method to get valid choices
    // TODO: method to get output sequence

    pub fn get_output_seq(&self, min_depth: usize, min_fraction: f64) -> String {
        let valid_seqs = self.get_valid_seqs(min_depth, min_fraction);
        if valid_seqs.len() == 1 {
            valid_seqs[0].clone()
        } else {
            self.original.to_string()
        }
    }

    fn get_valid_seqs(&self, min_depth: usize, min_fraction: f64) -> Vec<String> {
        // Use round-to-even logic so the behaviour matches previous Python implementation.
        let threshold = cmp::max(min_depth,
                                 half_to_even(self.depth * min_fraction, 0) as usize);
        let mut valid_seqs = Vec::new();
        for (seq, count) in &self.counts {
            if count >= &threshold {
                valid_seqs.push(seq.clone());
            }
        }
        valid_seqs
    }
}


#[derive(Debug)]
pub struct Pileup {
    pub bases: Vec<PileupBase>,
}

impl Pileup {
    pub fn new(seq: &str) -> Pileup {
        let mut bases = Vec::new();
        for b in seq.chars() {
            bases.push(PileupBase::new(b));
        }

        Pileup {
            bases: bases,
        }
    }

    pub fn add_alignment(&mut self, alignment: &Alignment, depth_contribution: f64) {
        let read_bases_per_pos = alignment.get_read_bases_for_each_target_base();
        let mut i = alignment.ref_start;
        for b in read_bases_per_pos {
            self.bases[i].add_seq(&b, depth_contribution);
            i += 1;
        }
    }
}