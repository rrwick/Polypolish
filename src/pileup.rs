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


#[derive(Debug)]
pub struct PileupBase {
    original: char,
    depth: f64,
    counts: HashMap<String, usize>,
}

impl PileupBase {
    fn new(original: char) -> PileupBase {
        PileupBase {
            original: original,
            depth: 0.0,
            counts: HashMap::new(),
        }
    }
    fn add_seq(&mut self, seq: &str, depth_contribution: f64) {
        *self.counts.entry(seq.to_string()).or_insert(0) += 1;
        self.depth += depth_contribution;
    }
    // TODO: method to get valid choices
    // TODO: method to get output sequence
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

    // TODO: method to add an alignment
}