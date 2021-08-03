// Copyright 2021 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Polypolish

//This file is part of Polypolish. Polypolish is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Polypolish
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Polypolish. If not, see <http://www.gnu.org/licenses/>.

use crate::alignment::Alignment;
use crate::misc::bankers_rounding;

use std::collections::HashMap;


pub enum BaseStatus {
    NoValidOptions,       // no sequences pass the threshold (not changed)
    MultipleValidOptions, // multiple sequences pass the threshold (not changed)
    OriginalBaseKept,     // one sequence passes the threshold and it matches the original base
    Changed,              // one sequence passes the threshold and it differs from the original base
}


#[derive(Debug)]
pub struct PileupBase {
    original: char,
    pub depth: f64,

    // A, C, G and T are the most common sequences, so we count them with integers (fast):
    count_a: u32,
    count_c: u32,
    count_g: u32,
    count_t: u32,

    // Everything else will be counted in a HashMap (slower but can handle any sequence):
    counts: HashMap<String, u32>,
}

impl PileupBase {
    fn new(original: char) -> PileupBase {
        PileupBase {
            original: original,
            depth: 0.0,
            count_a: 0,
            count_c: 0,
            count_g: 0,
            count_t: 0,
            counts: HashMap::new(),
        }
    }

    pub fn add_seq(&mut self, seq: &str, depth_contribution: f64) {
        match seq {
            "A" => {self.count_a += 1},
            "C" => {self.count_c += 1},
            "G" => {self.count_g += 1},
            "T" => {self.count_t += 1},
             _  => {*self.counts.entry(seq.to_string()).or_insert(0) += 1},
        }
        self.depth += depth_contribution;
    }

    pub fn get_polished_seq(&self, min_depth: u32, min_fraction: f64,
                            build_debug_line: bool) -> (String, BaseStatus, String) {
        let original = self.original.to_string();
        let threshold = std::cmp::max(min_depth, bankers_rounding(self.depth * min_fraction));
        let mut valid_seqs = Vec::new();
        if self.count_a >= threshold {valid_seqs.push("A".to_string());}
        if self.count_c >= threshold {valid_seqs.push("C".to_string());}
        if self.count_g >= threshold {valid_seqs.push("G".to_string());}
        if self.count_t >= threshold {valid_seqs.push("T".to_string());}
        for (seq, count) in &self.counts {
            if count >= &threshold {
                valid_seqs.push(seq.clone());
            }
        }

        let mut new_base = original.clone();
        let mut status = BaseStatus::OriginalBaseKept;

        if self.depth < min_depth as f64 {
            status = BaseStatus::NoValidOptions;
        } else if valid_seqs.len() == 1 {
            new_base = valid_seqs[0].clone();
            if new_base != original {
                status = BaseStatus::Changed;
            }
        } else if valid_seqs.len() == 0 {
            status = BaseStatus::NoValidOptions;
        } else {
            status = BaseStatus::MultipleValidOptions;
        }

        let debug_line = self.get_debug_line(build_debug_line, threshold, &status, &new_base);
        (new_base, status, debug_line)
    }

    /// Returns the sequence counts in string form (used in the debug output).
    fn get_count_str(&self) -> String {
        let mut counts = Vec::new();
        if self.count_a > 0 {counts.push(format!("Ax{}", self.count_a));}
        if self.count_c > 0 {counts.push(format!("Cx{}", self.count_c));}
        if self.count_g > 0 {counts.push(format!("Gx{}", self.count_g));}
        if self.count_t > 0 {counts.push(format!("Tx{}", self.count_t));}
        for (seq, count) in &self.counts {
            counts.push(format!("{}x{}", seq, count));
        }
        counts.sort();
        counts.join(",")
    }

    fn get_debug_line(&self, build_debug_line: bool, threshold: u32, status: &BaseStatus,
                      new_base: &str) -> String {
        if !build_debug_line {
            return String::new();
        }

        let status_str = match status {
            BaseStatus::OriginalBaseKept     => "kept",
            BaseStatus::Changed              => "changed",
            BaseStatus::NoValidOptions       => "none",
            BaseStatus::MultipleValidOptions => "multiple",
        };
        format!("{}\t{:.1}\t{}\t{}\t{}\t{}", self.original, self.depth, threshold,
                self.get_count_str(), status_str, new_base)
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
        let read_bases = alignment.get_read_bases_for_each_target_base();
        let mut i = alignment.ref_start;
        for (start, end) in read_bases {
            if start == end {
                self.bases[i].add_seq("-", depth_contribution);
            } else {
                self.bases[i].add_seq(&alignment.read_seq[start..end], depth_contribution);
            }
            i += 1;
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pileupbase_1() {
        let mut b = PileupBase::new('G');
        b.add_seq("A", 1.0);
        b.add_seq("T", 1.0);
        for _ in 0..50 {b.add_seq("G", 1.0);}
        assert_eq!(b.get_count_str(), "Ax1,Gx50,Tx1");
        let (polished, status, _) = b.get_polished_seq(5, 0.5, false);
        assert_eq!(polished, "G");
        assert!(matches!(status, BaseStatus::OriginalBaseKept));
    }

    #[test]
    fn test_pileupbase_2() {
        let mut b = PileupBase::new('T');
        b.add_seq("C", 1.0);
        for _ in 0..99 {b.add_seq("A", 1.0);}
        assert_eq!(b.get_count_str(), "Ax99,Cx1");
        let (polished, status, _) = b.get_polished_seq(5, 0.5, false);
        assert_eq!(polished, "A");
        assert!(matches!(status, BaseStatus::Changed));
    }

    #[test]
    fn test_pileupbase_3() {
        let mut b = PileupBase::new('A');
        b.add_seq("T", 1.0);
        b.add_seq("C", 1.0);
        b.add_seq("G", 1.0);
        assert_eq!(b.get_count_str(), "Cx1,Gx1,Tx1");
        let (polished, status, _) = b.get_polished_seq(5, 0.5, false);
        assert_eq!(polished, "A");
        assert!(matches!(status, BaseStatus::NoValidOptions));
    }

    #[test]
    fn test_pileupbase_4() {
        let mut b = PileupBase::new('C');
        for _ in 0..123 {b.add_seq("A", 0.1);}
        for _ in 0..321 {b.add_seq("T", 0.1);}
        assert_eq!(b.get_count_str(), "Ax123,Tx321");
        let (polished, status, _) = b.get_polished_seq(5, 0.5, false);
        assert_eq!(polished, "C");
        assert!(matches!(status, BaseStatus::MultipleValidOptions));
    }
}
