// Copyright 2021 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Polypolish

// This file is part of Polypolish. Polypolish is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Polypolish
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Polypolish. If not, see <http://www.gnu.org/licenses/>.

use std::path::PathBuf;
use std::collections::{HashMap, HashSet};
use std::time::Instant;
use std::fs::File;
use std::io;
use std::io::{prelude::*, BufReader, BufWriter};
use clap::crate_version;
use num_format::{Locale, ToFormattedString};

use crate::alignment::Alignment;
use crate::log;
use crate::misc::{quit_with_error, format_duration};


pub fn filter(in1: PathBuf, in2: PathBuf, out1: PathBuf, out2: PathBuf,
              orientation: String, low: f64, high: f64) {
    let start_time = Instant::now();
    check_inputs(&in1, &in2, &out1, &out2, low, high);
    starting_message(&in1, &in2, &out1, &out2, &orientation, low, high);
    let (alignments, before_count) = load_alignments(&in1, &in2);
    let (low, high, correct_orientation) = get_insert_size_thresholds(&alignments, &orientation,
                                                                      low, high);
    let after_count = filter_sams(&in1, &in2, &out1, &out2, &alignments, low, high,
                                  correct_orientation);
    finished_message(start_time, before_count, after_count)
}


fn check_inputs(in1: &PathBuf, in2: &PathBuf, out1: &PathBuf, out2: &PathBuf,
                low: f64, high: f64) {
    let mut files = HashSet::new();
    if !files.insert(in1.clone()) || !files.insert(in2.clone()) || 
        !files.insert(out1.clone()) || !files.insert(out2.clone()) {
        quit_with_error("--in1, --in2, --out1 and --out2 must all have unique values");
    }
    if low <= 0.0 || low >= 50.0 {
        quit_with_error("--low must be greater than 0 and less than 50")
    }
    if high <= 50.0 || high >= 100.0 {
        quit_with_error("--high must be greater than 50 and less than 100")
    }
}


fn starting_message(in1: &PathBuf, in2: &PathBuf, out1: &PathBuf, out2: &PathBuf,
                    orientation: &String, low: f64, high: f64) {
    log::section_header("Starting Polypolish filter");
    log::explanation("This runs a pre-processing filter on SAM alignments before they are used to \
                      polish. It looks at each read pair and flags alignments that do not seem to \
                      be part of a concordant pair. This can improve the accuracy Polypolish, \
                      especially near the edges of repeats.");
    eprintln!("Polypolish version: {}", crate_version!());
    eprintln!();
    eprintln!("Input alignments:");
    eprintln!("  {}", in1.display());
    eprintln!("  {}", in2.display());
    eprintln!();
    eprintln!("Output alignments:");
    eprintln!("  {}", out1.display());
    eprintln!("  {}", out2.display());
    eprintln!();
    eprintln!("Settings:");
    eprintln!("  --orientation {}", orientation);
    eprintln!("  --low {}", low);
    eprintln!("  --high {}", high);
    eprintln!();
}


fn finished_message(start_time: Instant, before_count: usize, after_count: usize) {
    log::section_header("Finished!");
    eprintln!("Alignments before filtering: {}", before_count.to_formatted_string(&Locale::en));
    eprintln!("Alignments after filtering:  {}", after_count.to_formatted_string(&Locale::en));
    eprintln!();
    eprintln!("Time to run: {}", format_duration(start_time.elapsed()));
    eprintln!();
}


fn load_alignments(sam_1: &PathBuf, sam_2: &PathBuf) -> (HashMap<String, Vec<Alignment>>, usize) {
    log::section_header("Loading alignments");
    let mut alignments = HashMap::new();
    let result_1 = load_alignments_one_file(sam_1, &mut alignments, "_1");
    match result_1 {
        Ok(()) => (),
        Err(_) => quit_with_error(&format!("unable to load alignments from {:?}", sam_1)),
    }
    let result_2 = load_alignments_one_file(sam_2, &mut alignments, "_2");
    match result_2 {
        Ok(()) => (),
        Err(_) => quit_with_error(&format!("unable to load alignments from {:?}", sam_2)),
    }
    eprintln!();
    let count = alignments.values().map(|v| v.len()).sum();
    (alignments, count)
}


fn load_alignments_one_file(sam_filename: &PathBuf,
                            alignments: &mut HashMap<String, Vec<Alignment>>,
                            read_name_suffix: &str) -> io::Result<()> {
    eprint!("{}: ", sam_filename.display());
    let sam_file = File::open(sam_filename)?;
    let reader = BufReader::new(sam_file);
    let mut alignment_count = 0;
    let mut read_names = HashSet::new();
    let mut line_count: usize = 0;
    for line in reader.lines() {
        line_count += 1;
        let sam_line = line?;
        if sam_line.starts_with('@') {
            continue;
        }
        let alignment_result = Alignment::new_quick(&sam_line);
        match alignment_result {
            Ok(_)  => (),
            Err(e) => quit_with_error(&format!("{} in {:?} (line {})",
                                               e, sam_filename, line_count)),
        }
        let mut alignment = alignment_result.unwrap();
        if !alignment.is_aligned() {continue;}
        alignment.read_name.push_str(read_name_suffix);
        read_names.insert(alignment.read_name.clone());
        alignments.entry(alignment.read_name.clone()).or_insert_with(Vec::new).push(alignment);
        alignment_count += 1;
    }
    eprintln!("{} alignments from {} reads",
              alignment_count.to_formatted_string(&Locale::en),
              read_names.len().to_formatted_string(&Locale::en));
    if alignments.is_empty() {
        quit_with_error(&format!("no alignments found in {:?}", sam_filename));
    }
    Ok(())
}


fn get_insert_size_thresholds(alignments: &HashMap<String, Vec<Alignment>>,
                              correct_orientation: &String,
                              low_percentile: f64, high_percentile: f64) -> (u32, u32, String) {
    log::section_header("Finding insert size thresholds");
    log::explanation("Read pairs with exactly one alignment per read are used to determine the \
                      orientation and insert size thresholds for the read set.");
    let mut insert_sizes: HashMap<String, Vec<u32>> = HashMap::new();
    for (name_1, alignments_1) in alignments {
        if !name_1.ends_with("_1") || alignments_1.len() != 1 {
            continue;
        }
        let name_2 = format!("{}_2", &name_1[..name_1.len() - 2]);
        if let Some(alignments_2) = alignments.get(&name_2) {
            if alignments_2.len() == 1 && alignments_1[0].ref_name == alignments_2[0].ref_name {
                let orientation = get_orientation(&alignments_1[0], &alignments_2[0]);
                let insert_size = get_insert_size(&alignments_1[0], &alignments_2[0]);
                insert_sizes.entry(orientation).or_default().push(insert_size);
            }
        }
    }
    if insert_sizes.is_empty() {
        quit_with_error("no one-alignment-per-read pairs available to determine orientation and \
                         insert size thresholds");
    }

    let correct_orientation = determine_correct_orientation(correct_orientation, &insert_sizes);
    let mut sizes = insert_sizes.remove(&correct_orientation).unwrap_or_else(Vec::new);
    if sizes.is_empty() {
        quit_with_error("no read pairs available to determine insert size thresholds");
    }
    sizes.sort_unstable();
    let low_threshold = get_percentile(&sizes, low_percentile);
    let high_threshold = get_percentile(&sizes, high_percentile);
    eprintln!("Low threshold:  {} ({})", low_threshold, get_percentile_name(low_percentile));
    eprintln!("High threshold: {} ({})", high_threshold, get_percentile_name(high_percentile));
    eprintln!();

    (low_threshold, high_threshold, correct_orientation)
}


fn get_orientation(a_1: &Alignment, a_2: &Alignment) -> String {
    let strand_1 = if a_1.is_on_forward_strand() { 'f' } else { 'r' };
    let strand_2 = if a_2.is_on_forward_strand() { 'f' } else { 'r' };

    // Get the read start positions, which is the ref-end position if on the negative strand.
    let a_1_pos = if a_1.is_on_forward_strand() { a_1.ref_start } else { a_1.get_ref_end() };
    let a_2_pos = if a_2.is_on_forward_strand() { a_2.ref_start } else { a_2.get_ref_end() };

    match (strand_1, strand_2) {
        ('f', 'r') | ('r', 'f') => {
            if a_1_pos < a_2_pos {
                format!("{}{}", strand_1, strand_2)
            } else {
                format!("{}{}", strand_2, strand_1)
            }
        },
        ('f', 'f') => { if a_1_pos < a_2_pos { "ff".to_string() } else { "rr".to_string() } },
        ('r', 'r') => { if a_2_pos < a_1_pos { "ff".to_string() } else { "rr".to_string() } },
        _ => unreachable!()
    }
}


fn get_insert_size(alignment_1: &Alignment, alignment_2: &Alignment) -> u32 {
    let positions = [alignment_1.ref_start, alignment_1.get_ref_end(),
                     alignment_2.ref_start, alignment_2.get_ref_end()];
    let insert_start = positions.iter().min().cloned().unwrap_or_default();
    let insert_end = positions.iter().max().cloned().unwrap_or_default();
    (insert_end - insert_start) as u32
}


fn determine_correct_orientation(correct_orientation: &str,
                                 insert_sizes: &HashMap<String, Vec<u32>>) -> String {
    for orientation in ["fr", "rf", "ff", "rr"].iter() {
        let count = insert_sizes.get(*orientation).map_or(0, |v| v.len());
        eprintln!("{}: {} pairs", orientation, count.to_formatted_string(&Locale::en));
    }
    if correct_orientation == "auto" {
        let auto_orientation = auto_determine_orientation(insert_sizes);
        eprintln!("\nAutomatically determined correct orientation: {}\n", auto_orientation);
        auto_orientation
    } else {
        eprintln!("\nUser-specified correct orientation: {}\n", correct_orientation);
        correct_orientation.to_string()
    }
}


fn auto_determine_orientation(insert_sizes: &HashMap<String, Vec<u32>>) -> String {
    let max_count = insert_sizes.values().map(|v| v.len()).max().unwrap_or(0);
    let orientations: Vec<&str> = ["fr", "rf", "ff", "rr"].iter()
        .filter(|&&o| insert_sizes.get(o).map_or(0, |v| v.len()) == max_count).cloned().collect();
    if orientations.len() != 1 {
        quit_with_error("could not automatically determine read pair orientation");
    }
    orientations[0].to_string()
}


fn get_percentile(sorted_list: &[u32], percentile: f64) -> u32 {
    // Returns a percentile of a list of numbers. Assumes the list has already been sorted.
    // Implements the nearest rank method:
    // https://en.wikipedia.org/wiki/Percentile#The_nearest-rank_method
    if sorted_list.is_empty() {
        return 0;
    }
    let fraction = percentile / 100.0;
    let rank = ((fraction * sorted_list.len() as f64).ceil() as usize).max(1);
    sorted_list.get(rank - 1).cloned().unwrap_or(0)
}


fn get_percentile_name(p: f64) -> String {
    let p_str = p.to_string();
    match p_str.as_str() {
        _ if p_str.ends_with("1") && p != 11.0 => format!("{}st percentile", p),
        _ if p_str.ends_with("2") && p != 12.0 => format!("{}nd percentile", p),
        _ if p_str.ends_with("3") && p != 13.0 => format!("{}rd percentile", p),
        _ => format!("{}th percentile", p),
    }
}


fn filter_sams(in1: &PathBuf, in2: &PathBuf, out1: &PathBuf, out2: &PathBuf,
               alignments: &HashMap<String, Vec<Alignment>>, low: u32, high: u32,
               correct_orientation: String) -> usize {
    log::section_header("Filtering SAM files");
    log::explanation("Read alignments that are part of a good pair (correct orientation and \
                      insert size) pass the filter and are written unaltered to the output file. \
                      Read alignments which are not part of good pair are written to the output \
                      file with a \"ZP:Z:fail\" tag so Polypolish will not use them.");
    let mut after_count = 0;
    let result_1 = filter_sam(&in1, &out1, &alignments, low, high, &correct_orientation, 1);
    match result_1 {
        Ok(count) => { after_count += count },
        Err(_) => quit_with_error(&format!("unable to write alignments to {:?}", out1)),
    }
    let result_2 = filter_sam(&in2, &out2, &alignments, low, high, &correct_orientation, 2);
    match result_2 {
        Ok(count) => { after_count += count },
        Err(_) => quit_with_error(&format!("unable to write alignments to {:?}", out2)),
    }
    after_count
}


fn filter_sam(in_filename: &PathBuf, out_filename: &PathBuf,
              alignments: &HashMap<String, Vec<Alignment>>, low: u32, high: u32,
              correct_orientation: &String, read_num: usize) -> io::Result<usize> {
    eprintln!("Filtering {}:", in_filename.display());
    let mut pass_count = 0;
    let mut fail_count = 0;

    let in_file = File::open(in_filename)?;
    let reader = io::BufReader::new(in_file);
    let out_file = File::create(out_filename)?;
    let mut writer = BufWriter::new(out_file);
    static NO_ALIGNMENTS: Vec<Alignment> = Vec::new();

    for line in reader.lines() {
        let sam_line = line?;
        if sam_line.starts_with('@') {
            writeln!(writer, "{}", sam_line)?;
            continue;
        }

        let a = Alignment::new_quick(&sam_line).unwrap();
        if !a.is_aligned() {
            writeln!(writer, "{}", sam_line)?;
            continue;
        }

        let (this_name, pair_name) = if read_num == 1 {
            (format!("{}_1", a.read_name), format!("{}_2", a.read_name))
        } else {
            (format!("{}_2", a.read_name), format!("{}_1", a.read_name))
        };

        let this_alignments = &alignments[&this_name];
        let pair_alignments = match alignments.get(&pair_name) {
            Some(alignments) => alignments,
            None => &NO_ALIGNMENTS,
        };

        if alignment_pass_qc(&a, this_alignments, pair_alignments, low, high, correct_orientation) {
            writeln!(writer, "{}", sam_line)?;
            pass_count += 1;
        } else {
            let mut parts: Vec<&str> = sam_line.split('\t').collect();
            parts.push("ZP:Z:fail");
            writeln!(writer, "{}", parts.join("\t"))?;
            fail_count += 1;
        }
    }

    eprintln!("  {} pass", pass_count.to_formatted_string(&Locale::en));
    eprintln!("  {} fail", fail_count.to_formatted_string(&Locale::en));
    eprintln!();
    Ok(pass_count)
}


fn alignment_pass_qc(a: &Alignment, this_alignments: &[Alignment], pair_alignments: &[Alignment],
                     low: u32, high: u32, correct_orientation: &str) -> bool {
    // Rules for whether an alignment passes or fails filtering:
    // * If there are no pair alignments, it passes. I.e. if we can't use read pairs to assess the
    //   alignment, we keep it.
    // * If there is exactly one alignment for this read, it passes. I.e. we're not going to throw
    //   out the only alignment for a read.
    // * If there are multiple alignments for this read and at least one pair alignment, then the
    //   alignment passes if it makes a good pair (same reference seq, good insert size and correct
    //   orientation) with any of the pair alignments.
    if pair_alignments.is_empty() {
        return true;
    }
    if this_alignments.len() == 1 {
        return true;
    }
    for pair_alignment in pair_alignments {
        let same_ref = a.ref_name == pair_alignment.ref_name;
        let insert = get_insert_size(a, pair_alignment);
        let orientation = get_orientation(a, pair_alignment);
        if same_ref && low <= insert && insert <= high && orientation == correct_orientation {
            return true;
        }
    }
    false
}


#[cfg(test)]
mod tests {
    use super::*;

    fn run_get_orientation_test(pos_1: i32, pos_2: i32,
                                strand_1: i32, strand_2: i32, result: &str) {
        let str_1 = format!("r_1\t{}\tx\t{}\t60\t150M\t*\t0\t0\tACTG\tKKKK\tNM:i:0",
                            strand_1, pos_1);
        let str_2 = format!("r_2\t{}\tx\t{}\t60\t150M\t*\t0\t0\tACTG\tKKKK\tNM:i:0",
                            strand_2, pos_2);
        let a_1 = Alignment::new_quick(&str_1).unwrap();
        let a_2 = Alignment::new_quick(&str_2).unwrap();
        assert_eq!(get_orientation(&a_1, &a_2), result);
    }

    #[test]
    fn test_get_orientation() {
        // 1------>
        //            <------2
        run_get_orientation_test(100000, 200000, 0, 16, "fr");

        // 2------>
        //            <------1
        run_get_orientation_test(200000, 100000, 16, 0, "fr");

        //            1------>
        // <------2
        run_get_orientation_test(200000, 100000, 0, 16, "rf");

        // <------1
        //            2------>
        run_get_orientation_test(100000, 200000, 16, 0, "rf");

        // 1------>   2------>
        run_get_orientation_test(100000, 200000, 0, 0, "ff");

        // <------2   <------1
        run_get_orientation_test(200000, 100000, 16, 16, "ff");

        // 2------>   1------>
        run_get_orientation_test(200000, 100000, 0, 0, "rr");

        // <------1   <------2
        run_get_orientation_test(100000, 200000, 16, 16, "rr");
    }

    #[test]
    fn test_auto_determine_orientation() {
        let insert_sizes: HashMap<String, Vec<u32>> = [
            ("fr", vec![100, 100, 100]), ("rf", vec![200]), ("ff", vec![300]), ("rr", vec![400])
        ].iter().map(|&(k, ref v)| (k.to_string(), v.clone())).collect();
        assert_eq!(auto_determine_orientation(&insert_sizes), "fr");

        let insert_sizes: HashMap<String, Vec<u32>> = [
            ("fr", vec![100]), ("rf", vec![200, 200, 200]), ("ff", vec![300]), ("rr", vec![400])
        ].iter().map(|&(k, ref v)| (k.to_string(), v.clone())).collect();
        assert_eq!(auto_determine_orientation(&insert_sizes), "rf");

        let insert_sizes: HashMap<String, Vec<u32>> = [
            ("fr", vec![100]), ("rf", vec![200]), ("ff", vec![300, 300, 300]), ("rr", vec![400])
        ].iter().map(|&(k, ref v)| (k.to_string(), v.clone())).collect();
        assert_eq!(auto_determine_orientation(&insert_sizes), "ff");

        let insert_sizes: HashMap<String, Vec<u32>> = [
            ("fr", vec![100]), ("rf", vec![200]), ("ff", vec![300]), ("rr", vec![400, 400, 400])
        ].iter().map(|&(k, ref v)| (k.to_string(), v.clone())).collect();
        assert_eq!(auto_determine_orientation(&insert_sizes), "rr");
    }

    #[test]
    fn test_get_percentile() {
        let nums: Vec<u32> = vec![15, 20, 35, 40, 50];
        assert_eq!(get_percentile(&nums, 00.1), 15);
        assert_eq!(get_percentile(&nums, 19.9), 15);
        assert_eq!(get_percentile(&nums, 20.1), 20);
        assert_eq!(get_percentile(&nums, 39.9), 20);
        assert_eq!(get_percentile(&nums, 40.1), 35);
        assert_eq!(get_percentile(&nums, 59.9), 35);
        assert_eq!(get_percentile(&nums, 60.1), 40);
        assert_eq!(get_percentile(&nums, 79.9), 40);
        assert_eq!(get_percentile(&nums, 80.1), 50);
        assert_eq!(get_percentile(&nums, 99.9), 50);
    }

    #[test]
    fn test_get_percentile_name() {
        assert_eq!(get_percentile_name(1.0), "1st percentile");
        assert_eq!(get_percentile_name(2.0), "2nd percentile");
        assert_eq!(get_percentile_name(3.0), "3rd percentile");
        assert_eq!(get_percentile_name(4.0), "4th percentile");
        assert_eq!(get_percentile_name(5.0), "5th percentile");
        assert_eq!(get_percentile_name(6.0), "6th percentile");
        assert_eq!(get_percentile_name(7.0), "7th percentile");
        assert_eq!(get_percentile_name(8.0), "8th percentile");
        assert_eq!(get_percentile_name(9.0), "9th percentile");
        assert_eq!(get_percentile_name(10.0), "10th percentile");
        assert_eq!(get_percentile_name(0.1), "0.1st percentile");
        assert_eq!(get_percentile_name(99.9), "99.9th percentile");
    }
}
