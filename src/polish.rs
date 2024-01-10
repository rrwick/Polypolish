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
use std::collections::HashMap;
use std::time::Instant;
use std::fs::File;
use std::io::prelude::*;
use clap::crate_version;
use num_format::{Locale, ToFormattedString};

use crate::alignment;
use crate::log;
use crate::misc;
use crate::pileup;


pub fn polish(debug: Option<PathBuf>, fraction_invalid: f64, fraction_valid: f64, max_errors: u32,
              min_depth: u32, assembly: PathBuf, sam: Vec<PathBuf>) {
    let start_time = Instant::now();
    check_option_values(&fraction_invalid, &fraction_valid);
    check_inputs_exist(&assembly, &sam);
    starting_message(&debug, &fraction_invalid, &fraction_valid, &max_errors, &min_depth,
                     &assembly, &sam);
    let (seq_names, mut pileups) = load_assembly(&assembly);
    load_alignments(&max_errors, &sam, &mut pileups);
    let new_lengths = polish_sequences(&debug, &fraction_invalid, &fraction_valid, &min_depth,
                                       &seq_names, &mut pileups);
    finished_message(&debug, new_lengths, start_time);
}


fn starting_message(debug: &Option<PathBuf>, fraction_invalid: &f64, fraction_valid: &f64,
                    max_errors: &u32, min_depth: &u32, assembly: &PathBuf, sam: &Vec<PathBuf>) {
    log::section_header("Starting Polypolish polish");
    log::explanation("Polypolish is a tool for polishing genome assemblies with short reads. \
                      Unlike other tools in this category, Polypolish uses SAM files where each \
                      read has been aligned to all possible locations (not just a single best \
                      location). This allows it to repair errors in repeat regions that other \
                      alignment-based polishers cannot fix.");
    eprintln!("Polypolish version: {}", crate_version!());
    eprintln!();
    eprintln!("Input assembly:");
    eprintln!("  {}", assembly.display());
    eprintln!();
    eprintln!("Input short-read alignments:");
    for s in sam {
        eprintln!("  {}", s.display());
    }
    eprintln!();
    eprintln!("Settings:");
    eprintln!("  --fraction_invalid {}", fraction_invalid);
    eprintln!("  --fraction_valid {}", fraction_valid);
    eprintln!("  --max_errors {}", max_errors);
    eprintln!("  --min_depth {}", min_depth);
    match debug {
        Some(filename) => eprintln!("  --debug {}", filename.display()),
        None           => eprintln!("  not logging debugging information"),
    }
    eprintln!();
}


fn finished_message(debug: &Option<PathBuf>, new_lengths: Vec<(String, usize)>,
                    start_time: Instant) {
    log::section_header("Finished!");
    eprintln!("Polished sequence (to stdout):");
    for (new_name, new_length) in new_lengths {
        eprintln!("  {}_polypolish ({} bp)", new_name, new_length.to_formatted_string(&Locale::en));
    }
    eprintln!();
    match debug {
        Some(filename) => eprintln!("Per-base debugging info written to {}", filename.display()),
        None           => {},
    }
    eprintln!("Time to run: {}", misc::format_duration(start_time.elapsed()));
    eprintln!();
}


fn load_assembly(assembly_filename: &PathBuf) -> (Vec<String>, HashMap<String, pileup::Pileup>) {
    log::section_header("Loading assembly");
    let fasta = misc::load_fasta(assembly_filename);
    let mut seq_names = Vec::new();
    let mut pileups = HashMap::new();
    for (name, sequence) in &fasta {
        eprintln!("{} ({} bp)", name, sequence.len().to_formatted_string(&Locale::en));
        seq_names.push(name.clone());
        pileups.insert(name.clone(), pileup::Pileup::new(sequence));
    }
    eprintln!();
    (seq_names, pileups)
}


fn load_alignments(max_errors: &u32, sam: &Vec<PathBuf>,
                   pileups: &mut HashMap<String, pileup::Pileup>) {
    log::section_header("Loading alignments");
    let mut alignment_total: usize = 0;
    let mut used_total: usize = 0;
    for s in sam {
        let (alignment_count, used_count, read_count) = alignment::process_sam(&s, pileups,
                                                                               *max_errors);
        eprintln!("{}: {} alignments from {} reads", s.display(),
                  alignment_count.to_formatted_string(&Locale::en),
                  read_count.to_formatted_string(&Locale::en));
        alignment_total += alignment_count;
        used_total += used_count;
    }
    let discarded_count = alignment_total - used_total;
    eprintln!();
    eprintln!("Filtering for high-quality end-to-end alignments:");
    eprintln!("  {} alignments kept", used_total.to_formatted_string(&Locale::en));
    eprintln!("  {} alignments discarded", discarded_count.to_formatted_string(&Locale::en));
    eprintln!();
}


fn polish_sequences(debug: &Option<PathBuf>, fraction_invalid: &f64, fraction_valid: &f64,
                    min_depth: &u32, seq_names: &Vec<String>,
                    pileups: &HashMap<String, pileup::Pileup>) -> Vec<(String, usize)>{
    log::section_header("Polishing assembly sequences");
    log::explanation("For each position in the assembly, Polypolish determines the read \
                     depth at that position and collects all aligned bases. It then polishes the \
                     assembly by looking for positions where the pileup unambiguously supports a \
                     different sequence than the assembly.");
    let mut debug_file = create_debug_file(&debug);
    let mut new_lengths = Vec::new();
    for name in seq_names {
        let pileup = pileups.get(name).unwrap();
        let new_length = polish_one_sequence(&debug, &fraction_invalid, &fraction_valid, &min_depth,
                                             name, pileup, &mut debug_file);
        new_lengths.push((name.clone(), new_length));
    }
    new_lengths
}


fn polish_one_sequence(debug: &Option<PathBuf>, fraction_invalid: &f64, fraction_valid: &f64,
                       min_depth: &u32, name: &str, pileup: &pileup::Pileup,
                       debug_file: &mut Option<File>) -> usize {
    let seq_len = pileup.bases.len();
    eprintln!("Polishing {} ({} bp):", name, seq_len.to_formatted_string(&Locale::en));

    let mut polished_seq: String = String::with_capacity(seq_len);
    let mut total_depth = 0.0;
    let mut zero_depth_count: usize = 0;
    let mut changed_count: usize = 0;
    let mut pos: usize = 0;
    let build_debug_str = match debug_file {Some(_) => true, None => false};

    for b in &pileup.bases {
        let (seq, status, debug_line) = b.get_polished_seq(*min_depth, *fraction_valid,
                                                           *fraction_invalid, build_debug_str);
        match status {
            pileup::BaseStatus::Changed => {changed_count += 1}
            _                           => {}
        }
        total_depth += b.depth;
        if b.depth == 0.0 {
            zero_depth_count += 1;
        }
        match debug_file {
            Some(file) => write_debug_line(file, name, pos, &debug_line, &debug),
            None       => {},
        }
        polished_seq.push_str(&seq);
        pos += 1;
    }
    polished_seq = polished_seq.replace("-", "");
    println!(">{}_polypolish", name);
    println!("{}", polished_seq);

    print_polishing_info(seq_len, total_depth, zero_depth_count, changed_count);

    polished_seq.len()
}


fn print_polishing_info(seq_len: usize, total_depth: f64, zero_depth_count: usize,
                        changed_count: usize) {
    let seq_len_f64 = seq_len as f64;
    let mean_depth = total_depth / seq_len_f64;
    eprintln!("  mean read depth: {:.1}x", mean_depth);

    let have = if zero_depth_count == 1 {"has"} else {"have"};
    let covered = seq_len - zero_depth_count;
    let coverage = 100.0 * (covered as f64) / seq_len_f64;
    eprintln!("  {} bp {} a depth of zero ({:.4}% coverage)",
              zero_depth_count.to_formatted_string(&Locale::en), have, coverage);

    let changed_percent = 100.0 * (changed_count as f64) / seq_len_f64;
    let estimated_accuracy = 100.0 - changed_percent;
    let positions = if changed_count == 1 {"position"} else {"positions"};
    eprintln!("  {} {} changed ({:.4}% of total positions)",
              changed_count.to_formatted_string(&Locale::en), positions, changed_percent);
    eprintln!("  estimated pre-polishing sequence accuracy: {:.4}%", estimated_accuracy);
    eprintln!();
}


fn create_debug_file(debug: &Option<PathBuf>) -> Option<File> {
    match debug {
        Some(_) => {},
        None    => {return None;},
    }
    let filename = debug.as_ref().unwrap();
    let create_result = File::create(filename);
    match create_result {
        Ok(_)  => (),
        Err(_) => misc::quit_with_error(&format!("unable to create {:?}", filename)),
    }
    let mut file = create_result.unwrap();
    write_debug_header(&mut file, filename);
    Some(file)
}


fn write_debug_header(file: &mut File, filename: &PathBuf) {
    let header = "name\tpos\tbase\tdepth\tinvalid\tvalid\tpileup\tstatus\tnew_base\n";
    let result = file.write_all(header.as_bytes());
    match result {
        Ok(_)  => (),
        Err(_) => misc::quit_with_error(&format!("unable to write to file {:?}", filename)),
    }
}


fn write_debug_line(file: &mut File, name: &str, pos: usize, debug_line: &str,
                    debug: &Option<PathBuf>) {
    let debug_line: String = format!("{}\t{}\t{}\n", name, pos, debug_line);
    let result = file.write_all(debug_line.as_bytes());
    match result {
        Ok(_)  => (),
        Err(_) => misc::quit_with_error(&format!("unable to write to file {:?}",
                                                 debug.as_ref().unwrap())),
    }
}


fn check_inputs_exist(assembly: &PathBuf, sam: &Vec<PathBuf>) {
    misc::check_if_file_exists(&assembly);
    for s in sam {
        misc::check_if_file_exists(&s);
    }
}


fn check_option_values(fraction_invalid: &f64, fraction_valid: &f64) {
    if *fraction_valid <= 0.0 || *fraction_valid >= 1.0 {
        misc::quit_with_error("--fraction_valid must be between 0 and 1 (exclusive)")
    }
    if *fraction_invalid <= 0.0 || *fraction_invalid >= 1.0 {
        misc::quit_with_error("--fraction_invalid must be between 0 and 1 (exclusive)")
    }
    if *fraction_invalid >= *fraction_valid {
        misc::quit_with_error("--fraction_invalid must be less than --fraction_valid")
    }
}
