// Copyright 2021 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Polypolish

//This file is part of Polypolish. Polypolish is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Polypolish
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Polypolish. If not, see <http://www.gnu.org/licenses/>.

mod log;
mod misc;
mod pileup;
mod alignment;

use std::path::PathBuf;
use std::collections::HashMap;
use std::time::Instant;
use clap::{AppSettings, Clap};
use num_format::{Locale, ToFormattedString};


#[derive(Clap)]
#[clap(name = log::ascii_art(), about = "short-read polishing of long-read assemblies\n\
                                         github.com/rrwick/Polypolish")]
#[clap(setting = AppSettings::ColoredHelp)]
struct Opts {
    /// Optional file to store per-base information for debugging purposes
    #[clap(long = "debug")]
    debug: Option<PathBuf>,

    /// Ignore alignments with more than this many mismatches and indels
    #[clap(short = 'm', long = "max_errors", default_value = "10")]
    max_errors: i32,

    /// A base must occur at least this many times in the pileup to be considered valid
    #[clap(short = 'd', long = "min_depth", default_value = "5")]
    min_depth: usize,

    /// A base must make up at least this fraction of the pileup to be considered valid
    #[clap(short = 'f', long = "min_fraction", default_value = "0.5")]
    min_fraction: f64,

    /// Assembly to polish (FASTA format)
    #[clap(parse(from_os_str), required = true)]
    assembly: PathBuf,

    /// Short read alignments (SAM format, one or more files)
    #[clap(parse(from_os_str), required = true)]
    sam: Vec<PathBuf>,
}


fn main() {
    let opts: Opts = Opts::parse();
    let start_time = Instant::now();
    check_inputs_exist(&opts);
    starting_message(&opts);
    let (seq_names, mut pileups) = load_assembly(&opts.assembly);
    load_alignments(&opts, &mut pileups);


    let mut new_lengths = Vec::new();
    for name in &seq_names {
        let pileup = pileups.get(name).unwrap();
        let mut polished_seq: String = String::with_capacity(pileup.bases.len());
        for b in &pileup.bases {
            let out_seq = b.get_output_seq(opts.min_depth, opts.min_fraction);
            polished_seq.push_str(&out_seq);
        }
        polished_seq = polished_seq.replace("-", "");
        println!(">{}_polypolish", name);
        println!("{}", polished_seq);
        new_lengths.push((name, polished_seq.len()));
    }
    finished_message(new_lengths, start_time);
}


fn starting_message(opts: &Opts) {
    log::section_header("Starting Polypolish");
    log::explanation("Polypolish is a tool for polishing genome assemblies with short reads. \
                      Unlike other tools in this category, Polypolish uses SAM files where each \
                      read has been aligned to all possible locations (not just a single best \
                      location). This allows it to repair errors in repeat regions that other \
                      alignment-based polishers cannot fix.");

    const VERSION: &'static str = env!("CARGO_PKG_VERSION");
    eprintln!("Polypolish version: {}", VERSION);
    eprintln!();
    eprintln!("Input assembly:");
    eprintln!("  {}", opts.assembly.display());
    eprintln!();
    eprintln!("Input short-read alignments:");
    for s in &opts.sam {
        eprintln!("  {}", s.display());
    }
    eprintln!();
    eprintln!("Settings:");
    eprintln!("  --max_errors {}", opts.max_errors);
    eprintln!("  --min_depth {}", opts.min_depth);
    eprintln!("  --min_fraction {}", opts.min_fraction);
    match &opts.debug {
        Some(v) => {eprintln!("  --debug {}", v.display());}
        None => {eprintln!("  not logging debugging information");},
    }
    eprintln!();
}


fn finished_message(new_lengths: Vec<(&String, usize)>, start_time: Instant) {
    log::section_header("Finished!");
    eprintln!("Polished sequence (to stdout):");
    for (new_name, new_length) in new_lengths {
        eprintln!("  {}_polypolish ({} bp)", new_name, new_length.to_formatted_string(&Locale::en));
    }
    eprintln!();
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


fn load_alignments(opts: &Opts, pileups: &mut HashMap<String, pileup::Pileup>) {
    log::section_header("Loading alignments");
    let mut alignment_total: usize = 0;
    let mut used_total: usize = 0;
    for s in &opts.sam {
        let (alignment_count, used_count, read_count) = alignment::process_sam(&s, pileups, opts.max_errors);
        eprintln!("{}: {} alignments from {} reads",
                  s.as_path().display().to_string(),
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


fn check_inputs_exist(opts: &Opts) {
    misc::check_if_file_exists(&opts.assembly);
    for s in &opts.sam {
        misc::check_if_file_exists(&s);
    }
}
