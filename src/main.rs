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

use std::path::PathBuf;
use clap::{AppSettings, Clap};
use bio::io::fasta;
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
    min_depth: i32,

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
    check_inputs_exist(&opts);
    starting_message(&opts);
    load_assembly(&opts.assembly);
}


fn check_inputs_exist(opts: &Opts) {
    misc::check_if_file_exists(&opts.assembly);
    for s in &opts.sam {
        misc::check_if_file_exists(&s);
    }
}


fn starting_message(opts: &Opts) {
    log::section_header("Starting Polypolish");
    log::explanation("Polypolish is a tool for polishing genome assemblies with short reads. \
                      Unlike other tools in this category, Polypolish uses SAM files where each \
                      read has been aligned to all possible locations (not just a single best \
                      location). This allows it to repair errors in repeat regions that other \
                      alignment-based polishers cannot fix.");

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


fn load_assembly(assembly_filename: &PathBuf) {
    log::section_header("Loading assembly");

    let result = fasta::Reader::from_file(assembly_filename);
    match result {
        Ok(_) => (),
        Err(ref e) => misc::quit_with_error(&e.to_string())
    }
    let mut records = result.unwrap().records();

    while let Some(Ok(record)) = records.next() {
        eprintln!("{} ({} bp)", record.id(), record.seq().len().to_formatted_string(&Locale::en));
        // TODO: save the sequence in a vector or something
    }
    // TODO: return the sequences
}
