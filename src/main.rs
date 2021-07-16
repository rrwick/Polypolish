// Copyright 2021 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Polypolish

//This file is part of Polypolish. Polypolish is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Polypolish
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Polypolish. If not, see <http://www.gnu.org/licenses/>.

use std::path::PathBuf;
use clap::{AppSettings, Clap};


#[derive(Clap)]
#[clap(name = "Polypolish", about = "short-read polishing of long-read assemblies\n\
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
    #[clap(parse(from_os_str))]
    assembly: PathBuf,

    /// Short read alignments (SAM format, one or more files)
    #[clap(parse(from_os_str))]
    sam: Vec<PathBuf>,
}


fn main() {
    let opts: Opts = Opts::parse();
    starting_message(&opts);
}



fn starting_message(opts: &Opts) {
    println!("Input assembly:");
    println!("  {}", opts.assembly.display());
    println!();
    println!("Input short-read alignments:");
    for s in &opts.sam {
        println!("  {}", s.display());
    }
    println!();
    println!("Settings:");
    println!("  --max_errors {}", opts.max_errors);
    println!("  --min_depth {}", opts.min_depth);
    println!("  --min_fraction {}", opts.min_fraction);
    match &opts.debug {
        Some(v) => {println!("  --debug {}", v.display());}
        None => {println!("  not logging debugging information");},
    }
    println!();
}
