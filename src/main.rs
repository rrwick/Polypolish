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
use structopt::StructOpt;


#[derive(Debug, StructOpt)]
#[structopt(name = "Polypolish", about = "short-read polishing of long-read assemblies\n\
                                          github.com/rrwick/Polypolish")]
struct Opt {
    #[structopt(long = "debug",
                help = "Optional file to store per-base information for debugging purposes")]
    debug: Option<PathBuf>,

    #[structopt(short = "m", long = "max_errors", default_value = "10",
                help = "Ignore alignments with more than this many mismatches and indels")]
    max_errors: i32,

    #[structopt(short = "d", long = "min_depth", default_value = "5",
                help = "A base must occur at least this many times in the pileup to be considered \
                        valid")]
    min_depth: i32,

    #[structopt(short = "f", long = "min_fraction", default_value = "0.5",
                help = "A base must make up at least this fraction of the pileup to be considered \
                        valid")]
    min_fraction: f64,

    #[structopt(parse(from_os_str),
                help = "Assembly to polish (FASTA format)")]
    assembly: PathBuf,

    #[structopt(parse(from_os_str),
                help = "Short read alignments (SAM format, one or more files)")]
    sam: Vec<PathBuf>,
}


fn main() {
    let opt = Opt::from_args();
    println!("{:?}", opt);
}
