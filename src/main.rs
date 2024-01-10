// Copyright 2021 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Polypolish

// This file is part of Polypolish. Polypolish is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Polypolish
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Polypolish. If not, see <http://www.gnu.org/licenses/>.

mod alignment;
mod filter;
mod log;
mod misc;
mod pileup;
mod polish;

use std::path::PathBuf;
use clap::{Parser, Subcommand, crate_version};


#[derive(Parser)]
#[clap(name = "Polypolish",
       version = concat!("v", crate_version!()),
       about = "short-read polishing of long-read assemblies\ngithub.com/rrwick/Polypolish",
       before_help = concat!(r#"  _____        _                       _  _       _     "#, "\n",
                             r#" |  __ \      | |                     | |(_)     | |    "#, "\n",
                             r#" | |__) |___  | | _   _  _ __    ___  | | _  ___ | |__  "#, "\n",
                             r#" |  ___// _ \ | || | | || '_ \  / _ \ | || |/ __|| '_ \ "#, "\n",
                             r#" | |   | (_) || || |_| || |_) || (_) || || |\__ \| | | |"#, "\n",
                             r#" |_|    \___/ |_| \__, || .__/  \___/ |_||_||___/|_| |_|"#, "\n",
                             r#"                   __/ || |                             "#, "\n",
                             r#"                  |___/ |_|                             "#))]
#[command(author, version, about, long_about = None, disable_help_subcommand = true,
          propagate_version = true)]
#[clap(subcommand_required = true)]
#[clap(arg_required_else_help = true)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// filter paired-end alignments based on insert size
    Filter {
        /// Input SAM file - first read in pairs
        #[clap(long = "in1")]
        in1: PathBuf,

        /// Input SAM file - first second in pairs
        #[clap(long = "in2")]
        in2: PathBuf,
    
        /// Output SAM file - first read in pairs
        #[clap(long = "out1")]
        out1: PathBuf,

        /// Output SAM file - first second in pairs
        #[clap(long = "out2")]
        out2: PathBuf,

        /// Expected pair orientation
        #[clap(long = "orientation", default_value = "auto")]
        orientation: String,

        /// Low percentile threshold
        #[clap(long = "low", default_value = "0.1")]
        low: f64,

        /// High percentile threshold
        #[clap(long = "high", default_value = "0.1")]
        high: f64,
    },

    /// polish a long-read assembly using short-read alignments
    Polish {
        /// Optional file to store per-base information for debugging purposes
        #[clap(long = "debug")]
        debug: Option<PathBuf>,

        /// A base must make up less than this fraction of the read depth to be considered invalid
        #[clap(short = 'i', long = "fraction_invalid", default_value = "0.2")]
        fraction_invalid: f64,

        /// A base must make up at least this fraction of the read depth to be considered valid
        #[clap(short = 'v', long = "fraction_valid", default_value = "0.5")]
        fraction_valid: f64,

        /// Ignore alignments with more than this many mismatches and indels
        #[clap(short = 'm', long = "max_errors", default_value = "10")]
        max_errors: u32,

        /// A base must occur at least this many times in the pileup to be considered valid
        #[clap(short = 'd', long = "min_depth", default_value = "5")]
        min_depth: u32,

        /// Assembly to polish (one file in FASTA format)
        assembly: PathBuf,

        /// Short read alignments (one or more files in SAM format)
        sam: Vec<PathBuf>,
    },
}


fn main() {
    let cli = Cli::parse();

    match cli.command {
        Some(Commands::Filter { in1, in2, out1, out2, orientation, low, high }) => {
            filter::filter(in1, in2, out1, out2, orientation, low, high);
        },
        Some(Commands::Polish { debug, fraction_invalid, fraction_valid, max_errors, min_depth,
                                assembly, sam}) => {
            polish::polish(debug, fraction_invalid, fraction_valid, max_errors, min_depth,
                           assembly, sam);
        },
        None => {}
    }
}
