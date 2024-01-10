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
use clap::{Parser, Subcommand, crate_version};
use num_format::{Locale, ToFormattedString};

use crate::alignment;
use crate::log;
use crate::misc;
use crate::pileup;


pub fn filter(in1: PathBuf, in2: PathBuf, out1: PathBuf, out2: PathBuf,
              orientation: String, low: f64, high: f64) {

}
