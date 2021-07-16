// Copyright 2021 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Polypolish

//This file is part of Polypolish. Polypolish is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Polypolish
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Polypolish. If not, see <http://www.gnu.org/licenses/>.

use std::path::Path;
use std::path::PathBuf;


pub fn check_if_file_exists(filename: &PathBuf) {
    if !Path::new(filename).exists() {
        let error_message = format!("{:?} does not exist", filename);
        quit_with_error(&error_message);
    }
}


pub fn quit_with_error(text: &String) {
    eprintln!();
    eprintln!("Error: {}", text);
    std::process::exit(1);
}