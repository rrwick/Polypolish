// Copyright 2021 Ryan Wick (rrwick@gmail.com)
// https://github.com/rrwick/Polypolish

//This file is part of Polypolish. Polypolish is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version. Polypolish
// is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details. You should have received a copy of the GNU General Public
// License along with Polypolish. If not, see <http://www.gnu.org/licenses/>.

use chrono::prelude::*;
use colored::Colorize;


pub fn section_header(text: &str) {
    let now = Local::now().format("%Y-%m-%d %H:%M:%S").to_string();
    let date = format!("({})", now);
    eprintln!();
    eprintln!("{} {}", text.bold().yellow().underline(), date.dimmed());
}


pub fn explanation(text: &str) {
    let mut term_width = 80;
    if let Some((w, _)) = term_size::dimensions_stderr() {
        term_width = w;
    }
    let indented_text = format!("    {}", text);
    eprintln!("{}", textwrap::fill(&indented_text, term_width).dimmed());
    eprintln!();
}


pub fn ascii_art() -> String {
    let s1 = r#"  _____        _                       _  _       _     "#;
    let s2 = r#" |  __ \      | |                     | |(_)     | |    "#;
    let s3 = r#" | |__) |___  | | _   _  _ __    ___  | | _  ___ | |__  "#;
    let s4 = r#" |  ___// _ \ | || | | || '_ \  / _ \ | || |/ __|| '_ \ "#;
    let s5 = r#" | |   | (_) || || |_| || |_) || (_) || || |\__ \| | | |"#;
    let s6 = r#" |_|    \___/ |_| \__, || .__/  \___/ |_||_||___/|_| |_|"#;
    let s7 = r#"                   __/ || |                             "#;
    let s8 = r#"                  |___/ |_|                             "#;
    format!("{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}", s1, s2, s3, s4, s5, s6, s7, s8)
}
