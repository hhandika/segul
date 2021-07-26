use std::io::{self, BufWriter, Result, Write};
use std::iter;

use chrono::NaiveTime;
use indicatif::{ProgressBar, ProgressStyle};
use num_format::{Locale, ToFormattedString};

pub fn fmt_num(num: &usize) -> String {
    num.to_formatted_string(&Locale::en)
}

fn parse_duration(duration: u64) -> String {
    let sec = (duration % 60) as u32;
    let min = ((duration / 60) % 60) as u32;
    let hours = ((duration / 60) / 60) as u32;
    let time = NaiveTime::from_hms(hours, min, sec);
    time.format("%H:%M:%S").to_string()
}

pub fn print_formatted_duration(duration: u64) {
    let time = parse_duration(duration);
    log::info!("Execution time (HH:MM:SS): {}\n", time);
}

pub fn set_spinner() -> ProgressBar {
    let spin = ProgressBar::new_spinner();
    spin.enable_steady_tick(150);
    spin.set_style(
        ProgressStyle::default_spinner()
            .tick_chars("🌑🌒🌓🌔🌕🌖🌗🌘")
            .template("{spinner} {msg}"),
    );
    spin
}

pub fn print_title(text: &str) {
    let sym = '=';
    let len = 50;
    let mut header = PrettyDivider::new(text, sym, len);
    header.print_header().unwrap();
}

pub fn print_divider() {
    let divider: String = iter::repeat('-').take(52).collect();
    let io = io::stdout();
    let mut writer = BufWriter::new(io);
    writeln!(writer, "\n\x1b[0;33m{}\x1b[0m", divider).expect("Failed writing to stdout");
}

struct PrettyDivider {
    text: String,
    sym: char,
    len: usize,
    text_len: usize,
    sym_len: usize,
    color: String,
}

impl PrettyDivider {
    fn new(text: &str, sym: char, len: usize) -> Self {
        Self {
            text: String::from(text),
            sym,
            len,
            text_len: 0,
            sym_len: 0,
            color: String::from("\x1b[0;33m"),
        }
    }

    fn print_header(&mut self) -> Result<()> {
        self.get_len();
        let io = io::stdout();
        let mut handle = BufWriter::new(io);
        write!(handle, "{}", self.color)?;
        if self.text_len > self.len {
            writeln!(handle, "{}", self.text)?;
        } else {
            self.print_with_symbol(&mut handle)?;
        }
        write!(handle, "\x1b[0m")?;
        Ok(())
    }

    fn print_with_symbol<W: Write>(&mut self, handle: &mut W) -> Result<()> {
        self.print_symbols(handle);
        write!(handle, " {} ", self.text)?;
        self.print_symbols(handle);

        if self.text_len % 2 != 0 {
            write!(handle, "{}", self.sym)?;
        }

        writeln!(handle)?;
        Ok(())
    }

    fn get_len(&mut self) {
        self.text_len = self.text.len();

        if self.len > self.text_len {
            self.sym_len = (self.len - self.text_len) / 2;
        } else {
            self.sym_len = self.len;
        }
    }

    fn print_symbols<W: Write>(&self, io: &mut W) {
        let sym: String = iter::repeat(self.sym).take(self.sym_len).collect();
        write!(io, "{}", sym).unwrap();
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn time_parsing_test() {
        let duration = 65;
        let duration_2 = 3600;
        let time = parse_duration(duration);
        let hours = parse_duration(duration_2);

        assert_eq!("00:01:05", time);
        assert_eq!("01:00:00", hours);
    }
}
