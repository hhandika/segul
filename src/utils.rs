use std::io::{self, Result, Write};
use std::iter;

use indicatif::{ProgressBar, ProgressStyle};
use num_format::{Locale, ToFormattedString};

pub fn fmt_num(num: &usize) -> String {
    num.to_formatted_string(&Locale::en)
}

pub fn set_spinner() -> ProgressBar {
    let spin = ProgressBar::new_spinner();
    spin.enable_steady_tick(150);
    spin.set_style(ProgressStyle::default_spinner().template("{spinner:.simpleDots} {msg}"));
    spin
}

pub fn print_title(text: &str) {
    let sym = '=';
    let len = 50;
    let mut header = PrettyDivider::new(text, sym, len);
    header.print_header().unwrap();
}

#[allow(dead_code)]
pub fn print_divider<W: Write>(writer: &mut W) -> Result<()> {
    let divider: String = iter::repeat('-').take(50).collect();
    writeln!(writer, "{}", divider)?;

    Ok(())
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
        let mut handle = io::BufWriter::new(io);
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
        (0..=self.sym_len).for_each(|_| {
            write!(io, "{}", self.sym).unwrap();
        });
    }
}
