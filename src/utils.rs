use std::io::{self, Result, Write};

use num_format::{Locale, ToFormattedString};

pub fn format_thousand_sep(num: &usize) -> String {
    num.to_formatted_string(&Locale::en)
}

pub fn print_divider(text: &str, len: usize) {
    let sym = '=';
    let mut header = PrettyDivider::new(text, sym, len);
    header.print_header().unwrap();
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
