use std::iter;

use ansi_term::Colour::Yellow;
use chrono::NaiveTime;
use indicatif::{ProgressBar, ProgressStyle};
use num_format::{Locale, ToFormattedString};

const DIVIDER_LEN: usize = 57;

pub fn fmt_num(num: &usize) -> String {
    num.to_formatted_string(&Locale::en)
}

pub fn parse_duration(duration: u64) -> String {
    let sec = (duration % 60) as u32;
    let min = ((duration / 60) % 60) as u32;
    let hours = ((duration / 60) / 60) as u32;
    let time = NaiveTime::from_hms(hours, min, sec);
    time.format("%H:%M:%S").to_string()
}

#[cfg(not(tarpaulin_include))]
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

#[cfg(not(tarpaulin_include))]
pub fn print_welcome_text(version: &str) {
    log::info!("{}", Yellow.paint(get_rep_str('=')));
    let text = format!("SEGUL v{}", version);
    log::info!("{}", Yellow.paint(text));
    log::info!("{}", Yellow.paint("An alignment tool for phylogenomics"));
    log::info!("{}", Yellow.paint(get_rep_str('-')));
}

pub fn print_divider() {
    let divider = get_rep_str('-');
    log::info!("{}", Yellow.paint(divider));
}

fn get_rep_str(sym: char) -> String {
    iter::repeat(sym).take(DIVIDER_LEN).collect()
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

    #[test]
    fn test_num_fmt() {
        let num = 1300;
        let fmt_num = fmt_num(&num);
        assert_eq!("1,300", fmt_num);
    }
}
