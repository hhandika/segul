//! Utility functions for the CLI.
use std::{iter, time::Duration};

use chrono::NaiveTime;
use colored::Colorize;
use indicatif::{ProgressBar, ProgressStyle};
use num_format::{Locale, ToFormattedString};

const DIVIDER_LEN: usize = 57;

pub fn fmt_num(num: &usize) -> String {
    num.to_formatted_string(&Locale::en)
}

pub fn print_execution_time(duration: Duration) {
    if duration.as_secs() < 60 {
        log::info!("{:18}: {:?}", "Execution time", duration);
    } else {
        let time = parse_duration(duration.as_secs());
        log::info!("{:18}: {}", "Execution time (HH:MM:SS)", time);
    }
}

pub fn parse_duration(duration: u64) -> String {
    let sec = (duration % 60) as u32;
    let min = ((duration / 60) % 60) as u32;
    let hours = ((duration / 60) / 60) as u32;
    let time = NaiveTime::from_hms_opt(hours, min, sec);
    match time {
        Some(t) => t.format("%H:%M:%S").to_string(),
        None => "00:00:00".to_string(),
    }
}

#[cfg(not(tarpaulin_include))]
pub fn set_spinner() -> ProgressBar {
    let spin = ProgressBar::new_spinner();
    let duration: Duration = Duration::from_millis(150);
    spin.enable_steady_tick(duration);
    spin.set_style(
        ProgressStyle::default_spinner()
            .tick_chars("🌑🌒🌓🌔🌕🌖🌗🌘")
            .template("{spinner} {msg}")
            .expect("Failed getting progress bar."),
    );
    spin
}

#[cfg(not(tarpaulin_include))]
pub fn print_welcome_text(version: &str) {
    use clap::crate_description;

    log::info!("{}", get_rep_str('=').yellow());
    let text = format!("SEGUL v{}", version);
    log::info!("{}", text.yellow());
    log::info!("{}", crate_description!().yellow());
    log::info!("{}", get_rep_str('-').yellow());
}

pub fn print_divider() {
    let divider = get_rep_str('-');
    log::info!("{}", divider.yellow());
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
