//! Interactive viewer of read reports in the terminal, matching the read
//! report of the SEGUI app: a summary of every sample with its quality
//! checks, then charts of the open sample by tab.
mod app;
mod metric;
mod spread;
mod ui;

use std::io::{self, stdout};

use crossterm::clipboard::CopyToClipboard;
use crossterm::event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyEventKind};
use crossterm::execute;
use segul_core::color::{Colormap, Theme};
use segul_core::stats::fastq::report::ReadReport;

use self::app::App;
use self::metric::samples_with_reads;

/// Shows [report], or says so when no sample has reads.
pub(super) fn view(report: ReadReport, file_name: String, theme: Theme, colormap: Colormap) {
    let qc_profile = report.qc_profile;
    let pairs = report.pairs.clone();
    let cohort = report.cohort.clone();
    let samples = samples_with_reads(report);
    if samples.is_empty() {
        eprintln!("{file_name} has no reads to show.");
        std::process::exit(1);
    }
    let mut app =
        App::new(samples, file_name, theme, colormap).with_report_meta(qc_profile, pairs, cohort);
    run(&mut app).expect("Failed running the read report viewer");
}

fn run(app: &mut App) -> io::Result<()> {
    // `ratatui::init` restores the screen on panic; also release the mouse.
    let hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let _ = execute!(stdout(), DisableMouseCapture);
        hook(info);
    }));
    let mut terminal = ratatui::init();
    execute!(stdout(), EnableMouseCapture)?;
    let result = (|| -> io::Result<()> {
        while !app.should_quit {
            terminal.draw(|frame| ui::render(frame, app))?;
            match event::read()? {
                Event::Key(key) if key.kind != KeyEventKind::Release => app.handle_key(key),
                Event::Mouse(mouse) => app.handle_mouse(mouse),
                _ => (),
            }
            if let Some(text) = app.clipboard.take() {
                execute!(stdout(), CopyToClipboard::to_clipboard_from(text))?;
            }
        }
        Ok(())
    })();
    execute!(stdout(), DisableMouseCapture)?;
    ratatui::restore();
    result
}
