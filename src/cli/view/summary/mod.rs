//! Interactive viewer of summary reports in the terminal.
mod app;
mod data;
mod ui;

use std::io::{self, stdout};
use std::path::Path;

use crossterm::event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyEventKind};
use crossterm::execute;
use segul_core::color::Theme;

use crate::cli::args::view::ViewSummaryArgs;

use self::app::App;
use self::data::SummaryData;

pub(in crate::cli) struct SummaryViewParser<'a> {
    args: &'a ViewSummaryArgs,
}

impl<'a> SummaryViewParser<'a> {
    pub(in crate::cli) fn new(args: &'a ViewSummaryArgs) -> Self {
        Self { args }
    }

    pub(in crate::cli) fn view(&self) {
        let path: &Path = &self.args.input;
        let data = match SummaryData::from_path(path) {
            Ok(data) => data,
            Err(e) => {
                eprintln!("{e:#}");
                std::process::exit(1);
            }
        };
        let file_name = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        let theme = self.args.theme.parse::<Theme>().unwrap_or(Theme::Dark);
        ui::set_accent(theme);
        let mut app = App::new(data, file_name);
        run(&mut app).expect("Failed running the summary viewer");
    }
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
        }
        Ok(())
    })();
    execute!(stdout(), DisableMouseCapture)?;
    ratatui::restore();
    result
}
