//! Interactive alignment viewer in the terminal.
mod app;
mod palette;
mod ui;

use std::io::{self, stdout};
use std::path::PathBuf;

use crossterm::clipboard::CopyToClipboard;
use crossterm::event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyEventKind};
use crossterm::execute;

use crate::{
    cli::{AlignSeqInput, InputCli, args::view::ViewAlignArgs, collect_paths},
    core::align::view::AlignmentView,
};

use self::app::App;
use self::palette::{Palette, Theme};

pub(in crate::cli) struct AlignViewParser<'a> {
    args: &'a ViewAlignArgs,
    input_dir: Option<PathBuf>,
}

impl InputCli for AlignViewParser<'_> {}
impl AlignSeqInput for AlignViewParser<'_> {}

impl<'a> AlignViewParser<'a> {
    pub(in crate::cli) fn new(args: &'a ViewAlignArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn view(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.input_fmt);
        let datatype = self.parse_datatype(&self.args.datatype);
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        if files.is_empty() {
            eprintln!("Failed finding alignment files to view");
            std::process::exit(1);
        }
        let view = match AlignmentView::open(&files[0], &input_fmt, &datatype) {
            Ok(view) => view,
            Err(e) => {
                eprintln!("{e}");
                std::process::exit(1);
            }
        };
        let theme = match self.args.theme.as_str() {
            "light" => Theme::Light,
            _ => Theme::Dark,
        };
        let mut app = App::new(&files, &input_fmt, &datatype, view, self.args.chunk);
        run(&mut app, theme).expect("Failed running the alignment viewer");
    }
}

fn run(app: &mut App, theme: Theme) -> io::Result<()> {
    let is_truecolor = palette::is_truecolor();
    // `ratatui::init` restores the screen on panic; also release the mouse.
    let hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let _ = execute!(stdout(), DisableMouseCapture);
        hook(info);
    }));
    let mut terminal = ratatui::init();
    execute!(stdout(), EnableMouseCapture)?;
    let result = (|| -> io::Result<()> {
        let mut palette = Palette::new(theme, app.view.is_dna(), is_truecolor);
        let mut is_dna = app.view.is_dna();
        while !app.should_quit {
            if app.view.is_dna() != is_dna {
                is_dna = app.view.is_dna();
                palette = Palette::new(theme, is_dna, is_truecolor);
            }
            terminal.draw(|frame| ui::render(frame, app, &palette))?;
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
