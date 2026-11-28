//! Interactive tree viewer in the terminal.
mod app;
mod layers;
mod ui;

use std::io::{self, stdout};
use std::path::{Path, PathBuf};
use std::time::Duration;

use crossterm::clipboard::CopyToClipboard;
use crossterm::event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyEventKind};
use crossterm::execute;
use segul_core::color::Theme;
use segul_core::core::tree::view::{TreeView, read_tree_at};
use segul_core::parser::tree::{TreeFileFinder, TreeFmt};

use crate::cli::{InputCli, args::view::ViewTreeArgs, collect_paths};

use self::app::App;
use self::ui::TreePalette;

// How often to check for the tree count while it is being counted.
const COUNT_POLL: Duration = Duration::from_millis(200);

pub(in crate::cli) struct TreeViewParser<'a> {
    args: &'a ViewTreeArgs,
    input_dir: Option<PathBuf>,
}

impl InputCli for TreeViewParser<'_> {}

impl<'a> TreeViewParser<'a> {
    pub(in crate::cli) fn new(args: &'a ViewTreeArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn view(&mut self) {
        let input_fmt = self.args.input_fmt;
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        if files.is_empty() {
            eprintln!("Failed finding tree files to view");
            std::process::exit(1);
        }
        let view = match read_tree_at(&files[0], input_fmt, 0) {
            Ok(tree) if !tree.is_empty() => TreeView::new(tree),
            Ok(_) => {
                eprintln!("Failed finding a tree in {}", files[0].display());
                std::process::exit(1);
            }
            Err(e) => {
                eprintln!("{e:#}");
                std::process::exit(1);
            }
        };
        let theme = self.args.theme.parse::<Theme>().unwrap_or(Theme::Dark);
        let palette = TreePalette::new(theme, self.args.colormap);
        let mut app = App::new(&files, input_fmt, view);
        run(&mut app, &palette).expect("Failed running the tree viewer");
    }

    fn glob_paths(&self, dir: &str, input_fmt: &TreeFmt) -> Vec<PathBuf> {
        TreeFileFinder::new(Path::new(dir))
            .find(input_fmt)
            .unwrap_or_else(|err| {
                eprintln!("{err:#}");
                std::process::exit(1);
            })
    }
}

fn run(app: &mut App, palette: &TreePalette) -> io::Result<()> {
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
            app.poll_count();
            terminal.draw(|frame| ui::render(frame, app, palette))?;
            // Wake up now and then while trees are counted, to show the count.
            if app.is_counting() && !event::poll(COUNT_POLL)? {
                continue;
            }
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
