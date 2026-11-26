//! Viewer state and key bindings.
//!
//! Key bindings follow VS Code. Drawing lives in `ui`; this module only
//! changes state, so it can be tested without a terminal.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

use crossterm::event::{KeyCode, KeyEvent, KeyModifiers, MouseButton, MouseEvent, MouseEventKind};
use ratatui::layout::{Position, Rect};

use crate::core::align::view::{AlignmentView, PartitionTrack};
use crate::helper::types::{DataType, InputFmt};

use super::palette::ColorMode;

/// Cell widths, in terminal columns per site.
pub(super) const ZOOM_LEVELS: [usize; 3] = [1, 2, 3];
const DEFAULT_ZOOM: usize = 1;
const JUMP_SITES: usize = 10;
const MIN_NAME_WIDTH: u16 = 6;
const MAX_NAME_WIDTH: u16 = 60;
pub(super) const DEFAULT_NAME_WIDTH: u16 = 16;

/// Popups and prompts over the alignment.
#[derive(Debug, PartialEq)]
pub(super) enum Overlay {
    None,
    /// Go-to prompt, with the partition suggestion picked by ↑ or ↓.
    GoTo {
        input: String,
        index: Option<usize>,
    },
    Find(String),
    Files {
        query: String,
        index: usize,
    },
    Help,
}

/// What the sidebar lists.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum SidebarTab {
    Windows,
    Partitions,
}

/// Screen areas from the last draw, used to place mouse clicks.
#[derive(Debug, Default, Clone, Copy)]
pub(super) struct Areas {
    pub(super) names: Rect,
    pub(super) grid: Rect,
    pub(super) overview: Rect,
    pub(super) divider: Rect,
    pub(super) track: Rect,
    /// Inside of the sidebar border.
    pub(super) sidebar: Rect,
}

/// Lines per sidebar entry: label, statistics, and bar.
pub(super) const SIDEBAR_LINES: usize = 3;

pub(super) struct App<'a> {
    pub(super) files: &'a [PathBuf],
    pub(super) file_index: usize,
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    pub(super) view: AlignmentView,
    /// Partitions of a concatenated alignment, when found.
    pub(super) partitions: Option<PartitionTrack>,
    pub(super) row: usize,
    pub(super) site: usize,
    pub(super) top_row: usize,
    pub(super) left_site: usize,
    pub(super) zoom: usize,
    pub(super) color: ColorMode,
    pub(super) letters: bool,
    pub(super) window_size: usize,
    pub(super) selected: BTreeSet<usize>,
    anchor: Option<usize>,
    pub(super) show_sidebar: bool,
    pub(super) sidebar_tab: SidebarTab,
    pub(super) name_width: u16,
    pub(super) overlay: Overlay,
    /// Last search, kept for F3 after the prompt closes.
    pub(super) find_query: String,
    /// A short note shown in the status bar until the next key.
    pub(super) message: Option<String>,
    // Waiting for the second key of a Ctrl+K chord.
    is_chord: bool,
    resizing: bool,
    /// FASTA text to send to the terminal clipboard.
    pub(super) clipboard: Option<String>,
    pub(super) should_quit: bool,
    /// Rows and columns of the grid in the last draw.
    pub(super) grid_rows: usize,
    pub(super) grid_cols: usize,
    pub(super) areas: Areas,
}

impl<'a> App<'a> {
    pub(super) fn new(
        files: &'a [PathBuf],
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        view: AlignmentView,
        window_size: usize,
    ) -> Self {
        let (partitions, message) = detect_partitions(&files[0], &view);
        Self {
            files,
            file_index: 0,
            input_fmt,
            datatype,
            view,
            sidebar_tab: if partitions.is_some() {
                SidebarTab::Partitions
            } else {
                SidebarTab::Windows
            },
            partitions,
            row: 0,
            site: 0,
            top_row: 0,
            left_site: 0,
            zoom: DEFAULT_ZOOM,
            color: ColorMode::Base,
            letters: true,
            window_size,
            selected: BTreeSet::new(),
            anchor: None,
            show_sidebar: false,
            name_width: DEFAULT_NAME_WIDTH,
            overlay: Overlay::None,
            find_query: String::new(),
            message,
            is_chord: false,
            resizing: false,
            clipboard: None,
            should_quit: false,
            grid_rows: 1,
            grid_cols: 1,
            areas: Areas::default(),
        }
    }

    pub(super) fn path(&self) -> &Path {
        &self.files[self.file_index]
    }

    /// Sites that fit across the grid.
    pub(super) fn grid_sites(&self) -> usize {
        (self.grid_cols / self.zoom).max(1)
    }

    pub(super) fn window_count(&self) -> usize {
        self.view.length().div_ceil(self.window_size).max(1)
    }

    pub(super) fn window(&self) -> usize {
        self.site / self.window_size
    }

    /// Sites the cursor and view may reach.
    pub(super) fn site_range(&self) -> (usize, usize) {
        (0, self.view.length().max(1))
    }

    /// Partition blocks suggested for go-to `input`: those whose name holds
    /// the text before any `:`. Exact and prefix matches come first.
    pub(super) fn goto_matches(&self, input: &str) -> Vec<usize> {
        let Some(track) = &self.partitions else {
            return Vec::new();
        };
        let name = input.split(':').next().unwrap_or("").trim().to_lowercase();
        if name.is_empty() || name.chars().all(|c| c.is_ascii_digit()) {
            return Vec::new();
        }
        let mut matches: Vec<(u8, usize)> = track
            .blocks
            .iter()
            .enumerate()
            .filter_map(|(i, b)| {
                let block = b.name.to_lowercase();
                let rank = if block == name {
                    0
                } else if block.starts_with(&name) {
                    1
                } else if block.contains(&name) {
                    2
                } else {
                    return None;
                };
                Some((rank, i))
            })
            .collect();
        matches.sort_unstable();
        matches.into_iter().map(|(_, i)| i).collect()
    }

    /// Partition block under the cursor.
    pub(super) fn block(&self) -> Option<usize> {
        self.partitions.as_ref()?.block_at(self.site)
    }

    /// Entries the sidebar lists: partition blocks or windows.
    pub(super) fn sidebar_count(&self) -> usize {
        match (&self.partitions, self.sidebar_tab) {
            (Some(track), SidebarTab::Partitions) => track.len(),
            _ => self.window_count(),
        }
    }

    /// First sidebar entry shown when `shown` fit, keeping the current one
    /// near the top.
    pub(super) fn sidebar_first(&self, shown: usize) -> usize {
        let current = match (&self.partitions, self.sidebar_tab) {
            (Some(track), SidebarTab::Partitions) => self
                .block()
                .or_else(|| track.prev_block(self.site))
                .unwrap_or(0),
            _ => self.window(),
        };
        current
            .saturating_sub(shown / 3)
            .min(self.sidebar_count().saturating_sub(shown))
    }

    /// Files whose names contain the picker query, as indices into `files`.
    pub(super) fn matching_files(&self, query: &str) -> Vec<usize> {
        let query = query.to_lowercase();
        self.files
            .iter()
            .enumerate()
            .filter(|(_, f)| f.to_string_lossy().to_lowercase().contains(&query))
            .map(|(i, _)| i)
            .collect()
    }

    /// Rows whose names contain the search text.
    pub(super) fn find_matches(&self, query: &str) -> Vec<usize> {
        if query.is_empty() {
            return Vec::new();
        }
        let query = query.to_lowercase();
        self.view
            .names()
            .iter()
            .enumerate()
            .filter(|(_, n)| n.to_lowercase().contains(&query))
            .map(|(i, _)| i)
            .collect()
    }

    pub(super) fn handle_key(&mut self, key: KeyEvent) {
        self.message = None;
        match self.overlay {
            Overlay::None => self.handle_view_key(key),
            Overlay::Help => {
                if matches!(key.code, KeyCode::Esc | KeyCode::F(1) | KeyCode::Enter)
                    || is_ctrl(&key, 'q')
                {
                    self.overlay = Overlay::None;
                }
            }
            Overlay::GoTo { .. } => self.handle_goto_key(key),
            Overlay::Find(_) => self.handle_find_key(key),
            Overlay::Files { .. } => self.handle_files_key(key),
        }
    }

    fn handle_view_key(&mut self, key: KeyEvent) {
        let ctrl = key.modifiers.contains(KeyModifiers::CONTROL);
        let alt = key.modifiers.contains(KeyModifiers::ALT);
        let shift = key.modifiers.contains(KeyModifiers::SHIFT);
        if self.is_chord {
            self.is_chord = false;
            match key.code {
                KeyCode::Char('c' | 'C') => self.color = self.color.next(),
                KeyCode::Char('l' | 'L') => self.letters = !self.letters,
                KeyCode::Char('t' | 'T') => self.toggle_sidebar_tab(),
                _ => self.message = Some(String::from("Unknown chord")),
            }
            return;
        }
        let page = self.grid_rows.saturating_sub(1).max(1);
        match key.code {
            KeyCode::Esc => self.should_quit = true,
            KeyCode::Char('q') if ctrl => self.should_quit = true,
            KeyCode::Char('k') if ctrl => {
                self.is_chord = true;
                self.message = Some(String::from("Ctrl+K: C color · L letters · T sidebar tab"));
            }
            KeyCode::Char('g') if ctrl => {
                self.overlay = Overlay::GoTo {
                    input: String::new(),
                    index: None,
                }
            }
            KeyCode::Char('f') if ctrl => self.overlay = Overlay::Find(String::new()),
            KeyCode::Char('p') if ctrl => {
                self.overlay = Overlay::Files {
                    query: String::new(),
                    index: self.file_index,
                }
            }
            KeyCode::Char('b') if ctrl => self.show_sidebar = !self.show_sidebar,
            KeyCode::Char('a') if ctrl => self.selected = (0..self.view.rows()).collect(),
            KeyCode::Char('c') if ctrl => self.copy_selection(),
            KeyCode::Char('s') if ctrl => self.save_selection(),
            KeyCode::Char('0') if ctrl => self.set_zoom(DEFAULT_ZOOM),
            KeyCode::Char('=' | '+') => self.zoom_by(1),
            KeyCode::Char('-' | '_') => self.zoom_by(-1),
            KeyCode::Char(' ') => self.toggle_row(self.row),
            // Option+← and Option+→ in macOS terminals that send Esc+b and Esc+f.
            KeyCode::Char('b') if alt => self.move_site(-(JUMP_SITES as isize)),
            KeyCode::Char('f') if alt => self.move_site(JUMP_SITES as isize),
            KeyCode::Char('[') => self.step_window(-1),
            KeyCode::Char(']') => self.step_window(1),
            KeyCode::Char('{') => self.step_partition(-1),
            KeyCode::Char('}') => self.step_partition(1),
            KeyCode::Char('<' | ',') => self.resize_names(-2),
            KeyCode::Char('>' | '.') => self.resize_names(2),
            KeyCode::Char('?') | KeyCode::F(1) => self.overlay = Overlay::Help,
            KeyCode::F(2) => self.color = self.color.next(),
            KeyCode::F(3) if shift => self.find_next(false),
            KeyCode::F(3) => self.find_next(true),
            KeyCode::F(4) => self.letters = !self.letters,
            KeyCode::Left if ctrl && shift => self.resize_names(-2),
            KeyCode::Right if ctrl && shift => self.resize_names(2),
            KeyCode::Left if alt && shift => self.move_site(-(self.grid_sites() as isize)),
            KeyCode::Right if alt && shift => self.move_site(self.grid_sites() as isize),
            // Ctrl+arrows on Windows and Linux, Option+arrows on macOS.
            KeyCode::Left if ctrl || alt => self.move_site(-(JUMP_SITES as isize)),
            KeyCode::Right if ctrl || alt => self.move_site(JUMP_SITES as isize),
            KeyCode::Left => self.move_site(-1),
            KeyCode::Right => self.move_site(1),
            KeyCode::Up if shift => self.extend_selection(-1),
            KeyCode::Down if shift => self.extend_selection(1),
            // Option+↑ and Option+↓, for Mac keyboards without page keys.
            KeyCode::Up if alt => self.move_row(-(page as isize)),
            KeyCode::Down if alt => self.move_row(page as isize),
            KeyCode::Up => self.move_row(-1),
            KeyCode::Down => self.move_row(1),
            KeyCode::Home if ctrl => self.set_cursor(0, self.site_range().0),
            KeyCode::End if ctrl => {
                let last_site = self.site_range().1 - 1;
                self.set_cursor(self.view.rows().saturating_sub(1), last_site)
            }
            KeyCode::Home => self.set_cursor(self.row, self.site_range().0),
            KeyCode::End => self.set_cursor(self.row, self.row_end()),
            KeyCode::PageUp if ctrl => self.step_file(-1),
            KeyCode::PageDown if ctrl => self.step_file(1),
            KeyCode::PageUp if alt => self.move_site(-(self.grid_sites() as isize)),
            KeyCode::PageDown if alt => self.move_site(self.grid_sites() as isize),
            KeyCode::PageUp => self.move_row(-(page as isize)),
            KeyCode::PageDown => self.move_row(page as isize),
            _ => (),
        }
    }

    fn handle_goto_key(&mut self, key: KeyEvent) {
        let matches = match &self.overlay {
            Overlay::GoTo { input, .. } => self.goto_matches(input),
            _ => return,
        };
        let Overlay::GoTo { input, index } = &mut self.overlay else {
            return;
        };
        let count = matches.len();
        match key.code {
            KeyCode::Esc => self.overlay = Overlay::None,
            KeyCode::Up if count > 0 => {
                *index = Some(index.map_or(count - 1, |i| (i + count - 1) % count));
            }
            KeyCode::Down if count > 0 => {
                *index = Some(index.map_or(0, |i| (i + 1) % count));
            }
            // Fill in the picked or first suggestion, ready for `:site`.
            KeyCode::Tab if count > 0 => {
                let block = matches[index.unwrap_or(0).min(count - 1)];
                if let Some(track) = &self.partitions {
                    *input = track.blocks[block].name.clone();
                }
                *index = None;
            }
            KeyCode::Enter => {
                let picked = index.and_then(|i| matches.get(i).copied());
                let input = input.trim().to_string();
                self.overlay = Overlay::None;
                match picked {
                    Some(block) => self.jump_to_block(block),
                    None => {
                        if let Err(msg) = self.go_to(&input) {
                            self.message = Some(msg);
                        }
                    }
                }
            }
            KeyCode::Backspace => {
                input.pop();
                *index = None;
            }
            KeyCode::Char(c) if !key.modifiers.contains(KeyModifiers::CONTROL) => {
                input.push(c);
                *index = None;
            }
            _ => (),
        }
    }

    fn handle_find_key(&mut self, key: KeyEvent) {
        let Overlay::Find(input) = &mut self.overlay else {
            return;
        };
        let shift = key.modifiers.contains(KeyModifiers::SHIFT);
        match key.code {
            KeyCode::Esc => self.overlay = Overlay::None,
            KeyCode::Enter | KeyCode::F(3) => self.find_next(!shift),
            KeyCode::Backspace => {
                input.pop();
                self.find_query = input.clone();
            }
            KeyCode::Char(c) if !key.modifiers.contains(KeyModifiers::CONTROL) => {
                input.push(c);
                self.find_query = input.clone();
                // Search as you type, from the current row.
                if let Some(&row) = self
                    .find_matches(&self.find_query)
                    .iter()
                    .find(|&&r| r >= self.row)
                    .or(self.find_matches(&self.find_query).first())
                {
                    self.set_cursor(row, self.site);
                }
            }
            _ => (),
        }
    }

    fn handle_files_key(&mut self, key: KeyEvent) {
        let matches = match &self.overlay {
            Overlay::Files { query, .. } => self.matching_files(query),
            _ => return,
        };
        let Overlay::Files { query, index } = &mut self.overlay else {
            return;
        };
        let position = matches.iter().position(|i| i == index);
        match key.code {
            KeyCode::Esc => self.overlay = Overlay::None,
            KeyCode::Up => {
                if let Some(p) = position.filter(|&p| p > 0) {
                    *index = matches[p - 1];
                } else if let Some(&last) = matches.last() {
                    *index = last;
                }
            }
            KeyCode::Down => match position {
                Some(p) if p + 1 < matches.len() => *index = matches[p + 1],
                _ => {
                    if let Some(&first) = matches.first() {
                        *index = first
                    }
                }
            },
            KeyCode::Enter => {
                let target = position.map(|_| *index).or(matches.first().copied());
                self.overlay = Overlay::None;
                if let Some(target) = target {
                    self.open_file(target);
                }
            }
            KeyCode::Backspace => {
                query.pop();
            }
            KeyCode::Char(c) if !key.modifiers.contains(KeyModifiers::CONTROL) => query.push(c),
            _ => (),
        }
    }

    pub(super) fn handle_mouse(&mut self, mouse: MouseEvent) {
        if self.overlay != Overlay::None {
            return;
        }
        let ctrl = mouse.modifiers.contains(KeyModifiers::CONTROL);
        let shift = mouse.modifiers.contains(KeyModifiers::SHIFT);
        let at = Position::new(mouse.column, mouse.row);
        match mouse.kind {
            MouseEventKind::ScrollUp if ctrl => self.zoom_by(1),
            MouseEventKind::ScrollDown if ctrl => self.zoom_by(-1),
            MouseEventKind::ScrollUp if shift => self.scroll_sites(-(JUMP_SITES as isize)),
            MouseEventKind::ScrollDown if shift => self.scroll_sites(JUMP_SITES as isize),
            MouseEventKind::ScrollLeft => self.scroll_sites(-(JUMP_SITES as isize)),
            MouseEventKind::ScrollRight => self.scroll_sites(JUMP_SITES as isize),
            MouseEventKind::ScrollUp => self.scroll_rows(-3),
            MouseEventKind::ScrollDown => self.scroll_rows(3),
            MouseEventKind::Down(MouseButton::Left) => {
                if self.areas.divider.contains(at) {
                    self.resizing = true;
                } else if self.areas.overview.contains(at) {
                    self.jump_overview(mouse.column);
                } else if self.areas.track.contains(at) {
                    self.click_track(at);
                } else if self.areas.sidebar.contains(at) {
                    self.click_sidebar(at);
                } else if self.areas.grid.contains(at) {
                    if let Some((row, site)) = self.cell_at(at) {
                        self.set_cursor(row, site);
                    }
                } else if self.areas.names.contains(at) {
                    self.click_name(at, ctrl, shift);
                }
            }
            MouseEventKind::Drag(MouseButton::Left) => {
                if self.resizing {
                    let width = mouse.column.saturating_sub(self.areas.names.x);
                    self.name_width = width.clamp(MIN_NAME_WIDTH, MAX_NAME_WIDTH);
                } else if self.areas.overview.contains(at) {
                    self.jump_overview(mouse.column);
                }
            }
            MouseEventKind::Up(_) => self.resizing = false,
            _ => (),
        }
    }

    /// Goes to `123` (a site), a partition, `partition:12` (a site within
    /// it), a taxon, or `taxon:123`. Partitions win over taxa of the same
    /// name; text matching only part of a partition name goes to the first
    /// such partition when no taxon matches.
    pub(super) fn go_to(&mut self, input: &str) -> Result<(), String> {
        if input.is_empty() {
            return Err(String::from("Enter a site, partition, or taxon."));
        }
        let (name, local) = match input.rsplit_once(':') {
            Some((name, site)) => (name.trim(), Some(site.trim())),
            None => (input, None),
        };
        if let Some(block) = self.partitions.as_ref().and_then(|t| t.find(name)) {
            return match local {
                Some(site) => self.go_to_partition_site(block, site),
                None => {
                    self.jump_to_block(block);
                    Ok(())
                }
            };
        }
        let length = self.view.length();
        let site_error = || format!("Enter a site from 1 to {}.", length);
        let (name, site) = match input.rsplit_once(':') {
            Some((name, site)) => (Some(name.trim()), Some(site.trim())),
            None if input.chars().all(|c| c.is_ascii_digit()) => (None, Some(input)),
            None => (Some(input), None),
        };
        let row = match name.filter(|n| !n.is_empty()) {
            Some(name) => match self.find_matches(name).first() {
                Some(&row) => row,
                None => {
                    if local.is_none()
                        && let Some(&block) = self.goto_matches(name).first()
                    {
                        self.jump_to_block(block);
                        return Ok(());
                    }
                    return Err(format!("No taxon or partition matches \"{}\".", name));
                }
            },
            None => self.row,
        };
        let site = match site {
            Some(site) => match site.parse::<usize>() {
                Ok(s) if (1..=length).contains(&s) => s - 1,
                _ => return Err(site_error()),
            },
            None => self.site,
        };
        self.jump_to(row, site);
        self.center_on_cursor();
        Ok(())
    }

    fn find_next(&mut self, forward: bool) {
        let matches = self.find_matches(&self.find_query);
        if matches.is_empty() {
            self.message = Some(if self.find_query.is_empty() {
                String::from("Press Ctrl+F to find a taxon.")
            } else {
                format!("No taxon matches \"{}\".", self.find_query)
            });
            return;
        }
        let next = if forward {
            matches.iter().find(|&&r| r > self.row).or(matches.first())
        } else {
            matches
                .iter()
                .rev()
                .find(|&&r| r < self.row)
                .or(matches.last())
        };
        if let Some(&row) = next {
            self.set_cursor(row, self.site);
        }
    }

    pub(super) fn move_row(&mut self, delta: isize) {
        let row = offset(self.row, delta, self.view.rows().saturating_sub(1));
        self.set_cursor(row, self.site);
    }

    // Goes to one-based `site` counted from the start of partition `block`.
    fn go_to_partition_site(&mut self, block: usize, site: &str) -> Result<(), String> {
        let Some(track) = &self.partitions else {
            return Ok(());
        };
        let block = &track.blocks[block];
        let sites = block.end - block.start;
        match site.parse::<usize>() {
            Ok(s) if (1..=sites).contains(&s) => {
                let target = block.start + s - 1;
                self.jump_to(self.row, target);
                self.center_on_cursor();
                Ok(())
            }
            _ => Err(format!(
                "Enter a site from 1 to {} in {}.",
                sites, block.name
            )),
        }
    }

    pub(super) fn move_site(&mut self, delta: isize) {
        let (start, end) = self.site_range();
        let site = offset(self.site, delta, end - 1).max(start);
        self.set_cursor(self.row, site);
    }

    /// Moves the cursor anywhere.
    fn jump_to(&mut self, row: usize, site: usize) {
        self.set_cursor(row, site);
    }

    /// Moves the cursor within the reachable sites and scrolls it into view.
    pub(super) fn set_cursor(&mut self, row: usize, site: usize) {
        let (start, end) = self.site_range();
        self.row = row.min(self.view.rows().saturating_sub(1));
        self.site = site.clamp(start, end.saturating_sub(1).max(start));
        self.anchor = None;
        self.scroll_to_cursor();
    }

    fn scroll_to_cursor(&mut self) {
        let rows = self.grid_rows.max(1);
        if self.row < self.top_row {
            self.top_row = self.row;
        } else if self.row >= self.top_row + rows {
            self.top_row = self.row + 1 - rows;
        }
        let sites = self.grid_sites();
        if self.site < self.left_site {
            self.left_site = self.site;
        } else if self.site >= self.left_site + sites {
            self.left_site = self.site + 1 - sites;
        }
        self.clamp_scroll();
    }

    fn center_on_cursor(&mut self) {
        self.top_row = self.row.saturating_sub(self.grid_rows / 2);
        self.left_site = self.site.saturating_sub(self.grid_sites() / 2);
        self.clamp_scroll();
    }

    /// Keeps the view inside the alignment.
    pub(super) fn clamp_scroll(&mut self) {
        let max_top = self.view.rows().saturating_sub(self.grid_rows.max(1));
        self.top_row = self.top_row.min(max_top);
        let (start, end) = self.site_range();
        let max_left = end.saturating_sub(self.grid_sites()).max(start);
        self.left_site = self.left_site.clamp(start, max_left);
    }

    fn scroll_rows(&mut self, delta: isize) {
        self.top_row = offset(self.top_row, delta, usize::MAX);
        self.clamp_scroll();
    }

    fn scroll_sites(&mut self, delta: isize) {
        self.left_site = offset(self.left_site, delta, usize::MAX);
        self.clamp_scroll();
    }

    fn row_end(&self) -> usize {
        let (start, end) = self.site_range();
        let len = self.view.row_len(self.row).max(1);
        (len - 1).clamp(start, end - 1)
    }

    pub(super) fn zoom_by(&mut self, step: isize) {
        let index = ZOOM_LEVELS
            .iter()
            .position(|&z| z == self.zoom)
            .unwrap_or(0);
        let index = offset(index, step, ZOOM_LEVELS.len() - 1);
        self.set_zoom(ZOOM_LEVELS[index]);
    }

    // Changes the cell width, keeping the center site in place.
    fn set_zoom(&mut self, zoom: usize) {
        let center = self.left_site + self.grid_sites() / 2;
        self.zoom = zoom;
        self.left_site = center.saturating_sub(self.grid_sites() / 2);
        self.clamp_scroll();
    }

    /// Moves to the start of the previous or next window.
    pub(super) fn step_window(&mut self, step: isize) {
        let window = offset(self.window(), step, self.window_count() - 1);
        let site = window * self.window_size;
        self.jump_to(self.row, site);
        self.left_site = site;
        self.clamp_scroll();
    }

    /// Moves to the start of the previous or next partition.
    fn step_partition(&mut self, step: isize) {
        let Some(track) = &self.partitions else {
            self.message = Some(String::from("No partitions found"));
            return;
        };
        let target = if step < 0 {
            track.prev_block(self.site)
        } else {
            track.next_block(self.site)
        };
        match target {
            Some(block) => self.jump_to_block(block),
            None if step < 0 => self.message = Some(String::from("At the first partition")),
            None => self.message = Some(String::from("At the last partition")),
        }
    }

    /// Moves the cursor to the first site of partition block `index`.
    pub(super) fn jump_to_block(&mut self, index: usize) {
        let Some(start) = self
            .partitions
            .as_ref()
            .and_then(|t| t.blocks.get(index))
            .map(|b| b.start)
        else {
            return;
        };
        self.jump_to(self.row, start);
        self.left_site = start;
        self.clamp_scroll();
    }

    fn toggle_sidebar_tab(&mut self) {
        if self.partitions.is_none() {
            self.message = Some(String::from("No partitions found"));
            return;
        }
        self.sidebar_tab = match self.sidebar_tab {
            SidebarTab::Windows => SidebarTab::Partitions,
            SidebarTab::Partitions => SidebarTab::Windows,
        };
        self.show_sidebar = true;
    }

    fn click_track(&mut self, at: Position) {
        let site = self.left_site + (at.x - self.areas.track.x) as usize / self.zoom;
        if let Some(block) = self.partitions.as_ref().and_then(|t| t.block_at(site)) {
            self.jump_to_block(block);
        }
    }

    fn click_sidebar(&mut self, at: Position) {
        let area = self.areas.sidebar;
        let shown = (area.height as usize / SIDEBAR_LINES).max(1);
        let index = self.sidebar_first(shown) + (at.y - area.y) as usize / SIDEBAR_LINES;
        if index >= self.sidebar_count() {
            return;
        }
        match (&self.partitions, self.sidebar_tab) {
            (Some(_), SidebarTab::Partitions) => self.jump_to_block(index),
            _ => {
                let site = index * self.window_size;
                self.jump_to(self.row, site);
                self.left_site = site;
                self.clamp_scroll();
            }
        }
    }

    fn resize_names(&mut self, delta: i16) {
        self.name_width = self
            .name_width
            .saturating_add_signed(delta)
            .clamp(MIN_NAME_WIDTH, MAX_NAME_WIDTH);
    }

    fn toggle_row(&mut self, row: usize) {
        if !self.selected.remove(&row) {
            self.selected.insert(row);
        }
    }

    // Shift+Up/Down: select from the anchor row to the cursor.
    fn extend_selection(&mut self, delta: isize) {
        let anchor = self.anchor.unwrap_or(self.row);
        self.selected.insert(self.row);
        self.move_row(delta);
        self.anchor = Some(anchor);
        let (from, to) = (anchor.min(self.row), anchor.max(self.row));
        self.selected.extend(from..=to);
    }

    fn click_name(&mut self, at: Position, ctrl: bool, shift: bool) {
        let Some((row, _)) = self.cell_at(Position::new(self.areas.grid.x, at.y)) else {
            return;
        };
        if ctrl {
            self.toggle_row(row);
        } else if shift {
            let (from, to) = (self.row.min(row), self.row.max(row));
            self.selected.extend(from..=to);
        } else {
            self.selected = BTreeSet::from([row]);
        }
        self.set_cursor(row, self.site);
    }

    // Grid row and site under a screen position.
    fn cell_at(&self, at: Position) -> Option<(usize, usize)> {
        let grid = self.areas.grid;
        if at.y < grid.y || at.x < grid.x {
            return None;
        }
        let row = self.top_row + (at.y - grid.y) as usize;
        let site = self.left_site + (at.x - grid.x) as usize / self.zoom;
        (row < self.view.rows() && site < self.site_range().1).then_some((row, site))
    }

    fn jump_overview(&mut self, column: u16) {
        let area = self.areas.overview;
        if area.width == 0 {
            return;
        }
        let x = column.saturating_sub(area.x) as usize;
        let site = x * self.view.length() / area.width as usize;
        self.jump_to(self.row, site);
        self.center_on_cursor();
    }

    // Rows to copy or save: the selection, or the cursor row.
    fn chosen_rows(&self) -> Vec<usize> {
        if self.selected.is_empty() {
            vec![self.row]
        } else {
            self.selected.iter().copied().collect()
        }
    }

    fn copy_selection(&mut self) {
        let rows = self.chosen_rows();
        self.clipboard = Some(self.view.fasta(&rows));
        self.message = Some(format!("Copied {} sequence(s) as FASTA", rows.len()));
    }

    fn save_selection(&mut self) {
        let rows = self.chosen_rows();
        let stem = self
            .path()
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| String::from("alignment"));
        let output = PathBuf::from(format!("{}-selected.fas", stem));
        if output.exists() {
            self.message = Some(format!("{} already exists", output.display()));
            return;
        }
        self.message = Some(match fs::write(&output, self.view.fasta(&rows)) {
            Ok(()) => format!("Saved {} sequence(s) to {}", rows.len(), output.display()),
            Err(e) => format!("Failed saving {}: {}", output.display(), e),
        });
    }

    fn step_file(&mut self, step: isize) {
        if self.files.len() < 2 {
            self.message = Some(String::from("Only one file is open"));
            return;
        }
        let count = self.files.len() as isize;
        let index = (self.file_index as isize + step).rem_euclid(count) as usize;
        self.open_file(index);
    }

    /// Opens another input file, keeping the view settings.
    pub(super) fn open_file(&mut self, index: usize) {
        if index == self.file_index {
            return;
        }
        match AlignmentView::open(&self.files[index], self.input_fmt, self.datatype) {
            Ok(view) => {
                let (partitions, message) = detect_partitions(&self.files[index], &view);
                if partitions.is_some() {
                    self.sidebar_tab = SidebarTab::Partitions;
                } else {
                    self.sidebar_tab = SidebarTab::Windows;
                }
                self.partitions = partitions;
                self.message = message;
                self.view = view;
                self.file_index = index;
                self.selected.clear();
                self.anchor = None;
                self.row = 0;
                self.site = 0;
                self.top_row = 0;
                self.left_site = 0;
            }
            Err(e) => self.message = Some(e.to_string()),
        }
    }
}

// Partitions of `path`, or a note when its partition file cannot be read.
fn detect_partitions(
    path: &Path,
    view: &AlignmentView,
) -> (Option<PartitionTrack>, Option<String>) {
    match PartitionTrack::detect(path, view.length()) {
        Ok(track) => (track, None),
        Err(e) => (None, Some(format!("Failed reading partitions: {e}"))),
    }
}

fn is_ctrl(key: &KeyEvent, c: char) -> bool {
    key.modifiers.contains(KeyModifiers::CONTROL) && key.code == KeyCode::Char(c)
}

// `value + delta`, kept within 0 and `max`.
fn offset(value: usize, delta: isize, max: usize) -> usize {
    value.saturating_add_signed(delta).min(max)
}

#[cfg(test)]
mod test {
    use super::*;

    const FILE: &str = "crates/segul-core/tests/files/simple.fas";

    fn key(code: KeyCode) -> KeyEvent {
        KeyEvent::new(code, KeyModifiers::NONE)
    }

    fn with(code: KeyCode, modifiers: KeyModifiers) -> KeyEvent {
        KeyEvent::new(code, modifiers)
    }

    fn app<'a>(files: &'a [PathBuf]) -> App<'a> {
        let view = AlignmentView::open(&files[0], &InputFmt::Auto, &DataType::Ignore)
            .expect("Failed opening fixture");
        let mut app = App::new(files, &InputFmt::Auto, &DataType::Ignore, view, 2);
        app.grid_rows = 2;
        app.grid_cols = 3;
        app
    }

    #[test]
    fn test_cursor_moves_and_clamps() {
        let files = [PathBuf::from(FILE)];
        let mut app = app(&files);
        app.handle_key(key(KeyCode::Left));
        app.handle_key(key(KeyCode::Up));
        assert_eq!((app.row, app.site), (0, 0));
        app.handle_key(key(KeyCode::Down));
        app.handle_key(key(KeyCode::Down));
        assert_eq!(app.row, 1);
        app.handle_key(key(KeyCode::End));
        assert_eq!(app.site, 5);
        // The view scrolls to keep the cursor visible.
        assert_eq!(app.left_site, 3);
        app.handle_key(with(KeyCode::Home, KeyModifiers::CONTROL));
        assert_eq!((app.row, app.site, app.left_site), (0, 0, 0));
        app.handle_key(with(KeyCode::Right, KeyModifiers::CONTROL));
        assert_eq!(app.site, 5);
    }

    #[test]
    fn test_mac_keys() {
        let files = [PathBuf::from(FILE)];
        let mut app = app(&files);
        // Option+→ as an Alt arrow and as Esc+f.
        app.handle_key(with(KeyCode::Right, KeyModifiers::ALT));
        assert_eq!(app.site, 5);
        app.handle_key(with(KeyCode::Char('b'), KeyModifiers::ALT));
        assert_eq!(app.site, 0);
        app.handle_key(with(KeyCode::Down, KeyModifiers::ALT));
        assert_eq!(app.row, 1);
        app.handle_key(key(KeyCode::Char('>')));
        assert_eq!(app.name_width, DEFAULT_NAME_WIDTH + 2);
        app.handle_key(key(KeyCode::Char('?')));
        assert_eq!(app.overlay, Overlay::Help);
    }

    #[test]
    fn test_go_to() {
        let files = [PathBuf::from(FILE)];
        let mut app = app(&files);
        assert!(app.go_to("4").is_ok());
        assert_eq!(app.site, 3);
        assert!(app.go_to("cdef:2").is_ok());
        assert_eq!((app.row, app.site), (1, 1));
        assert!(app.go_to("abcd").is_ok());
        assert_eq!(app.row, 0);
        assert!(app.go_to("0").is_err());
        assert!(app.go_to("7").is_err());
        assert!(app.go_to("zzz:1").is_err());
    }

    #[test]
    fn test_go_to_prompt() {
        let files = [PathBuf::from(FILE)];
        let mut app = app(&files);
        app.handle_key(with(KeyCode::Char('g'), KeyModifiers::CONTROL));
        app.handle_key(key(KeyCode::Char('9')));
        app.handle_key(key(KeyCode::Enter));
        assert_eq!(app.overlay, Overlay::None);
        assert_eq!(app.message.as_deref(), Some("Enter a site from 1 to 6."));
    }

    #[test]
    fn test_find_wraps() {
        let files = [PathBuf::from(FILE)];
        let mut app = app(&files);
        app.handle_key(with(KeyCode::Char('f'), KeyModifiers::CONTROL));
        app.handle_key(key(KeyCode::Char('d')));
        // Both names contain "d"; typing keeps the current match.
        assert_eq!(app.row, 0);
        app.handle_key(key(KeyCode::Enter));
        assert_eq!(app.row, 1);
        app.handle_key(key(KeyCode::Enter));
        assert_eq!(app.row, 0);
        app.handle_key(key(KeyCode::Esc));
        app.handle_key(with(KeyCode::F(3), KeyModifiers::SHIFT));
        assert_eq!(app.row, 1);
    }

    #[test]
    fn test_zoom_keeps_center() {
        let files = [PathBuf::from(FILE)];
        let mut app = app(&files);
        app.grid_cols = 2;
        app.left_site = 3;
        // Two sites wide at zoom 1: center site 4.
        app.handle_key(key(KeyCode::Char('+')));
        assert_eq!(app.zoom, 2);
        assert_eq!(app.left_site, 4);
        app.handle_key(with(KeyCode::Char('0'), KeyModifiers::CONTROL));
        assert_eq!(app.zoom, DEFAULT_ZOOM);
    }

    #[test]
    fn test_windows() {
        let files = [PathBuf::from(FILE)];
        let mut app = app(&files);
        assert_eq!(app.window_count(), 3);
        app.handle_key(key(KeyCode::Char(']')));
        assert_eq!((app.window(), app.site), (1, 2));
        app.handle_key(key(KeyCode::Char(']')));
        assert_eq!(app.site, 4);
        app.handle_key(key(KeyCode::Char('[')));
        app.handle_key(key(KeyCode::Char('[')));
        assert_eq!(app.site, 0);
    }

    #[test]
    fn test_selection() {
        let files = [PathBuf::from(FILE)];
        let mut app = app(&files);
        app.handle_key(with(KeyCode::Down, KeyModifiers::SHIFT));
        assert_eq!(app.selected, BTreeSet::from([0, 1]));
        app.handle_key(key(KeyCode::Char(' ')));
        assert_eq!(app.selected, BTreeSet::from([0]));
        app.handle_key(with(KeyCode::Char('c'), KeyModifiers::CONTROL));
        assert_eq!(app.clipboard.as_deref(), Some(">ABCD\nAGTATG\n"));
    }

    #[test]
    fn test_switch_files() {
        let files = [
            PathBuf::from(FILE),
            PathBuf::from("crates/segul-core/tests/files/simple_aa.fas"),
        ];
        let mut app = app(&files);
        app.handle_key(with(KeyCode::PageDown, KeyModifiers::CONTROL));
        assert_eq!(app.file_index, 1);
        assert!(!app.view.is_dna());
        app.handle_key(with(KeyCode::PageDown, KeyModifiers::CONTROL));
        assert_eq!(app.file_index, 0);
    }

    #[test]
    fn test_partitions() {
        let files = [PathBuf::from(
            "crates/segul-core/tests/files/view/concat_charset.nex",
        )];
        let mut app = app(&files);
        assert_eq!(app.sidebar_tab, SidebarTab::Partitions);
        assert_eq!(app.block(), Some(0));
        app.handle_key(key(KeyCode::Char('}')));
        assert_eq!((app.site, app.block()), (5, Some(1)));
        app.handle_key(key(KeyCode::Right));
        app.handle_key(key(KeyCode::Char('{')));
        assert_eq!(app.site, 5);
        app.handle_key(key(KeyCode::Char('{')));
        assert_eq!(app.site, 0);
        app.handle_key(key(KeyCode::Char('{')));
        assert_eq!(app.message.as_deref(), Some("At the first partition"));
        assert!(app.go_to("locus3").is_ok());
        assert_eq!(app.site, 10);
        assert!(app.go_to("LOCUS2:3").is_ok());
        assert_eq!(app.site, 7);
        assert_eq!(
            app.go_to("locus2:9").unwrap_err(),
            "Enter a site from 1 to 5 in locus2."
        );
        // Part of a name goes to the first match when no taxon matches.
        assert!(app.go_to("cus3").is_ok());
        assert_eq!(app.site, 10);
        // A taxon still wins over a partial partition match.
        assert!(app.go_to("abce").is_ok());
        assert_eq!(app.row, 1);
        app.handle_key(with(KeyCode::Char('k'), KeyModifiers::CONTROL));
        app.handle_key(key(KeyCode::Char('t')));
        assert_eq!(app.sidebar_tab, SidebarTab::Windows);
        assert!(app.show_sidebar);
    }

    #[test]
    fn test_go_to_suggestions() {
        let files = [PathBuf::from(
            "crates/segul-core/tests/files/view/concat_charset.nex",
        )];
        let mut app = app(&files);
        assert_eq!(app.goto_matches("locus2"), vec![1]);
        assert_eq!(app.goto_matches("2"), Vec::<usize>::new());
        assert_eq!(app.goto_matches("cus"), vec![0, 1, 2]);
        assert!(app.goto_matches("12").is_empty());
        app.handle_key(with(KeyCode::Char('g'), KeyModifiers::CONTROL));
        for c in "loc".chars() {
            app.handle_key(key(KeyCode::Char(c)));
        }
        app.handle_key(key(KeyCode::Down));
        app.handle_key(key(KeyCode::Down));
        app.handle_key(key(KeyCode::Enter));
        assert_eq!(app.overlay, Overlay::None);
        assert_eq!(app.block(), Some(1));
        // Tab fills the name, then a site within the partition.
        app.handle_key(with(KeyCode::Char('g'), KeyModifiers::CONTROL));
        for c in "us3".chars() {
            app.handle_key(key(KeyCode::Char(c)));
        }
        app.handle_key(key(KeyCode::Tab));
        app.handle_key(key(KeyCode::Char(':')));
        app.handle_key(key(KeyCode::Char('2')));
        app.handle_key(key(KeyCode::Enter));
        assert_eq!(app.site, 11);
    }

    #[test]
    fn test_no_partitions() {
        let files = [PathBuf::from(FILE)];
        let mut app = app(&files);
        assert!(app.partitions.is_none());
        app.handle_key(key(KeyCode::Char('}')));
        assert_eq!(app.message.as_deref(), Some("No partitions found"));
    }

    #[test]
    fn test_quit() {
        let files = [PathBuf::from(FILE)];
        let mut app = app(&files);
        app.handle_key(key(KeyCode::F(1)));
        app.handle_key(key(KeyCode::Esc));
        assert!(!app.should_quit);
        app.handle_key(with(KeyCode::Char('q'), KeyModifiers::CONTROL));
        assert!(app.should_quit);
    }
}
