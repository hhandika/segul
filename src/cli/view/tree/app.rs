//! Viewer state and key bindings.
//!
//! Key bindings follow the alignment viewer. Drawing lives in `ui`; this
//! module only changes state, so it can be tested without a terminal.

use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, Receiver};
use std::thread;

use crossterm::event::{KeyCode, KeyEvent, KeyModifiers, MouseButton, MouseEvent, MouseEventKind};
use ratatui::layout::{Position, Rect};
use segul_core::core::tree::view::{Segment, TreeScale, TreeView, count_trees, read_tree_at};
use segul_core::parser::tree::{NodeId, TreeFmt};

use super::layers::Layers;

/// Columns per unit of tree width, as multiples of the screen width.
pub(super) const ZOOM_LEVELS: [usize; 6] = [1, 2, 4, 8, 16, 32];
const JUMP_COLUMNS: usize = 10;
const MIN_NAME_WIDTH: u16 = 6;
const MAX_NAME_WIDTH: u16 = 60;
pub(super) const DEFAULT_NAME_WIDTH: u16 = 24;

/// Popups and prompts over the tree.
#[derive(Debug, PartialEq)]
pub(super) enum Overlay {
    None,
    Find(String),
    Files { query: String, index: usize },
    Help,
}

/// What the sidebar lists.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum SidebarTab {
    Trees,
    Clades,
}

/// Screen areas from the last draw, used to place mouse clicks.
#[derive(Debug, Default, Clone)]
pub(super) struct Areas {
    pub(super) canvas: Rect,
    pub(super) overview: Rect,
    /// Inside of the sidebar border.
    pub(super) sidebar: Rect,
    /// Layer bar entries and the slot each toggles.
    pub(super) layers: Vec<(Rect, usize)>,
}

pub(super) struct App<'a> {
    pub(super) files: &'a [PathBuf],
    pub(super) file_index: usize,
    input_fmt: TreeFmt,
    pub(super) view: TreeView,
    pub(super) layers: Layers,
    /// Zero-based tree in the file.
    pub(super) tree_index: usize,
    /// Trees in the file, once counted.
    pub(super) tree_count: Option<usize>,
    counter: Option<Receiver<usize>>,
    /// Named clades, for the sidebar.
    pub(super) clades: Vec<NodeId>,
    pub(super) cursor: NodeId,
    pub(super) top_row: usize,
    pub(super) left_col: usize,
    pub(super) zoom: usize,
    pub(super) scale: TreeScale,
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
    /// Newick text to send to the terminal clipboard.
    pub(super) clipboard: Option<String>,
    pub(super) should_quit: bool,
    /// Rows and columns of the canvas in the last draw.
    pub(super) canvas_rows: usize,
    pub(super) canvas_cols: usize,
    /// Changes whenever rows change, so cached drawings can be refreshed.
    pub(super) layout_version: u64,
    /// Overview bins with the width and layout version they were made for.
    pub(super) overview: Option<(usize, u64, Vec<f64>)>,
    pub(super) areas: Areas,
}

impl<'a> App<'a> {
    pub(super) fn new(files: &'a [PathBuf], input_fmt: TreeFmt, view: TreeView) -> Self {
        let layers = Layers::for_tree(view.tree(), None);
        let scale = default_scale(&view);
        let mut app = Self {
            files,
            file_index: 0,
            input_fmt,
            layers,
            tree_index: 0,
            tree_count: None,
            counter: None,
            clades: view.named_clades(),
            cursor: view.row_node(0),
            view,
            top_row: 0,
            left_col: 0,
            zoom: ZOOM_LEVELS[0],
            scale,
            show_sidebar: false,
            sidebar_tab: SidebarTab::Trees,
            name_width: DEFAULT_NAME_WIDTH,
            overlay: Overlay::None,
            find_query: String::new(),
            message: None,
            is_chord: false,
            clipboard: None,
            should_quit: false,
            canvas_rows: 1,
            canvas_cols: 1,
            layout_version: 0,
            overview: None,
            areas: Areas::default(),
        };
        app.count_trees();
        app
    }

    pub(super) fn path(&self) -> &Path {
        &self.files[self.file_index]
    }

    /// Whether trees are still being counted.
    pub(super) fn is_counting(&self) -> bool {
        self.counter.is_some()
    }

    /// Picks up the tree count when the background count finishes.
    pub(super) fn poll_count(&mut self) {
        if let Some(rx) = &self.counter {
            match rx.try_recv() {
                Ok(count) => {
                    self.tree_count = Some(count);
                    self.counter = None;
                }
                Err(mpsc::TryRecvError::Disconnected) => self.counter = None,
                Err(mpsc::TryRecvError::Empty) => (),
            }
        }
    }

    // Counts trees on another thread, so large files open at once.
    fn count_trees(&mut self) {
        let (tx, rx) = mpsc::channel();
        let path = self.path().to_path_buf();
        let fmt = self.input_fmt;
        thread::spawn(move || {
            if let Ok(count) = count_trees(&path, fmt) {
                let _ = tx.send(count);
            }
        });
        self.tree_count = None;
        self.counter = Some(rx);
    }

    /// Canvas columns for branches; tip names use the rest.
    pub(super) fn tree_cols(&self) -> usize {
        self.canvas_cols
            .saturating_sub(self.name_width as usize)
            .max(4)
    }

    /// Width of the whole tree in columns at the current zoom.
    pub(super) fn virtual_cols(&self) -> usize {
        self.tree_cols() * self.zoom
    }

    /// Screen column of position `x`, before scrolling. Column 0 holds the
    /// root stub.
    pub(super) fn column(&self, x: f64) -> isize {
        let (lo, hi) = self.view.x_range(self.scale);
        let span = hi - lo;
        let usable = self.virtual_cols().saturating_sub(2) as f64;
        let fraction = if span > 0.0 { (x - lo) / span } else { 1.0 };
        1 + (fraction * usable).round() as isize - self.left_col as isize
    }

    pub(super) fn cursor_row(&self) -> usize {
        self.view.row_of(self.cursor).unwrap_or(0)
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

    /// Entries in the sidebar list.
    pub(super) fn sidebar_count(&self) -> usize {
        match self.sidebar_tab {
            SidebarTab::Trees => self.tree_count.unwrap_or(self.tree_index + 1),
            SidebarTab::Clades => self.clades.len(),
        }
    }

    /// First sidebar entry shown, keeping the current one in view.
    pub(super) fn sidebar_first(&self, shown: usize) -> usize {
        let current = match self.sidebar_tab {
            SidebarTab::Trees => self.tree_index,
            SidebarTab::Clades => self
                .clades
                .iter()
                .position(|&c| c == self.cursor)
                .unwrap_or(0),
        };
        current
            .saturating_sub(shown / 2)
            .min(self.sidebar_count().saturating_sub(shown))
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
                KeyCode::Char('c' | 'C') => self.cycle_color(),
                KeyCode::Char('s' | 'S') => self.toggle_scale(),
                KeyCode::Char('t' | 'T') => self.toggle_sidebar_tab(),
                _ => self.message = Some(String::from("Unknown chord")),
            }
            return;
        }
        let page = self.canvas_rows.saturating_sub(1).max(1) as isize;
        match key.code {
            KeyCode::Esc => self.should_quit = true,
            KeyCode::Char('q') if ctrl => self.should_quit = true,
            KeyCode::Char('k') if ctrl => {
                self.is_chord = true;
                self.message = Some(String::from(
                    "Ctrl+K: C color branches · S scale · T sidebar list",
                ));
            }
            KeyCode::Char('f') if ctrl => self.overlay = Overlay::Find(String::new()),
            KeyCode::Char('p') if ctrl => {
                self.overlay = Overlay::Files {
                    query: String::new(),
                    index: self.file_index,
                }
            }
            KeyCode::Char('b') if ctrl => self.show_sidebar = !self.show_sidebar,
            KeyCode::Char('c') if ctrl => self.copy_subtree(),
            KeyCode::Char('0') if ctrl => self.set_zoom(ZOOM_LEVELS[0]),
            KeyCode::Char('=' | '+') => self.zoom_by(1),
            KeyCode::Char('-' | '_') => self.zoom_by(-1),
            KeyCode::Char(' ') | KeyCode::Enter => self.toggle_fold(),
            KeyCode::Char('*') => self.unfold_all(),
            KeyCode::Char('[') => self.step_tree(-1),
            KeyCode::Char(']') => self.step_tree(1),
            KeyCode::Char('<' | ',') => self.resize_names(-2),
            KeyCode::Char('>' | '.') => self.resize_names(2),
            KeyCode::Char(c @ '1'..='9') => self.toggle_layer(c as usize - '1' as usize),
            KeyCode::Char('?') | KeyCode::F(1) => self.overlay = Overlay::Help,
            KeyCode::F(2) => self.cycle_color(),
            KeyCode::F(3) if shift => self.find_next(false),
            KeyCode::F(3) => self.find_next(true),
            KeyCode::F(4) => self.toggle_scale(),
            KeyCode::Left if ctrl && shift => self.resize_names(-2),
            KeyCode::Right if ctrl && shift => self.resize_names(2),
            // Ctrl or Option with arrows, or Shift, scroll sideways.
            KeyCode::Left if ctrl || alt || shift => self.scroll_cols(-(JUMP_COLUMNS as isize)),
            KeyCode::Right if ctrl || alt || shift => self.scroll_cols(JUMP_COLUMNS as isize),
            KeyCode::Left => self.go_parent(),
            KeyCode::Right => self.go_child(),
            // Option+↑ and Option+↓, for Mac keyboards without page keys.
            KeyCode::Up if alt => self.move_row(-page),
            KeyCode::Down if alt => self.move_row(page),
            KeyCode::Up => self.move_row(-1),
            KeyCode::Down => self.move_row(1),
            KeyCode::Home => self.set_cursor(self.view.row_node(0)),
            KeyCode::End => self.set_cursor(self.view.row_node(self.view.row_count() - 1)),
            KeyCode::PageUp if ctrl => self.step_file(-1),
            KeyCode::PageDown if ctrl => self.step_file(1),
            KeyCode::PageUp => self.move_row(-page),
            KeyCode::PageDown => self.move_row(page),
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
                // Search as you type, from the first match.
                if let Some(&tip) = self.view.find_tips(&self.find_query).first() {
                    self.go_to_node(tip);
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
            MouseEventKind::ScrollUp if shift => self.scroll_cols(-(JUMP_COLUMNS as isize)),
            MouseEventKind::ScrollDown if shift => self.scroll_cols(JUMP_COLUMNS as isize),
            MouseEventKind::ScrollLeft => self.scroll_cols(-(JUMP_COLUMNS as isize)),
            MouseEventKind::ScrollRight => self.scroll_cols(JUMP_COLUMNS as isize),
            MouseEventKind::ScrollUp => self.scroll_rows(-3),
            MouseEventKind::ScrollDown => self.scroll_rows(3),
            MouseEventKind::Down(MouseButton::Left) => {
                if self.areas.overview.contains(at) {
                    self.jump_overview(mouse.column);
                } else if self.areas.sidebar.contains(at) {
                    self.click_sidebar(at);
                } else if self.areas.canvas.contains(at) {
                    self.click_canvas(at);
                } else if let Some(&(_, slot)) =
                    self.areas.layers.iter().find(|(area, _)| area.contains(at))
                {
                    self.toggle_layer(slot);
                }
            }
            MouseEventKind::Drag(MouseButton::Left) if self.areas.overview.contains(at) => {
                self.jump_overview(mouse.column);
            }
            _ => (),
        }
    }

    pub(super) fn move_row(&mut self, delta: isize) {
        let last = self.view.row_count().saturating_sub(1);
        let row = offset(self.cursor_row(), delta, last);
        self.set_cursor(self.view.row_node(row));
    }

    fn go_parent(&mut self) {
        match self.view.tree().parent(self.cursor) {
            Some(parent) => self.set_cursor(parent),
            None => self.message = Some(String::from("At the root")),
        }
    }

    // Opens a folded clade, or moves to the first child.
    fn go_child(&mut self) {
        if self.view.is_folded(self.cursor) {
            self.toggle_fold();
        } else if let Some(child) = self.view.tree().first_child(self.cursor) {
            self.set_cursor(child);
        }
    }

    pub(super) fn set_cursor(&mut self, id: NodeId) {
        self.cursor = id;
        self.scroll_to_cursor();
    }

    /// Unfolds clades hiding `id`, then moves to it and centers it.
    pub(super) fn go_to_node(&mut self, id: NodeId) {
        if self.view.row_of(id).is_none() {
            self.view.reveal(id);
            self.layout_version += 1;
        }
        self.cursor = id;
        let row = self.cursor_row();
        self.top_row = row.saturating_sub(self.canvas_rows / 2);
        self.clamp_scroll();
        self.scroll_to_cursor();
    }

    /// Folds or unfolds the clade at the cursor. On a tip, folds its parent.
    pub(super) fn toggle_fold(&mut self) {
        let tree = self.view.tree();
        let target = if tree.is_leaf(self.cursor) {
            match tree.parent(self.cursor) {
                Some(parent) => parent,
                None => return,
            }
        } else {
            self.cursor
        };
        self.view.toggle_fold(target);
        self.layout_version += 1;
        self.set_cursor(target);
        self.clamp_scroll();
    }

    fn unfold_all(&mut self) {
        self.view.unfold_all();
        self.layout_version += 1;
        self.clamp_scroll();
        self.scroll_to_cursor();
    }

    fn scroll_to_cursor(&mut self) {
        let row = self.cursor_row();
        if row < self.top_row {
            self.top_row = row;
        } else if row >= self.top_row + self.canvas_rows {
            self.top_row = row + 1 - self.canvas_rows;
        }
        let col = self.column(self.view.x(self.cursor, self.scale));
        let visible = self.tree_cols() as isize;
        if col < 0 {
            self.left_col = self.left_col.saturating_sub((-col) as usize + 1);
        } else if col >= visible {
            self.left_col += (col - visible) as usize + 1;
        }
        self.clamp_scroll();
    }

    pub(super) fn clamp_scroll(&mut self) {
        let rows = self.view.row_count();
        self.top_row = self.top_row.min(rows.saturating_sub(self.canvas_rows));
        let max_left = self.virtual_cols().saturating_sub(self.tree_cols());
        self.left_col = self.left_col.min(max_left);
    }

    fn scroll_rows(&mut self, delta: isize) {
        self.top_row = self.top_row.saturating_add_signed(delta);
        self.clamp_scroll();
    }

    fn scroll_cols(&mut self, delta: isize) {
        self.left_col = self.left_col.saturating_add_signed(delta);
        self.clamp_scroll();
    }

    pub(super) fn zoom_by(&mut self, step: isize) {
        let at = ZOOM_LEVELS
            .iter()
            .position(|&z| z == self.zoom)
            .unwrap_or(0);
        let next = offset(at, step, ZOOM_LEVELS.len() - 1);
        self.set_zoom(ZOOM_LEVELS[next]);
    }

    // Keeps the cursor in view while zooming.
    fn set_zoom(&mut self, zoom: usize) {
        self.zoom = zoom;
        self.left_col = 0;
        self.clamp_scroll();
        self.scroll_to_cursor();
    }

    fn toggle_scale(&mut self) {
        if !self.view.has_lengths() {
            self.message = Some(String::from("The tree has no branch lengths"));
            return;
        }
        self.scale = match self.scale {
            TreeScale::Length => TreeScale::Cladogram,
            TreeScale::Cladogram => TreeScale::Length,
        };
        self.layout_version += 1;
        self.scroll_to_cursor();
    }

    fn toggle_layer(&mut self, slot: usize) {
        if !self.layers.toggle(slot) {
            self.message = Some(format!("No layer {}", slot + 1));
        }
    }

    fn cycle_color(&mut self) {
        self.layers.cycle_color(self.view.tree());
        self.message = Some(match self.layers.color_name() {
            Some(name) => format!("Branches colored by {name}"),
            None => String::from("Branch colors off"),
        });
    }

    fn toggle_sidebar_tab(&mut self) {
        self.show_sidebar = true;
        self.sidebar_tab = match self.sidebar_tab {
            SidebarTab::Trees => SidebarTab::Clades,
            SidebarTab::Clades => SidebarTab::Trees,
        };
    }

    fn resize_names(&mut self, delta: i16) {
        self.name_width = self
            .name_width
            .saturating_add_signed(delta)
            .clamp(MIN_NAME_WIDTH, MAX_NAME_WIDTH);
        self.clamp_scroll();
    }

    fn find_next(&mut self, forward: bool) {
        let matches = self.view.find_tips(&self.find_query);
        if matches.is_empty() {
            self.message = Some(String::from("No results"));
            return;
        }
        // Matches are in tip order, so their position marks the current one.
        let at = matches.iter().position(|&m| m == self.cursor);
        let count = matches.len();
        let next = match (at, forward) {
            (Some(i), true) => (i + 1) % count,
            (Some(i), false) => (i + count - 1) % count,
            (None, true) => 0,
            (None, false) => count - 1,
        };
        self.go_to_node(matches[next]);
    }

    fn copy_subtree(&mut self) {
        let newick = self.view.subtree_newick(self.cursor);
        let tips = self.view.tips_below(self.cursor);
        self.clipboard = Some(newick);
        self.message = Some(format!("Copied subtree of {tips} tip(s) as Newick"));
    }

    // The node on the clicked row nearest the click. Clicking the cursor
    // again folds or unfolds it.
    fn click_canvas(&mut self, at: Position) {
        let canvas = self.areas.canvas;
        let row = self.top_row + (at.y - canvas.y) as usize;
        if row >= self.view.row_count() {
            return;
        }
        let col = (at.x - canvas.x) as isize;
        let nearest = self
            .view
            .segments(row, 1, self.scale)
            .into_iter()
            .filter_map(|s| match s {
                Segment::Branch { node, to, .. } => Some((node, self.column(to))),
                Segment::Connector { .. } => None,
            })
            .min_by_key(|&(node, node_col)| {
                // Labels follow tips, so clicks on them pick the tip.
                let distance = if self.view.tree().is_leaf(node) && col > node_col {
                    0
                } else {
                    (col - node_col).abs()
                };
                (distance, -node_col)
            })
            .map(|(node, _)| node);
        if let Some(node) = nearest {
            if node == self.cursor && !self.view.tree().is_leaf(node) {
                self.toggle_fold();
            } else {
                self.set_cursor(node);
            }
        }
    }

    fn click_sidebar(&mut self, at: Position) {
        let area = self.areas.sidebar;
        let shown = area.height as usize;
        let index = self.sidebar_first(shown) + (at.y - area.y) as usize;
        if index >= self.sidebar_count() {
            return;
        }
        match self.sidebar_tab {
            SidebarTab::Trees => self.open_tree(index),
            SidebarTab::Clades => self.go_to_node(self.clades[index]),
        }
    }

    fn jump_overview(&mut self, column: u16) {
        let area = self.areas.overview;
        if area.width == 0 {
            return;
        }
        let x = column.saturating_sub(area.x) as usize;
        let rows = self.view.row_count();
        let row = (x * rows / area.width as usize).min(rows.saturating_sub(1));
        self.cursor = self.view.row_node(row);
        self.top_row = row.saturating_sub(self.canvas_rows / 2);
        self.clamp_scroll();
    }

    fn step_tree(&mut self, step: isize) {
        let index = self.tree_index.saturating_add_signed(step);
        if step < 0 && self.tree_index == 0 {
            self.message = Some(String::from("At the first tree"));
            return;
        }
        if self.tree_count.is_some_and(|n| index >= n) {
            self.message = Some(String::from("At the last tree"));
            return;
        }
        self.open_tree(index);
    }

    /// Opens tree `index` of the current file, keeping the view settings.
    pub(super) fn open_tree(&mut self, index: usize) {
        if index == self.tree_index {
            return;
        }
        match read_tree_at(self.path(), self.input_fmt, index) {
            Ok(tree) => {
                self.tree_index = index;
                self.show(TreeView::new(tree));
            }
            Err(_) => {
                self.message = Some(String::from("At the last tree"));
                if self.tree_count.is_none() {
                    self.tree_count = Some(self.tree_index + 1);
                }
            }
        }
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

    /// Opens the first tree of another input file, keeping the view settings.
    pub(super) fn open_file(&mut self, index: usize) {
        if index == self.file_index {
            return;
        }
        match read_tree_at(&self.files[index], self.input_fmt, 0) {
            Ok(tree) => {
                self.file_index = index;
                self.tree_index = 0;
                self.show(TreeView::new(tree));
                self.count_trees();
            }
            Err(e) => self.message = Some(format!("{e:#}")),
        }
    }

    fn show(&mut self, view: TreeView) {
        self.layers = Layers::for_tree(view.tree(), Some(&self.layers));
        if !view.has_lengths() {
            self.scale = TreeScale::Cladogram;
        }
        self.clades = view.named_clades();
        self.cursor = view.row_node(0);
        self.view = view;
        self.top_row = 0;
        self.left_col = 0;
        self.layout_version += 1;
    }
}

fn default_scale(view: &TreeView) -> TreeScale {
    if view.has_lengths() {
        TreeScale::Length
    } else {
        TreeScale::Cladogram
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

    use segul_core::parser::tree::Tree;

    const SIMPLE: &str = "tests/files/trees/simple.nwk";
    const MULTI: &str = "tests/files/trees/multi.nwk";

    fn key(code: KeyCode) -> KeyEvent {
        KeyEvent::new(code, KeyModifiers::NONE)
    }

    fn with(code: KeyCode, modifiers: KeyModifiers) -> KeyEvent {
        KeyEvent::new(code, modifiers)
    }

    fn open(path: &Path) -> Tree {
        read_tree_at(path, TreeFmt::Auto, 0).expect("Failed reading fixture")
    }

    fn app(files: &[PathBuf]) -> App<'_> {
        let mut app = App::new(files, TreeFmt::Auto, TreeView::new(open(&files[0])));
        app.canvas_rows = 3;
        app.canvas_cols = 60;
        app
    }

    fn cursor_label(app: &App) -> String {
        app.view.label(app.cursor).unwrap_or("").to_string()
    }

    #[test]
    fn test_cursor_moves() {
        let files = [PathBuf::from(SIMPLE)];
        let mut app = app(&files);
        assert_eq!(cursor_label(&app), "Taxon_A");
        app.handle_key(key(KeyCode::Down));
        assert_eq!(cursor_label(&app), "Taxon_B");
        app.handle_key(key(KeyCode::End));
        assert_eq!(cursor_label(&app), "Taxon_E");
        // The canvas holds three rows, so the last row scrolls into view.
        assert_eq!(app.top_row, 2);
        app.handle_key(key(KeyCode::Down));
        assert_eq!(cursor_label(&app), "Taxon_E");
        app.handle_key(key(KeyCode::Left));
        assert_eq!(cursor_label(&app), "85");
        app.handle_key(key(KeyCode::Right));
        assert_eq!(cursor_label(&app), "Taxon_C");
        app.handle_key(key(KeyCode::Home));
        assert_eq!(app.top_row, 0);
    }

    #[test]
    fn test_fold() {
        let files = [PathBuf::from(SIMPLE)];
        let mut app = app(&files);
        app.handle_key(key(KeyCode::Char(' ')));
        assert_eq!(cursor_label(&app), "98");
        assert!(app.view.is_folded(app.cursor));
        assert_eq!(app.view.row_count(), 4);
        // Right opens the folded clade again.
        app.handle_key(key(KeyCode::Right));
        assert_eq!(app.view.row_count(), 5);
        app.handle_key(key(KeyCode::Enter));
        app.handle_key(key(KeyCode::Char('*')));
        assert_eq!(app.view.row_count(), 5);
    }

    #[test]
    fn test_find_reveals() {
        let files = [PathBuf::from(SIMPLE)];
        let mut app = app(&files);
        app.handle_key(key(KeyCode::Char(' ')));
        app.handle_key(with(KeyCode::Char('f'), KeyModifiers::CONTROL));
        app.handle_key(key(KeyCode::Char('b')));
        assert_eq!(cursor_label(&app), "Taxon_B");
        assert_eq!(app.view.row_count(), 5);
        app.handle_key(key(KeyCode::Esc));
        app.find_query = String::from("taxon");
        app.handle_key(key(KeyCode::F(3)));
        assert_eq!(cursor_label(&app), "Taxon_C");
        app.handle_key(with(KeyCode::F(3), KeyModifiers::SHIFT));
        assert_eq!(cursor_label(&app), "Taxon_B");
    }

    #[test]
    fn test_trees_in_file() {
        let files = [PathBuf::from(MULTI)];
        let mut app = app(&files);
        app.handle_key(key(KeyCode::Char('[')));
        assert_eq!(app.message.as_deref(), Some("At the first tree"));
        app.handle_key(key(KeyCode::Char(']')));
        app.handle_key(key(KeyCode::Char(']')));
        assert_eq!(app.tree_index, 2);
        assert_eq!(cursor_label(&app), "TaxonE");
        app.handle_key(key(KeyCode::Char(']')));
        assert_eq!(app.tree_index, 2);
        assert_eq!(app.message.as_deref(), Some("At the last tree"));
    }

    #[test]
    fn test_count_trees() {
        let files = [PathBuf::from(MULTI)];
        let mut app = app(&files);
        while app.is_counting() {
            app.poll_count();
            thread::yield_now();
        }
        assert_eq!(app.tree_count, Some(3));
        assert_eq!(app.sidebar_count(), 3);
    }

    #[test]
    fn test_switch_files() {
        let files = [PathBuf::from(SIMPLE), PathBuf::from(MULTI)];
        let mut app = app(&files);
        app.handle_key(key(KeyCode::Char('2')));
        app.handle_key(with(KeyCode::PageDown, KeyModifiers::CONTROL));
        assert_eq!(app.file_index, 1);
        assert_eq!(cursor_label(&app), "TaxonA");
        // Layer choices carry over.
        assert!(app.layers.slots[1].visible);
        app.handle_key(with(KeyCode::PageUp, KeyModifiers::CONTROL));
        assert_eq!(app.file_index, 0);
    }

    #[test]
    fn test_zoom_and_scale() {
        let files = [PathBuf::from(SIMPLE)];
        let mut app = app(&files);
        app.handle_key(key(KeyCode::Char('+')));
        assert_eq!(app.zoom, 2);
        assert_eq!(app.virtual_cols(), 72);
        app.handle_key(with(KeyCode::Char('0'), KeyModifiers::CONTROL));
        assert_eq!(app.zoom, 1);
        app.handle_key(key(KeyCode::F(4)));
        assert_eq!(app.scale, TreeScale::Cladogram);
        // Tips line up on the right in a cladogram.
        assert_eq!(app.column(app.view.x(app.cursor, app.scale)), 35);
    }

    #[test]
    fn test_layers_and_copy() {
        let files = [PathBuf::from(SIMPLE)];
        let mut app = app(&files);
        app.handle_key(key(KeyCode::Char('9')));
        assert_eq!(app.message.as_deref(), Some("No layer 9"));
        app.handle_key(key(KeyCode::F(2)));
        assert_eq!(app.message.as_deref(), Some("Branches colored by Support"));
        app.handle_key(key(KeyCode::Left));
        app.handle_key(with(KeyCode::Char('c'), KeyModifiers::CONTROL));
        assert_eq!(
            app.clipboard.as_deref(),
            Some("(Taxon_A:0.12,Taxon_B:0.34)98;")
        );
    }

    #[test]
    fn test_quit() {
        let files = [PathBuf::from(SIMPLE)];
        let mut app = app(&files);
        app.handle_key(key(KeyCode::F(1)));
        assert_eq!(app.overlay, Overlay::Help);
        app.handle_key(key(KeyCode::Esc));
        assert!(!app.should_quit);
        app.handle_key(key(KeyCode::Esc));
        assert!(app.should_quit);
    }
}
