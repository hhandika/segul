//! Draws the viewer, laid out like the alignment viewer: toolbar, overview,
//! layer bar, sidebar and tree canvas, then the status bar.
//!
//! The canvas is drawn in passes: branches, then layer labels, then tip
//! names. Only the rows in view are drawn.

use ratatui::Frame;
use ratatui::buffer::Buffer;
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::style::{Color, Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, BorderType, Clear, Paragraph, Wrap};
use segul_core::color::residue::{ResiduePalette, ResidueScheme};
use segul_core::color::{Colormap, Theme};
use segul_core::core::tree::view::{Segment, TreeScale};
use segul_core::helper::types::DataType;
use segul_core::parser::tree::Rootedness;

use crate::cli::view::align::palette::{is_truecolor, to_color};
use crate::cli::view::{centered, count, truncate};

use super::app::{App, Overlay, SidebarTab};

const MIN_WIDTH: u16 = 40;
const MIN_HEIGHT: u16 = 10;
const SIDEBAR_WIDTH: u16 = 24;
// Columns between ruler labels, at least.
const RULER_GAP: usize = 10;
// Bar heights for the overview, in eighths of a cell.
const BARS: [char; 9] = [' ', '▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];
// Cells in the branch color swatch on the layer bar.
const SWATCH: usize = 5;

// Line directions in a canvas cell.
const UP: u8 = 1;
const DOWN: u8 = 2;
const LEFT: u8 = 4;
const RIGHT: u8 = 8;

/// Colors for the tree viewer.
pub(super) struct TreePalette {
    accent: Color,
    tertiary: Color,
    colormap: Colormap,
    is_truecolor: bool,
}

impl TreePalette {
    /// The shared accent for `theme`. Terminals without 24-bit color get
    /// ANSI cyan on dark themes, which every terminal theme keeps readable.
    pub(super) fn new(theme: Theme, colormap: Colormap) -> Self {
        let is_truecolor = is_truecolor();
        let source = ResiduePalette::new(&ResidueScheme::Clustal, &theme, &DataType::Dna);
        let accent = match (theme, is_truecolor) {
            (Theme::Dark, false) => Color::Cyan,
            _ => to_color(source.accent, is_truecolor),
        };
        Self {
            accent,
            tertiary: to_color(source.tertiary, is_truecolor),
            colormap,
            is_truecolor,
        }
    }

    fn colormap_at(&self, t: f64) -> Color {
        to_color(self.colormap.at(t.clamp(0.0, 1.0)), self.is_truecolor)
    }
}

pub(super) fn render(frame: &mut Frame, app: &mut App, palette: &TreePalette) {
    let area = frame.area();
    if area.width < MIN_WIDTH || area.height < MIN_HEIGHT {
        let note = format!("Enlarge the terminal to at least {MIN_WIDTH}x{MIN_HEIGHT}");
        frame.render_widget(Paragraph::new(note).wrap(Wrap { trim: true }), area);
        return;
    }
    let show_overview = area.height >= 16;
    let [toolbar, overview, layer_bar, main, status] = Layout::vertical([
        Constraint::Length(1),
        Constraint::Length(if show_overview { 3 } else { 0 }),
        Constraint::Length(1),
        Constraint::Min(3),
        Constraint::Length(1),
    ])
    .areas(area);
    let [sidebar, tree_area] = Layout::horizontal([
        Constraint::Length(if app.show_sidebar { SIDEBAR_WIDTH } else { 0 }),
        Constraint::Min(1),
    ])
    .areas(main);
    let [ruler_top, canvas, ruler_bottom] = Layout::vertical([
        Constraint::Length(1),
        Constraint::Min(1),
        Constraint::Length(1),
    ])
    .areas(tree_area);

    app.canvas_rows = canvas.height as usize;
    app.canvas_cols = canvas.width as usize;
    app.clamp_scroll();
    app.areas.canvas = canvas;
    app.areas.overview = overview;
    app.areas.sidebar = if app.show_sidebar {
        sidebar_block().inner(sidebar)
    } else {
        Rect::default()
    };

    render_toolbar(frame, toolbar, app);
    if show_overview {
        render_overview(frame.buffer_mut(), overview, app, palette);
    }
    render_layer_bar(frame.buffer_mut(), layer_bar, app, palette);
    if app.show_sidebar {
        render_sidebar(frame, sidebar, app, palette);
    }
    render_ruler(frame.buffer_mut(), ruler_top, app, Edge::Top);
    render_ruler(frame.buffer_mut(), ruler_bottom, app, Edge::Bottom);
    render_canvas(frame.buffer_mut(), canvas, app, palette);
    render_status(frame, status, app, palette);

    match &app.overlay {
        Overlay::Help => render_help(frame, area),
        Overlay::Files { query, index } => render_files(frame, area, app, query, *index, palette),
        _ => (),
    }
}

fn render_toolbar(frame: &mut Frame, area: Rect, app: &App) {
    let view = &app.view;
    let name = app
        .path()
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let mut summary = format!("{} tips", count(view.tip_count()));
    match view.rootedness() {
        Rootedness::Rooted => summary.push_str(" · rooted"),
        Rootedness::Unrooted => summary.push_str(" · unrooted"),
        Rootedness::Unknown => (),
    }
    if view.is_rerooted() {
        summary.push_str(" · rerooted");
    }
    if !view.has_lengths() {
        summary.push_str(" · no lengths");
    }
    if view.tree().is_network() {
        summary.push_str(" · network");
    }
    let mut spans = vec![Span::from(name).bold()];
    if app.files.len() > 1 {
        spans.push(
            Span::from(format!(
                " ({}/{})",
                count(app.file_index + 1),
                count(app.files.len())
            ))
            .dim(),
        );
    }
    let total = app.tree_count.map_or_else(|| String::from("?"), count);
    let mut tree = format!("Tree {} of {}", count(app.tree_index + 1), total);
    if let Some(tree_name) = &view.tree().metadata.name {
        tree.push_str(&format!(" · {tree_name}"));
    }
    spans.extend([
        Span::from("  "),
        Span::from(summary),
        Span::from("  │ ").dim(),
        Span::from(match app.scale {
            TreeScale::Length => "Lengths",
            TreeScale::Cladogram => "Cladogram",
        }),
        Span::from(format!(" · Zoom {}×", app.zoom)),
        Span::from("  │ ").dim(),
        Span::from(tree),
    ]);
    frame.render_widget(Line::from(spans), area);
}

// Mean tip distance from the root for runs of rows, as bars two cells tall,
// then a line marking the rows in view.
fn render_overview(buf: &mut Buffer, area: Rect, app: &mut App, palette: &TreePalette) {
    let width = area.width as usize;
    let is_fresh =
        matches!(&app.overview, Some((w, v, _)) if *w == width && *v == app.layout_version);
    if !is_fresh {
        let bins = app.view.overview(width, app.scale);
        app.overview = Some((width, app.layout_version, bins));
    }
    let Some((_, _, bins)) = &app.overview else {
        return;
    };
    for (x, &value) in bins.iter().enumerate() {
        // At least one eighth, so every run of rows shows.
        let level = ((value * 16.0).round() as usize).clamp(1, 16);
        let top = if level > 8 { BARS[level - 8] } else { ' ' };
        let bottom = BARS[level.min(8)];
        let col = area.x + x as u16;
        buf[(col, area.y)].set_char(top).set_fg(palette.tertiary);
        buf[(col, area.y + 1)]
            .set_char(bottom)
            .set_fg(palette.tertiary);
    }
    let rows = app.view.row_count().max(1);
    let from = app.top_row * width / rows;
    let to = ((app.top_row + app.canvas_rows).min(rows) * width).div_ceil(rows);
    let y = area.y + 2;
    for x in 0..width {
        let cell = &mut buf[(area.x + x as u16, y)];
        if x >= from && x < to.max(from + 1) {
            cell.set_char('━').set_fg(palette.accent);
        } else {
            cell.set_char('─').set_style(Style::new().dim());
        }
    }
}

// Each layer with its number and whether it is shown, then the layer
// coloring the branches. Remembers where each entry is for clicks.
fn render_layer_bar(buf: &mut Buffer, area: Rect, app: &mut App, palette: &TreePalette) {
    app.areas.layers.clear();
    let right = area.x + area.width;
    let mut x = area.x;
    let put = |buf: &mut Buffer, x: &mut u16, text: &str, style: Style| {
        let room = right.saturating_sub(*x) as usize;
        let (end, _) = buf.set_stringn(*x, area.y, text, room, style);
        *x = end;
    };
    put(buf, &mut x, "Layers ", Style::new().dim());
    for (i, slot) in app.layers.slots.iter().enumerate() {
        let start = x;
        let mark = if slot.visible { '■' } else { '□' };
        let number = if i < 9 {
            format!(" {}", i + 1)
        } else {
            String::from(" ")
        };
        put(buf, &mut x, &number, Style::new().dim());
        let style = if slot.visible {
            Style::new().fg(palette.accent)
        } else {
            Style::new().dim()
        };
        put(buf, &mut x, &format!("{mark}{} ", slot.layer.name()), style);
        if x > start {
            let hit = Rect::new(start, area.y, x - start, 1);
            app.areas.layers.push((hit, i));
        }
    }
    if let Some(name) = app.layers.color_name() {
        put(
            buf,
            &mut x,
            &format!(" │ Color {name} "),
            Style::new().dim(),
        );
        for i in 0..SWATCH {
            let color = palette.colormap_at(i as f64 / (SWATCH - 1) as f64);
            put(buf, &mut x, "█", Style::new().fg(color));
        }
    }
}

fn sidebar_block() -> Block<'static> {
    Block::bordered().border_type(BorderType::Rounded)
}

fn render_sidebar(frame: &mut Frame, area: Rect, app: &App, palette: &TreePalette) {
    let is_clades = app.sidebar_tab == SidebarTab::Clades;
    let tab = |label: &'static str, active: bool| {
        if active {
            Span::from(label).fg(palette.accent).bold()
        } else {
            Span::from(label).dim()
        }
    };
    let title = Line::from(vec![
        Span::from(" "),
        tab("Trees", !is_clades),
        Span::from(" · ").dim(),
        tab("Clades", is_clades),
        Span::from(" "),
    ]);
    let block = sidebar_block().title(title);
    let inner = block.inner(area);
    frame.render_widget(block, area);
    let shown = inner.height as usize;
    let first = app.sidebar_first(shown);
    let width = inner.width as usize;
    let entries: Vec<(String, bool)> = if is_clades {
        app.clades
            .iter()
            .skip(first)
            .take(shown)
            .map(|&id| {
                let label = app.view.label(id).unwrap_or("");
                let tips = count(app.view.tips_below(id));
                (format!("{label} · {tips}"), id == app.cursor)
            })
            .collect()
    } else {
        (first..(first + shown).min(app.sidebar_count()))
            .map(|i| {
                let mut label = format!("Tree {}", count(i + 1));
                if i == app.tree_index
                    && let Some(name) = &app.view.tree().metadata.name
                {
                    label.push_str(&format!(" · {name}"));
                }
                (label, i == app.tree_index)
            })
            .collect()
    };
    let mut lines: Vec<Line> = entries
        .into_iter()
        .map(|(label, is_current)| {
            let marker = if is_current { "▌" } else { " " };
            let text = truncate(&format!("{marker}{label}"), width);
            if is_current {
                Line::from(text).fg(palette.accent).bold()
            } else {
                Line::from(text)
            }
        })
        .collect();
    if is_clades && app.clades.is_empty() {
        lines.push(Line::from(" No named clades").dim());
    } else if !is_clades && app.is_counting() && lines.len() < shown {
        lines.push(Line::from(" Counting trees…").dim());
    }
    frame.render_widget(Paragraph::new(lines), inner);
}

/// Which side of the canvas a ruler sits on.
#[derive(Clone, Copy, PartialEq)]
enum Edge {
    Top,
    Bottom,
}

// Branch length scale with round steps, ticks pointing at the canvas. For
// cladograms, a note on the top ruler only.
fn render_ruler(buf: &mut Buffer, area: Rect, app: &App, edge: Edge) {
    let width = area.width as usize;
    if app.scale == TreeScale::Cladogram {
        if edge == Edge::Top {
            buf.set_stringn(area.x, area.y, "Cladogram", width, Style::new().dim());
        }
        return;
    }
    let tick_mark = match edge {
        Edge::Top => '╵',
        Edge::Bottom => '╷',
    };
    let (lo, hi) = app.view.x_range(app.scale);
    let usable = app.virtual_cols().saturating_sub(2).max(1) as f64;
    let per_col = (hi - lo) / usable;
    if per_col <= 0.0 || !per_col.is_finite() {
        return;
    }
    let step = nice_step(per_col * RULER_GAP as f64);
    let decimals = (-step.log10().floor()).max(0.0) as usize;
    let mut tick = (lo / step).ceil() * step;
    let mut free_from = 0isize;
    while tick <= hi + step * 1e-9 {
        let col = app.column(tick);
        if col >= width as isize {
            break;
        }
        if col >= free_from && col >= 0 {
            let label = format!("{tick_mark}{:.*}", decimals, tick);
            let room = width - col as usize;
            buf.set_stringn(
                area.x + col as u16,
                area.y,
                &label,
                room,
                Style::new().dim(),
            );
            free_from = col + label.chars().count() as isize + 1;
        }
        tick += step;
    }
}

// Smallest of 1, 2, or 5 times a power of ten at or above `value`.
fn nice_step(value: f64) -> f64 {
    let power = 10f64.powf(value.log10().floor());
    [1.0, 2.0, 5.0, 10.0]
        .iter()
        .map(|m| m * power)
        .find(|&s| s >= value)
        .unwrap_or(10.0 * power)
}

struct Canvas {
    width: usize,
    height: usize,
    lines: Vec<u8>,
    styles: Vec<Style>,
}

impl Canvas {
    fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            lines: vec![0; width * height],
            styles: vec![Style::new(); width * height],
        }
    }

    fn add(&mut self, col: isize, y: usize, line: u8, style: Style) {
        if col < 0 || col as usize >= self.width || y >= self.height {
            return;
        }
        let i = y * self.width + col as usize;
        self.lines[i] |= line;
        self.styles[i] = style;
    }

    fn draw(&self, buf: &mut Buffer, area: Rect) {
        for y in 0..self.height {
            for x in 0..self.width {
                let i = y * self.width + x;
                if let Some(c) = line_char(self.lines[i]) {
                    buf[(area.x + x as u16, area.y + y as u16)]
                        .set_char(c)
                        .set_style(self.styles[i]);
                }
            }
        }
    }
}

fn line_char(lines: u8) -> Option<char> {
    let c = match lines {
        0 => return None,
        l if l == UP | DOWN || l == UP || l == DOWN => '│',
        l if l == DOWN | RIGHT => '┌',
        l if l == UP | RIGHT => '└',
        l if l == DOWN | LEFT => '┐',
        l if l == UP | LEFT => '┘',
        l if l == UP | DOWN | RIGHT => '├',
        l if l == UP | DOWN | LEFT => '┤',
        l if l == LEFT | RIGHT | DOWN => '┬',
        l if l == LEFT | RIGHT | UP => '┴',
        l if l == UP | DOWN | LEFT | RIGHT => '┼',
        _ => '─',
    };
    Some(c)
}

fn render_canvas(buf: &mut Buffer, area: Rect, app: &App, palette: &TreePalette) {
    let view = &app.view;
    let tree = view.tree();
    let first = app.top_row;
    let rows = (area.height as usize).min(view.row_count().saturating_sub(first));
    let width = area.width as usize;
    let mut canvas = Canvas::new(width, rows);
    let segments = view.segments(first, rows, app.scale);
    let style_of = |node| {
        let mut style = match app.layers.color_at(tree, node) {
            Some(t) => Style::new().fg(palette.colormap_at(t)),
            None => Style::new(),
        };
        if node == app.cursor {
            style = style.fg(palette.accent).bold();
        }
        style
    };

    // Branches and connectors.
    for segment in &segments {
        match *segment {
            Segment::Connector {
                node,
                x,
                top,
                bottom,
                ..
            } => {
                let col = app.column(x);
                let style = style_of(node);
                for row in top.max(first)..=bottom.min(first + rows - 1) {
                    let mut line = 0;
                    if row > top {
                        line |= UP;
                    }
                    if row < bottom {
                        line |= DOWN;
                    }
                    canvas.add(col, row - first, line, style);
                }
            }
            Segment::Branch {
                node,
                row,
                from,
                to,
                ..
            } => {
                let end = app.column(to);
                // The root gets a one-column stub.
                let start = if tree.is_root(node) {
                    end - 1
                } else {
                    app.column(from)
                };
                let style = style_of(node);
                for col in start.max(0)..=end.min(width as isize - 1) {
                    let mut line = 0;
                    if col > start {
                        line |= LEFT;
                    }
                    if col < end {
                        line |= RIGHT;
                    }
                    if line != 0 {
                        canvas.add(col, row - first, line, style);
                    }
                }
            }
        }
    }
    canvas.draw(buf, area);

    // Layer values on internal branches, right before the node, when they
    // fit on the branch.
    let label_style = Style::new().fg(palette.tertiary);
    for segment in &segments {
        let Segment::Branch {
            node,
            row,
            from,
            to,
            ..
        } = *segment
        else {
            continue;
        };
        if tree.is_leaf(node) || tree.is_root(node) {
            continue;
        }
        let Some(text) = app.layers.labels(tree, node) else {
            continue;
        };
        let (start, end) = (app.column(from), app.column(to));
        let len = text.chars().count() as isize;
        if end - start > len {
            put_text(buf, area, end - len, row - first, &text, label_style);
        }
    }

    // Tip names and folded clades.
    let matches = match &app.overlay {
        Overlay::Find(_) => view.find_tips(&app.find_query),
        _ => Vec::new(),
    };
    for y in 0..rows {
        let node = view.row_node(first + y);
        let col = app.column(view.x(node, app.scale)) + 1;
        let label = view.label(node).unwrap_or("");
        let mut name_style = Style::new();
        if node == app.cursor {
            name_style = name_style.fg(palette.accent).bold();
        }
        if matches.contains(&node) {
            name_style = name_style.underlined();
        }
        let mut col = if view.is_folded(node) {
            let tips = format!("▸ [+{} tips]", count(view.tips_below(node)));
            let text = if label.is_empty() {
                tips
            } else {
                format!("{tips} {label}")
            };
            put_text(buf, area, col, y, &text, name_style.fg(palette.accent))
        } else {
            put_text(buf, area, col, y, &format!(" {label}"), name_style)
        };
        if tree.is_leaf(node)
            && let Some(values) = app.layers.labels(tree, node)
        {
            col = put_text(buf, area, col, y, "  ", Style::new());
            put_text(buf, area, col, y, &values, label_style);
        }
    }
}

// Writes `text` at canvas column `col` (which may be left of the canvas),
// clipped to the canvas, and returns the column after it.
fn put_text(buf: &mut Buffer, area: Rect, col: isize, y: usize, text: &str, style: Style) -> isize {
    let width = area.width as isize;
    for (x, c) in (col..).zip(text.chars()) {
        if x >= width {
            break;
        }
        if x >= 0 {
            buf[(area.x + x as u16, area.y + y as u16)]
                .set_char(c)
                .set_style(style);
        }
    }
    col + text.chars().count() as isize
}

fn render_status(frame: &mut Frame, area: Rect, app: &App, palette: &TreePalette) {
    let left = match &app.overlay {
        Overlay::Find(input) => {
            let matches = app.view.find_tips(input);
            let at = matches.iter().position(|&m| m == app.cursor);
            let found = match (matches.len(), at) {
                (0, _) if !input.is_empty() => String::from("  No results"),
                (0, _) => String::new(),
                (n, Some(i)) => format!("  {} of {}", count(i + 1), count(n)),
                (n, None) => format!("  {} matches", count(n)),
            };
            Line::from(vec![
                Span::from(" Find tip: ").fg(palette.accent),
                Span::from(format!("{input}▏")),
                Span::from(found).dim(),
                Span::from("  Enter next · Shift+Enter previous · Esc close").dim(),
            ])
        }
        _ => match &app.message {
            Some(message) => Line::from(format!(" {message}")).fg(palette.accent),
            None => Line::from(cursor_status(app)),
        },
    };
    let mut right = Vec::new();
    if app.files.len() > 1 {
        right.push(Span::from("Ctrl+P Files · ").dim());
    }
    right.push(Span::from("? Help ").dim());
    frame.render_widget(left, area);
    frame.render_widget(Line::from(right).right_aligned(), area);
}

fn cursor_status(app: &App) -> String {
    let view = &app.view;
    let tree = view.tree();
    let node = app.cursor;
    let mut parts = Vec::new();
    let label = view.label(node).filter(|l| !l.is_empty());
    if tree.is_leaf(node) {
        parts.push(label.unwrap_or("Unnamed tip").to_string());
    } else {
        let kind = if tree.is_root(node) { "Root" } else { "Clade" };
        // Numeric internal labels are support values, listed below.
        match label.filter(|l| l.trim().parse::<f64>().is_err()) {
            Some(name) => parts.push(format!("{kind} {name}")),
            None => parts.push(kind.to_string()),
        }
        let mut tips = format!("{} tips", count(view.tips_below(node)));
        if view.is_folded(node) {
            tips.push_str(" (folded)");
        }
        parts.push(tips);
    }
    if let Some(info) = tree.hybrid_info(node) {
        parts.push(format!("hybrid {}", info.tag));
    }
    parts.extend(app.layers.describe(tree, node));
    format!(" {}", parts.join(" · "))
}

// Keys for Windows and Linux, keys for macOS, and the action. macOS
// terminals do not pass Command to programs, so the Mac keys use Option and
// fn instead.
const HELP: &[(&str, &str, &str)] = &[
    ("↑ / ↓", "↑ / ↓", "Previous / next row"),
    ("← / →", "← / →", "Parent / first child"),
    (
        "Space or Enter, *",
        "Space or Return, *",
        "Fold clade, unfold all",
    ),
    ("PgUp / PgDn", "Option+↑ / ↓", "Page rows"),
    ("Home / End", "fn+← / →", "First / last row"),
    ("Shift+← / →", "Shift+← / →", "Scroll sideways"),
    (
        "+ / - / 0, Ctrl+wheel",
        "+ / - / 0, Option+wheel",
        "Zoom in / out / fit",
    ),
    (
        "r, R, Ctrl+R",
        "r, R, Ctrl+R",
        "Root on branch, at midpoint, original",
    ),
    ("[ / ]", "[ / ]", "Previous / next tree in file"),
    (
        "Ctrl+F, F3",
        "Ctrl+F, fn+F3",
        "Find tip; Enter or F3 next, Shift back",
    ),
    ("1 to 9", "1 to 9", "Show or hide a layer"),
    ("F2 or Ctrl+K C", "Ctrl+K C", "Color branches by a layer"),
    ("F4 or Ctrl+K S", "Ctrl+K S", "Lengths or cladogram"),
    (
        "Ctrl+B, Ctrl+K T",
        "Ctrl+B, Ctrl+K T",
        "Toggle sidebar / switch its list",
    ),
    ("< / >", "< / >", "Narrow / widen names"),
    ("Ctrl+C", "Ctrl+C", "Copy subtree as Newick"),
    (
        "Ctrl+PgUp / PgDn, Ctrl+P",
        "Ctrl+fn+↑ / ↓, Ctrl+P",
        "Previous / next / pick file",
    ),
    (
        "Mouse",
        "Mouse",
        "Click node · click again folds · wheel scrolls",
    ),
    ("F1 or ?", "?", "Show this help"),
    ("Ctrl+Q or Esc", "Ctrl+Q or Esc", "Quit"),
];

fn render_help(frame: &mut Frame, area: Rect) {
    let width = 78.min(area.width);
    let height = (HELP.len() as u16 + 2).min(area.height);
    let popup = centered(area, width, height);
    let lines: Vec<Line> = HELP
        .iter()
        .map(|(keys, mac_keys, action)| {
            let keys = if cfg!(target_os = "macos") {
                mac_keys
            } else {
                keys
            };
            Line::from(vec![
                Span::from(format!(" {:25}", keys)).bold(),
                Span::from(*action),
            ])
        })
        .collect();
    let block = Block::bordered()
        .border_type(BorderType::Rounded)
        .title(" Keyboard shortcuts ")
        .title_bottom(Line::from(" Esc close ").right_aligned());
    frame.render_widget(Clear, popup);
    frame.render_widget(Paragraph::new(lines).block(block), popup);
}

fn render_files(
    frame: &mut Frame,
    area: Rect,
    app: &App,
    query: &str,
    index: usize,
    palette: &TreePalette,
) {
    let width = 78.min(area.width);
    let height = 16.min(area.height);
    let popup = centered(area, width, height);
    let block = Block::bordered()
        .border_type(BorderType::Rounded)
        .title(format!(" Open file ({}) ", count(app.files.len())))
        .title_bottom(
            Line::from(" Type to filter · ↑↓ choose · Enter open · Esc close ").right_aligned(),
        );
    let inner = block.inner(popup);
    frame.render_widget(Clear, popup);
    frame.render_widget(block, popup);
    let [input, list] = Layout::vertical([Constraint::Length(1), Constraint::Min(1)]).areas(inner);
    frame.render_widget(
        Line::from(vec![
            Span::from("> ").fg(palette.accent),
            Span::from(format!("{query}▏")),
        ]),
        input,
    );
    let matches = app.matching_files(query);
    let shown = list.height as usize;
    let at = matches.iter().position(|&i| i == index).unwrap_or(0);
    let first = at.saturating_sub(shown.saturating_sub(1));
    let lines: Vec<Line> = matches
        .iter()
        .skip(first)
        .take(shown)
        .map(|&i| {
            let path = app.files[i].display().to_string();
            let line = Line::from(format!(" {}", truncate(&path, list.width as usize - 1)));
            if i == index {
                line.fg(palette.accent).bold().reversed()
            } else if i == app.file_index {
                line.bold()
            } else {
                line
            }
        })
        .collect();
    frame.render_widget(Paragraph::new(lines), list);
}

#[cfg(test)]
mod test {
    use super::*;

    use std::path::PathBuf;

    use ratatui::Terminal;
    use ratatui::backend::TestBackend;
    use segul_core::core::tree::view::{TreeView, read_tree_at};
    use segul_core::parser::tree::TreeFmt;

    fn app(files: &[PathBuf]) -> App<'_> {
        let tree = read_tree_at(&files[0], TreeFmt::Auto, 0).expect("Failed reading fixture");
        App::new(files, TreeFmt::Auto, TreeView::new(tree))
    }

    fn palette() -> TreePalette {
        TreePalette::new(Theme::Dark, Colormap::Viridis)
    }

    fn screen(app: &mut App, palette: &TreePalette) -> String {
        let mut terminal =
            Terminal::new(TestBackend::new(80, 24)).expect("Failed creating terminal");
        terminal
            .draw(|frame| render(frame, app, palette))
            .expect("Failed drawing");
        let buffer = terminal.backend().buffer();
        (0..buffer.area.height)
            .map(|y| {
                (0..buffer.area.width)
                    .map(|x| buffer[(x, y)].symbol())
                    .collect::<String>()
            })
            .collect::<Vec<_>>()
            .join("\n")
    }

    #[test]
    fn test_render() {
        let files = [PathBuf::from("tests/files/trees/simple.nwk")];
        let mut app = app(&files);
        let palette = palette();
        let text = screen(&mut app, &palette);
        assert!(text.contains("simple.nwk"));
        assert!(text.contains("5 tips"));
        assert!(text.contains("Layers  1■Support  2□Lengths"));
        assert!(text.contains(" Taxon_A"));
        assert!(text.contains(" Taxon_E"));
        // Support values sit on their branches.
        assert!(text.contains("98┴────────────────────────── Taxon_B"));
        // The root sits midway between its clades.
        assert!(text.contains("─┤                              ┌─────── Taxon_C"));
        assert!(text.contains("85┼"));
        assert!(text.contains("╵0.0"));
        assert!(text.contains("╷0.0"));
        assert!(text.contains("Taxon_A · Lengths 0.12"));
        assert!(text.contains("? Help"));

        app.view.root_on(app.cursor).expect("Failed rerooting");
        let text = screen(&mut app, &palette);
        assert!(text.contains("rerooted"));

        app.show_sidebar = true;
        let text = screen(&mut app, &palette);
        assert!(text.contains("Trees · Clades"));
        assert!(text.contains("▌Tree 1"));

        // The whole help fits a standard 80x24 terminal.
        app.overlay = Overlay::Help;
        let text = screen(&mut app, &palette);
        assert!(text.contains("Keyboard shortcuts"));
        assert!(text.contains("Quit"));
    }

    #[test]
    fn test_render_fold_and_layers() {
        let files = [PathBuf::from("tests/files/trees/simple.nwk")];
        let mut app = app(&files);
        let palette = palette();
        app.toggle_fold();
        let text = screen(&mut app, &palette);
        assert!(text.contains("▸ [+2 tips] 98"));
        assert!(text.contains("Clade · 2 tips (folded) · Support 98"));

        app.layers.toggle(0);
        app.layers.toggle(1);
        let text = screen(&mut app, &palette);
        assert!(!text.contains("85┼"));
        assert!(text.contains("Taxon_C  0.1"));
    }

    #[test]
    fn test_render_large_tree() {
        let files = [PathBuf::from("tests/files/trees/caterpillar.nwk")];
        let mut app = app(&files);
        let palette = palette();
        app.move_row(1500);
        let text = screen(&mut app, &palette);
        assert!(text.contains("2,000 tips"));
        assert!(text.contains(" t1501"));
        assert!(text.contains('┤') || text.contains('├'));
    }

    #[test]
    fn test_render_files() {
        let files = [
            PathBuf::from("tests/files/trees/simple.nwk"),
            PathBuf::from("tests/files/trees/multi.nwk"),
        ];
        let mut app = app(&files);
        let palette = palette();
        let text = screen(&mut app, &palette);
        assert!(text.contains("(1/2)"));
        assert!(text.contains("Ctrl+P Files"));
        app.overlay = Overlay::Files {
            query: String::new(),
            index: 0,
        };
        let text = screen(&mut app, &palette);
        assert!(text.contains("Open file (2)"));
    }

    #[test]
    fn test_line_chars() {
        assert_eq!(line_char(UP | DOWN | LEFT), Some('┤'));
        assert_eq!(line_char(DOWN | RIGHT), Some('┌'));
        assert_eq!(line_char(LEFT), Some('─'));
        assert_eq!(line_char(0), None);
        assert_eq!(nice_step(0.07), 0.1);
        assert_eq!(nice_step(0.3), 0.5);
    }
}
