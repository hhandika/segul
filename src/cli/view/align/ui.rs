//! Draws the viewer, laid out like the SEGUI alignment viewer: toolbar,
//! overview, legend, names and residue grid, then the status bar.

use ratatui::Frame;
use ratatui::buffer::Buffer;
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::style::{Color, Modifier, Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, BorderType, Clear, Paragraph, Wrap};

use crate::cli::view::thousands;
use crate::core::align::view::{WindowStat, is_gap};

use super::app::{App, Overlay, SIDEBAR_LINES, SidebarTab};
use super::palette::{ColorMode, Palette};

const MIN_WIDTH: u16 = 40;
const MIN_HEIGHT: u16 = 10;
// Rows of residues sampled by the overview; two per terminal row.
const OVERVIEW_ROWS: usize = 4;
const SIDEBAR_WIDTH: u16 = 24;
// Columns between ruler labels, at least.
const RULER_GAP: usize = 8;
// Partition suggestions shown above the go-to prompt, at most.
const GOTO_MATCHES: usize = 8;
const RULER_STEPS: [usize; 13] = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000];

pub(super) fn render(frame: &mut Frame, app: &mut App, palette: &Palette) {
    let area = frame.area();
    if area.width < MIN_WIDTH || area.height < MIN_HEIGHT {
        let note = format!("Enlarge the terminal to at least {MIN_WIDTH}x{MIN_HEIGHT}");
        frame.render_widget(Paragraph::new(note).wrap(Wrap { trim: true }), area);
        return;
    }
    let show_overview = area.height >= 16;
    let [toolbar, overview, legend, main, status] = Layout::vertical([
        Constraint::Length(1),
        Constraint::Length(if show_overview { 3 } else { 0 }),
        Constraint::Length(1),
        Constraint::Min(3),
        Constraint::Length(1),
    ])
    .areas(area);

    let [sidebar, names, divider, grid] = Layout::horizontal([
        Constraint::Length(if app.show_sidebar { SIDEBAR_WIDTH } else { 0 }),
        Constraint::Length(app.name_width.min(main.width / 3)),
        Constraint::Length(1),
        Constraint::Min(1),
    ])
    .areas(main);
    // The grid column starts with the partition track, when there is one,
    // and the ruler.
    let track_height = if app.partitions.is_some() { 1 } else { 0 };
    let [track, ruler, cells] = Layout::vertical([
        Constraint::Length(track_height),
        Constraint::Length(1),
        Constraint::Min(1),
    ])
    .areas(grid);
    let [track_header, names_header, name_rows] = Layout::vertical([
        Constraint::Length(track_height),
        Constraint::Length(1),
        Constraint::Min(1),
    ])
    .areas(names);

    app.grid_rows = cells.height as usize;
    app.grid_cols = cells.width as usize;
    app.clamp_scroll();
    app.areas.grid = cells;
    app.areas.names = name_rows;
    app.areas.divider = divider;
    app.areas.overview = overview;
    app.areas.track = track;
    app.areas.sidebar = if app.show_sidebar {
        sidebar_block().inner(sidebar)
    } else {
        Rect::default()
    };

    render_toolbar(frame, toolbar, app);
    if show_overview {
        render_overview(frame.buffer_mut(), overview, app, palette);
    }
    render_legend(frame, legend, palette);
    if app.show_sidebar {
        render_sidebar(frame, sidebar, app, palette);
    }
    let header = format!("Taxa ({})", count(app.view.rows()));
    frame.render_widget(Span::from(header).dim(), names_header);
    render_names(frame.buffer_mut(), name_rows, app, palette);
    render_divider(frame.buffer_mut(), divider);
    if let Some(partitions) = &app.partitions {
        let header = format!("Partitions ({})", count(partitions.len()));
        frame.render_widget(Span::from(header).dim(), track_header);
        render_track(frame.buffer_mut(), track, app, palette);
    }
    render_ruler(frame.buffer_mut(), ruler, app, palette);
    render_grid(frame.buffer_mut(), cells, app, palette);
    render_status(frame, status, app, palette);

    match &app.overlay {
        Overlay::Help => render_help(frame, area),
        Overlay::Files { query, index } => render_files(frame, area, app, query, *index, palette),
        Overlay::GoTo { input, index } => {
            render_goto_matches(frame, status, app, input, *index, palette)
        }
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
    let mut summary = format!(
        "{} taxa · {} sites · {}",
        count(view.rows()),
        count(view.length()),
        if view.is_dna() { "DNA" } else { "Amino acid" }
    );
    if !view.is_aligned() {
        summary.push_str(" · unaligned");
    }
    let start = app.window() * app.window_size;
    let window_end = (start + app.window_size).min(view.length());
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
    spans.extend([
        Span::from("  "),
        Span::from(summary),
        Span::from("  │ ").dim(),
        Span::from(app.color.label()),
        Span::from(if app.letters { " · Letters" } else { "" }),
        Span::from(format!(" · Zoom {}", app.zoom)),
        Span::from("  │ ").dim(),
        Span::from(format!(
            "Window {}–{} ({} of {})",
            count(start + 1),
            count(window_end),
            count(app.window() + 1),
            count(app.window_count()),
        )),
    ]);
    frame.render_widget(Line::from(spans), area);
}

// Two terminal rows of half-block cells, each showing two sampled rows,
// then a line marking the sites in view.
fn render_overview(buf: &mut Buffer, area: Rect, app: &App, palette: &Palette) {
    let view = &app.view;
    let width = area.width as usize;
    let rows = view.rows();
    let picked: Vec<usize> = (0..OVERVIEW_ROWS.min(rows))
        .map(|i| i * rows / OVERVIEW_ROWS.min(rows))
        .collect();
    let sample = view.sample(&picked, width);
    let color_at = |row: usize, x: usize| -> Color {
        match sample.get(row * width + x) {
            Some(&code) if code != 0 => palette.residue(code),
            _ => Color::Reset,
        }
    };
    for line in 0..2 {
        for x in 0..width {
            let top = color_at(line * 2, x);
            let bottom = if line * 2 + 1 < picked.len() {
                color_at(line * 2 + 1, x)
            } else {
                Color::Reset
            };
            buf[(area.x + x as u16, area.y + line as u16)]
                .set_char('▀')
                .set_fg(top)
                .set_bg(bottom);
        }
    }
    let length = view.length().max(1);
    let from = app.left_site * width / length;
    let to = ((app.left_site + app.grid_sites()) * width).div_ceil(length);
    let y = area.y + 2;
    for x in 0..width {
        let in_view = x >= from && x < to.max(from + 1);
        let cell = &mut buf[(area.x + x as u16, y)];
        if in_view {
            cell.set_char('━').set_fg(palette.accent);
        } else {
            cell.set_char('─').set_style(Style::new().dim());
        }
    }
}

fn render_legend(frame: &mut Frame, area: Rect, palette: &Palette) {
    let mut spans = Vec::new();
    let mut push = |label: &'static str, color: Color| {
        spans.push(Span::from("  ").bg(color));
        spans.push(Span::from(format!(" {label}  ")));
    };
    for &(label, color) in &palette.legend {
        push(label, color);
    }
    push("Gap", palette.gap);
    push("Other", palette.fallback);
    frame.render_widget(Line::from(spans), area);
}

fn sidebar_block() -> Block<'static> {
    Block::bordered().border_type(BorderType::Rounded)
}

fn render_sidebar(frame: &mut Frame, area: Rect, app: &App, palette: &Palette) {
    let is_partitions = app.partitions.is_some() && app.sidebar_tab == SidebarTab::Partitions;
    let title = match &app.partitions {
        Some(_) => {
            let tab = |label: &'static str, active: bool| {
                if active {
                    Span::from(label).fg(palette.accent).bold()
                } else {
                    Span::from(label).dim()
                }
            };
            Line::from(vec![
                Span::from(" "),
                tab("Windows", !is_partitions),
                Span::from(" · ").dim(),
                tab("Partitions", is_partitions),
                Span::from(" "),
            ])
        }
        None => Line::from(" Windows "),
    };
    let block = sidebar_block().title(title);
    let inner = block.inner(area);
    frame.render_widget(block, area);
    let shown = (inner.height as usize / SIDEBAR_LINES).max(1);
    let first = app.sidebar_first(shown);
    let entries: Vec<(String, bool, WindowStat)> = match &app.partitions {
        Some(track) if is_partitions => {
            let current = app.block();
            track
                .blocks
                .iter()
                .enumerate()
                .skip(first)
                .take(shown)
                .map(|(i, b)| {
                    let label = format!("{} · {}", b.name, count(b.end - b.start));
                    (
                        label,
                        Some(i) == current,
                        app.view.range_stats(b.start, b.end),
                    )
                })
                .collect()
        }
        _ => app
            .view
            .window_stats(app.window_size, first, shown)
            .into_iter()
            .enumerate()
            .map(|(i, stat)| {
                let label = format!("{}–{}", count(stat.start + 1), count(stat.end));
                (label, first + i == app.window(), stat)
            })
            .collect(),
    };
    let bar_width = inner.width.saturating_sub(1) as usize;
    let mut lines = Vec::with_capacity(entries.len() * SIDEBAR_LINES);
    for (label, is_current, stat) in entries {
        let marker = if is_current { "▌" } else { " " };
        let label = truncate(&format!("{marker}{label}"), inner.width as usize);
        lines.push(if is_current {
            Line::from(label).fg(palette.accent).bold()
        } else {
            Line::from(label)
        });
        lines.push(
            Line::from(format!(
                " Gaps {:.0}% · Var {:.0}%",
                stat.gap_fraction * 100.0,
                stat.variable_fraction * 100.0
            ))
            .dim(),
        );
        let filled = (stat.variable_fraction * bar_width as f64).round() as usize;
        lines.push(Line::from(vec![
            Span::from(" "),
            Span::from("█".repeat(filled)).fg(palette.tertiary),
            Span::from("░".repeat(bar_width - filled)).dim(),
        ]));
    }
    frame.render_widget(Paragraph::new(lines), inner);
}

// Partition bands in alternating colors, each named where it fits.
fn render_track(buf: &mut Buffer, area: Rect, app: &App, palette: &Palette) {
    let Some(track) = &app.partitions else {
        return;
    };
    let zoom = app.zoom;
    let sites = (area.width as usize / zoom).min(app.site_range().1 - app.left_site);
    let (start, end) = (app.left_site, app.left_site + sites);
    let current = app.block();
    let first = track.blocks.partition_point(|b| b.end <= start);
    for (offset, block) in track.blocks_in(start, end).iter().enumerate() {
        let index = first + offset;
        let bg = if index.is_multiple_of(2) {
            palette.mono
        } else {
            palette.gap
        };
        let mut style = Style::new().bg(bg).fg(palette.on_residue);
        if Some(index) == current {
            style = style.bold();
        }
        let from = (block.start.max(start) - start) * zoom;
        let to = (block.end.min(end) - start) * zoom;
        let width = to - from;
        let mut label = String::from(if Some(index) == current { "▸" } else { "" });
        label.push_str(&block.name);
        let text = format!("{:width$}", truncate(&label, width), width = width);
        buf.set_stringn(area.x + from as u16, area.y, text, width, style);
    }
}

fn render_names(buf: &mut Buffer, area: Rect, app: &App, palette: &Palette) {
    let matches = match &app.overlay {
        Overlay::Find(_) => app.find_matches(&app.find_query),
        _ => Vec::new(),
    };
    let width = area.width as usize;
    for y in 0..area.height as usize {
        let row = app.top_row + y;
        if row >= app.view.rows() {
            break;
        }
        let mut style = Style::new();
        if app.selected.contains(&row) {
            style = style.bg(palette.wash(palette.gap, 0.35));
        }
        if row == app.row {
            style = style.bold().fg(palette.accent);
        }
        if matches.binary_search(&row).is_ok() {
            style = style.underlined();
        }
        let name = truncate(app.view.name(row), width);
        let line = format!("{:width$}", name, width = width);
        buf.set_stringn(area.x, area.y + y as u16, line, width, style);
    }
}

fn render_divider(buf: &mut Buffer, area: Rect) {
    for y in area.top()..area.bottom() {
        buf[(area.x, y)].set_char('│').set_style(Style::new().dim());
    }
}

fn render_ruler(buf: &mut Buffer, area: Rect, app: &App, palette: &Palette) {
    // Room for the widest label, its tick, and a space.
    let gap = RULER_GAP.max(count(app.view.length()).chars().count() + 2);
    let step = RULER_STEPS
        .iter()
        .copied()
        .find(|s| s * app.zoom >= gap)
        .unwrap_or(*RULER_STEPS.last().unwrap_or(&1));
    let sites = (area.width as usize / app.zoom).min(app.site_range().1 - app.left_site);
    let mut free_from = 0usize;
    for i in 0..sites {
        let site = app.left_site + i;
        let x = i * app.zoom;
        // Labels drawn before take precedence over boundary marks. Partition
        // starts replace window marks when the alignment has partitions.
        let is_boundary = match &app.partitions {
            Some(track) => track
                .blocks
                .binary_search_by_key(&site, |b| b.start)
                .is_ok(),
            None => site.is_multiple_of(app.window_size),
        };
        if site > 0 && is_boundary && x >= free_from {
            buf[(area.x + x as u16, area.y)]
                .set_char('┃')
                .set_fg(palette.tertiary);
            free_from = free_from.max(x + 1);
        }
        let position = site + 1;
        if (position.is_multiple_of(step) || site == 0) && x >= free_from {
            let label = format!("╵{}", count(position));
            let room = area.width as usize - x;
            buf.set_stringn(area.x + x as u16, area.y, &label, room, Style::new().dim());
            free_from = x + label.chars().count() + 1;
        }
    }
}

fn render_grid(buf: &mut Buffer, area: Rect, app: &App, palette: &Palette) {
    let view = &app.view;
    let zoom = app.zoom;
    let sites = (area.width as usize / zoom).min(app.site_range().1 - app.left_site);
    let rows = (area.height as usize).min(view.rows().saturating_sub(app.top_row));
    let tile = view.tile(app.top_row, rows, app.left_site, sites);
    let reference = match app.color {
        ColorMode::Reference => view.tile(0, 1, app.left_site, sites),
        _ => vec![0; sites],
    };
    let letter_at = zoom / 2;
    let starts: Vec<usize> = match &app.partitions {
        Some(track) => track
            .blocks_in(app.left_site, app.left_site + sites)
            .iter()
            .map(|b| b.start)
            .filter(|&s| s > 0)
            .collect(),
        None => Vec::new(),
    };
    for y in 0..rows {
        let row = app.top_row + y;
        let is_selected = app.selected.contains(&row);
        for i in 0..sites {
            let code = tile[y * sites + i];
            let site = app.left_site + i;
            let (mut bg, mut fg, mut letter) =
                cell_style(app, palette, code, row, &reference[i..=i]);
            if code == 0 {
                letter = ' ';
            }
            if is_selected {
                bg = palette.wash(bg, 0.25);
            }
            if site == app.site && row != app.row {
                bg = palette.wash(bg, 0.18);
            }
            if !app.letters && letter != '.' {
                letter = ' ';
            }
            let is_cursor = row == app.row && site == app.site;
            if is_cursor {
                (bg, fg) = (fg, bg);
            }
            let is_start = starts.binary_search(&site).is_ok();
            for dx in 0..zoom {
                let x = area.x + (i * zoom + dx) as u16;
                let cell = &mut buf[(x, area.y + y as u16)];
                cell.set_char(if dx == letter_at { letter } else { ' ' })
                    .set_bg(bg)
                    .set_fg(fg);
                // Zoomed cells have room left of the letter for a boundary.
                if is_start && dx == 0 && dx != letter_at {
                    cell.set_char('▏').set_fg(palette.tertiary);
                }
                if is_cursor {
                    cell.modifier.insert(Modifier::BOLD);
                }
            }
        }
    }
}

// Background, foreground, and letter of one residue.
fn cell_style(
    app: &App,
    palette: &Palette,
    code: u8,
    row: usize,
    reference: &[u8],
) -> (Color, Color, char) {
    if code == 0 {
        return (Color::Reset, Color::Reset, ' ');
    }
    let letter = code as char;
    match app.color {
        ColorMode::Mono if is_gap(code) => (palette.gap, palette.on_residue, letter),
        ColorMode::Mono => (palette.mono, palette.on_residue, letter),
        // The first row itself stays in full color.
        ColorMode::Reference if row > 0 && !is_gap(code) && reference == [code] => {
            (palette.identity, palette.on_identity, '.')
        }
        _ => (palette.residue(code), palette.on_residue, letter),
    }
}

fn render_status(frame: &mut Frame, area: Rect, app: &App, palette: &Palette) {
    let left = match &app.overlay {
        Overlay::GoTo { input, .. } => Line::from(vec![
            Span::from(" Go to: ").fg(palette.accent),
            Span::from(format!("{input}▏")),
            Span::from(if app.partitions.is_some() {
                "  site · partition · partition:site · taxon"
            } else {
                "  site · taxon · taxon:site"
            })
            .dim(),
        ]),
        Overlay::Find(input) => {
            let matches = app.find_matches(input);
            let at = matches.iter().position(|&r| r == app.row);
            let count = match (matches.len(), at) {
                (0, _) if !input.is_empty() => String::from("  No results"),
                (0, _) => String::new(),
                (n, Some(i)) => format!("  {} of {}", count(i + 1), count(n)),
                (n, None) => format!("  {} matches", count(n)),
            };
            Line::from(vec![
                Span::from(" Find taxon: ").fg(palette.accent),
                Span::from(format!("{input}▏")),
                Span::from(count).dim(),
                Span::from("  Enter next · Shift+Enter previous · Esc close").dim(),
            ])
        }
        _ => match &app.message {
            Some(message) => Line::from(format!(" {message}")).fg(palette.accent),
            None => Line::from(cursor_status(app)),
        },
    };
    let mut right = Vec::new();
    if !app.selected.is_empty() {
        right.push(Span::from(format!(
            "{} selected · ",
            count(app.selected.len())
        )));
    }
    if app.files.len() > 1 {
        right.push(Span::from("Ctrl+P Files · ").dim());
    }
    right.push(Span::from("? Help ").dim());
    let right = Line::from(right).right_aligned();
    frame.render_widget(left, area);
    frame.render_widget(right, area);
}

fn cursor_status(app: &App) -> String {
    let view = &app.view;
    if view.rows() == 0 {
        return String::from(" No sequences");
    }
    let residue = view.residue(app.row, app.site);
    let position = match (residue, view.ungapped_position(app.row, app.site)) {
        (None, _) => String::from("past end"),
        (_, Some(p)) => format!("ungapped {}", count(p)),
        (_, None) => String::from("gap"),
    };
    let residue = residue.map(|c| (c as char).to_string()).unwrap_or_default();
    let mut status = format!(
        " {} · site {} · {} · {}",
        view.name(app.row),
        count(app.site + 1),
        position,
        residue
    );
    if let (Some(track), Some(index)) = (&app.partitions, app.block()) {
        let block = &track.blocks[index];
        // Name the exact subset when several share the block.
        let name = match track.partition_at(app.site) {
            Some(p) if block.partitions.len() > 1 => &track.partitions[p].gene,
            _ => &block.name,
        };
        status.push_str(&format!(
            " · {} ({}/{}) site {} of {}",
            name,
            count(index + 1),
            count(track.len()),
            count(app.site - block.start + 1),
            count(block.end - block.start)
        ));
    }
    status
}

// Keys for Windows and Linux, keys for macOS, and the action. macOS
// terminals do not pass Command to programs, and Control+arrows switch
// spaces, so the Mac keys use Option and fn instead.
const HELP: &[(&str, &str, &str)] = &[
    ("←↑→↓", "←↑→↓", "Move"),
    ("Ctrl+← / →", "Option+← / →", "Jump 10 sites"),
    ("Home / End", "fn+← / →", "Row start / end"),
    ("Ctrl+Home / End", "Ctrl+fn+← / →", "First / last taxon"),
    ("PgUp / PgDn", "Option+↑ / ↓", "Page rows"),
    ("Alt+PgUp / PgDn", "Option+Shift+← / →", "Page sites"),
    (
        "[ / ], { / }",
        "[ / ], { / }",
        "Previous / next window, partition",
    ),
    (
        "Ctrl+G",
        "Ctrl+G",
        "Go to site, partition, partition:site, or taxon",
    ),
    (
        "Ctrl+F, F3",
        "Ctrl+F, fn+F3",
        "Find taxon; Enter or F3 next, Shift back",
    ),
    ("+ / - / Ctrl+0", "+ / - / Ctrl+0", "Zoom in / out / reset"),
    ("F2 or Ctrl+K C", "Ctrl+K C", "Cycle color mode"),
    ("F4 or Ctrl+K L", "Ctrl+K L", "Toggle letters"),
    (
        "Ctrl+B, Ctrl+K T",
        "Ctrl+B, Ctrl+K T",
        "Toggle sidebar / switch its list",
    ),
    ("< / >", "< / >", "Narrow / widen names"),
    (
        "Space, Shift+↑↓, Ctrl+A",
        "Space, Shift+↑↓, Ctrl+A",
        "Select row / extend / all",
    ),
    (
        "Ctrl+C / Ctrl+S",
        "Ctrl+C / Ctrl+S",
        "Copy / save selected as FASTA",
    ),
    (
        "Ctrl+PgUp / PgDn, Ctrl+P",
        "Ctrl+fn+↑ / ↓, Ctrl+P",
        "Previous / next / pick file",
    ),
    ("Mouse", "Mouse", "Wheel rows · Shift sites · Ctrl zoom"),
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
    palette: &Palette,
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

// Partitions matching the go-to text, in a box just above the status bar.
fn render_goto_matches(
    frame: &mut Frame,
    status: Rect,
    app: &App,
    input: &str,
    index: Option<usize>,
    palette: &Palette,
) {
    let Some(track) = &app.partitions else {
        return;
    };
    let matches = app.goto_matches(input);
    if matches.is_empty() {
        return;
    }
    let shown = matches
        .len()
        .min(GOTO_MATCHES)
        .min((status.y as usize).saturating_sub(2));
    let width = 60.min(status.width);
    let height = shown as u16 + 2;
    let popup = Rect::new(status.x, status.y.saturating_sub(height), width, height);
    let block = Block::bordered()
        .border_type(BorderType::Rounded)
        .title(format!(" Partitions ({}) ", count(matches.len())))
        .title_bottom(Line::from(" ↑↓ choose · Tab fill · Enter go ").right_aligned());
    let inner = block.inner(popup);
    frame.render_widget(Clear, popup);
    frame.render_widget(block, popup);
    // Keep the picked match in view.
    let first = index.map_or(0, |i| i.saturating_sub(shown.saturating_sub(1)));
    let name_width = (inner.width as usize).saturating_sub(26).max(8);
    let lines: Vec<Line> = matches
        .iter()
        .enumerate()
        .skip(first)
        .take(shown)
        .map(|(i, &b)| {
            let block = &track.blocks[b];
            let range = format!("{}–{}", count(block.start + 1), count(block.end));
            let sites = format!("{} sites", count(block.end - block.start));
            let name = truncate(&block.name, name_width);
            let line = Line::from(format!(" {name:name_width$} {range:>13} {sites:>9}"));
            if Some(i) == index {
                line.fg(palette.accent).bold().reversed()
            } else {
                line
            }
        })
        .collect();
    frame.render_widget(Paragraph::new(lines), inner);
}

fn centered(area: Rect, width: u16, height: u16) -> Rect {
    Rect::new(
        area.x + (area.width - width) / 2,
        area.y + (area.height - height) / 2,
        width,
        height,
    )
}

// Site and taxon counts with thousands separators.
fn count(value: usize) -> String {
    thousands(value as i128)
}

// Cuts `text` to `width` characters, ending with an ellipsis when cut.
fn truncate(text: &str, width: usize) -> String {
    if text.chars().count() <= width {
        return text.to_string();
    }
    let mut cut: String = text.chars().take(width.saturating_sub(1)).collect();
    cut.push('…');
    cut
}

#[cfg(test)]
mod test {
    use super::*;

    use std::path::PathBuf;

    use ratatui::Terminal;
    use ratatui::backend::TestBackend;

    use crate::core::align::view::AlignmentView;
    use crate::helper::types::{DataType, InputFmt};

    use segul_core::color::residue::ResidueScheme;

    use super::super::palette::Theme;

    fn screen(app: &mut App, palette: &Palette) -> String {
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
        let files = [PathBuf::from("tests/files/simple.fas")];
        let view = AlignmentView::open(&files[0], &InputFmt::Auto, &DataType::Ignore)
            .expect("Failed opening fixture");
        let mut app = App::new(&files, &InputFmt::Auto, &DataType::Ignore, view, 1000);
        let palette = Palette::new(ResidueScheme::Clustal, Theme::Dark, true, true);
        let text = screen(&mut app, &palette);
        assert!(text.contains("simple.fas"));
        assert!(text.contains("2 taxa · 6 sites · DNA"));
        assert!(text.contains("ABCD"));
        assert!(text.contains("AGTATG"));
        assert!(text.contains("ABCD · site 1 · ungapped 1 · A"));
        assert!(text.contains("? Help"));

        app.show_sidebar = true;
        let text = screen(&mut app, &palette);
        assert!(text.contains("Windows"));

        // The whole help fits a standard 80x24 terminal.
        app.overlay = Overlay::Help;
        let text = screen(&mut app, &palette);
        assert!(text.contains("Keyboard shortcuts"));
        assert!(text.contains("Quit"));
    }

    #[test]
    fn test_render_partitions() {
        let files = [PathBuf::from("tests/files/view/concat_charset.nex")];
        let view = AlignmentView::open(&files[0], &InputFmt::Auto, &DataType::Ignore)
            .expect("Failed opening fixture");
        let mut app = App::new(&files, &InputFmt::Auto, &DataType::Ignore, view, 1000);
        let palette = Palette::new(ResidueScheme::Clustal, Theme::Dark, true, true);
        let text = screen(&mut app, &palette);
        assert!(text.contains("Partitions (3)"));
        // Five-site bands fit only part of each name.
        assert!(text.contains("▸loc…locu…locu…"));
        assert!(text.contains("locus1 (1/3) site 1 of 5"));

        app.show_sidebar = true;
        let text = screen(&mut app, &palette);
        assert!(text.contains("Windows · Partitions"));
        assert!(text.contains("▌locus1 · 5"));

        app.overlay = Overlay::GoTo {
            input: String::from("loc"),
            index: Some(1),
        };
        let text = screen(&mut app, &palette);
        assert!(text.contains("Go to: loc"));
        assert!(text.contains("6–10"));
        assert!(text.contains("Enter go"));
    }

    #[test]
    fn test_render_files() {
        let files = [
            PathBuf::from("tests/files/simple.fas"),
            PathBuf::from("tests/files/simple_aa.fas"),
        ];
        let view = AlignmentView::open(&files[0], &InputFmt::Auto, &DataType::Ignore)
            .expect("Failed opening fixture");
        let mut app = App::new(&files, &InputFmt::Auto, &DataType::Ignore, view, 1000);
        let palette = Palette::new(ResidueScheme::Clustal, Theme::Dark, true, true);
        let text = screen(&mut app, &palette);
        assert!(text.contains("(1/2)"));
        assert!(text.contains("Ctrl+P Files"));

        app.overlay = Overlay::Files {
            query: String::new(),
            index: 0,
        };
        let text = screen(&mut app, &palette);
        assert!(text.contains("Open file (2)"));
        assert!(text.contains("Enter open"));
    }

    #[test]
    fn test_truncate() {
        assert_eq!(truncate("sample", 10), "sample");
        assert_eq!(truncate("sample_long", 6), "sampl…");
    }
}
