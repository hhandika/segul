//! Draws the summary viewer: toolbar, section tabs, the overview or a
//! table as rows or a plot, then the status bar.

use std::sync::OnceLock;

use ratatui::Frame;
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::style::{Color, Modifier, Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Bar, BarChart, BarGroup, Block, BorderType, Clear, Paragraph, Wrap};
use segul_core::color::Theme;
use segul_core::color::residue::{ResiduePalette, ResidueScheme};
use segul_core::helper::types::DataType;

use super::app::{App, Mode, PlotKind};
use super::data::{Column, Format, SummaryTable};
use crate::cli::view::align::palette::{is_truecolor, to_color};

const MIN_WIDTH: u16 = 40;
const MIN_HEIGHT: u16 = 10;
/// Accent for labels, bars, and messages, set once by [`set_accent`].
static ACCENT: OnceLock<Color> = OnceLock::new();
const NAME_WIDTH: [usize; 2] = [6, 28];
const MAX_BINS: usize = 20;
/// Widest histogram bar, so a few bins do not fill the screen.
const MAX_BAR_WIDTH: usize = 9;
// Eighths of a cell, for smooth bar ends.
const PARTIAL: [&str; 8] = ["", "▏", "▎", "▍", "▌", "▋", "▊", "▉"];

const HELP: &[(&str, &str)] = &[
    ("Tab / Shift+Tab", "Next or previous section"),
    ("1 to 9", "Go to a section"),
    ("v or F2", "Switch between table and plot"),
    ("h", "Switch the plot between bars and histogram"),
    ("Left / Right", "Choose a column"),
    ("s", "Sort by the column: highest, lowest, report order"),
    ("Up / Down", "Move the selection"),
    ("PgUp / PgDn", "Move a page"),
    ("Home / End", "First or last item"),
    ("Mouse wheel", "Scroll"),
    ("? or F1", "This help"),
    ("q or Esc", "Quit"),
];

/// Use the accent of the shared color library for `theme`. It reaches 4.5:1
/// on black (dark) or white (light) backgrounds. Terminals without 24-bit
/// color get the nearest xterm-256 color on light themes and ANSI cyan on
/// dark ones, which every terminal theme keeps readable.
pub(super) fn set_accent(theme: Theme) {
    let source = ResiduePalette::new(&ResidueScheme::Clustal, &theme, &DataType::Dna);
    let accent = match (theme, is_truecolor()) {
        (_, true) => to_color(source.accent, true),
        (Theme::Light, false) => to_color(source.accent, false),
        (Theme::Dark, false) => Color::Cyan,
    };
    let _ = ACCENT.set(accent);
}

fn accent() -> Color {
    ACCENT.get().copied().unwrap_or(Color::Cyan)
}

pub(super) fn render(frame: &mut Frame, app: &mut App) {
    let area = frame.area();
    if area.width < MIN_WIDTH || area.height < MIN_HEIGHT {
        let note = format!("Enlarge the terminal to at least {MIN_WIDTH}x{MIN_HEIGHT}");
        frame.render_widget(Paragraph::new(note).wrap(Wrap { trim: true }), area);
        return;
    }
    let [toolbar, tabs, main, status] = Layout::vertical([
        Constraint::Length(1),
        Constraint::Length(1),
        Constraint::Min(3),
        Constraint::Length(1),
    ])
    .areas(area);
    render_toolbar(frame, toolbar, app);
    render_tabs(frame, tabs, app);
    match (app.section, app.mode) {
        (0, _) => render_overview(frame, main, app),
        (_, Mode::Table) => render_table(frame, main, app),
        (_, Mode::Plot) => match app.plot {
            PlotKind::Bars => render_bars(frame, main, app),
            PlotKind::Histogram => render_histogram(frame, main, app),
        },
    }
    render_status(frame, status, app);
    if app.show_help {
        render_help(frame, area);
    }
}

fn render_toolbar(frame: &mut Frame, area: Rect, app: &App) {
    let line = Line::from(vec![
        Span::from(format!(" {}", app.file_name)).bold(),
        Span::from("  "),
        Span::from(app.data.kind),
        Span::from(" · ").dim(),
        Span::from(app.data.headline.as_str()),
    ]);
    frame.render_widget(Paragraph::new(line), area);
}

fn render_tabs(frame: &mut Frame, area: Rect, app: &App) {
    let mut spans = vec![Span::from(" ")];
    let mut tab = |index: usize, label: String| {
        let text = format!(" {} {} ", index + 1, label);
        spans.push(if index == app.section {
            Span::from(text).add_modifier(Modifier::REVERSED | Modifier::BOLD)
        } else {
            Span::from(text)
        });
        spans.push(Span::from(" "));
    };
    tab(0, String::from("Overview"));
    for (i, table) in app.data.tables.iter().enumerate() {
        tab(i + 1, format!("{} ({})", table.name, table.names.len()));
    }
    let [left, right] =
        Layout::horizontal([Constraint::Min(1), Constraint::Length(24)]).areas(area);
    frame.render_widget(Paragraph::new(Line::from(spans)), left);
    if app.section > 0 {
        let pick = |label: &'static str, on: bool| {
            if on {
                Span::from(label).fg(accent()).bold()
            } else {
                Span::from(label).dim()
            }
        };
        let mut mode = vec![
            pick("Table", app.mode == Mode::Table),
            Span::from(" │ ").dim(),
            pick("Plot", app.mode == Mode::Plot),
        ];
        if app.mode == Mode::Plot {
            mode.push(Span::from(match app.plot {
                PlotKind::Bars => " (bars)",
                PlotKind::Histogram => " (histogram)",
            }));
        }
        frame.render_widget(Paragraph::new(Line::from(mode).right_aligned()), right);
    }
}

fn render_overview(frame: &mut Frame, area: Rect, app: &App) {
    let data = &app.data;
    let key_width = data
        .overview
        .iter()
        .map(|(k, _)| k.chars().count())
        .max()
        .unwrap_or(0);
    let value_width = data
        .overview
        .iter()
        .map(|(_, v)| v.chars().count())
        .max()
        .unwrap_or(0);
    // Wide enough for the longest line and the borders, leaving the
    // charts at least a third of the screen.
    let width = ((key_width + value_width + 5) as u16).min(area.width * 2 / 3);
    let [values, charts] =
        Layout::horizontal([Constraint::Length(width), Constraint::Min(10)]).areas(area);
    let lines: Vec<Line> = data
        .overview
        .iter()
        .map(|(key, value)| {
            Line::from(vec![
                Span::from(format!(" {:width$}  ", key, width = key_width)).dim(),
                Span::from(value.as_str()),
            ])
        })
        .collect();
    let block = Block::bordered()
        .border_type(BorderType::Rounded)
        .title(" Overview ");
    frame.render_widget(
        Paragraph::new(lines)
            .block(block)
            .wrap(Wrap { trim: false }),
        values,
    );
    if data.charts.is_empty() {
        return;
    }
    let areas = Layout::vertical(
        data.charts
            .iter()
            .map(|c| Constraint::Max(c.values.len() as u16 + 2)),
    )
    .split(charts);
    for (chart, &chart_area) in data.charts.iter().zip(areas.iter()) {
        let block = Block::bordered()
            .border_type(BorderType::Rounded)
            .title(format!(" {} ", chart.title));
        let inner = block.inner(chart_area);
        let values: Vec<Option<f64>> = chart.values.iter().map(|&v| Some(v)).collect();
        let labels: Vec<&str> = chart.labels.iter().map(String::as_str).collect();
        let lines = bar_lines(&labels, &values, chart.format, inner.width as usize, None);
        frame.render_widget(Paragraph::new(lines).block(block), chart_area);
    }
}

fn render_table(frame: &mut Frame, area: Rect, app: &mut App) {
    let [title, header, rows] = Layout::vertical([
        Constraint::Length(1),
        Constraint::Length(1),
        Constraint::Min(1),
    ])
    .areas(area);
    app.page_rows = rows.height as usize;
    app.keep_row_visible();
    let Some((table, state)) = app.table() else {
        return;
    };
    let name_width = name_width(table);
    let widths: Vec<usize> = table.columns.iter().map(column_width).collect();
    let first = first_visible_column(&widths, name_width + 1, area.width as usize, state.column);
    let last = last_visible_column(&widths, name_width + 1, area.width as usize, first);
    let hidden = widths.len().saturating_sub(last + 1) + first;
    frame.render_widget(Paragraph::new(table_title(app, table, hidden)), title);
    let mut head = vec![Span::from(format!(" {:name_width$}", table.item)).bold()];
    for (c, column) in table.columns.iter().enumerate().take(last + 1).skip(first) {
        let arrow = match state.sort {
            Some((s, true)) if s == c => "↓",
            Some((s, false)) if s == c => "↑",
            _ => "",
        };
        let label = format!(
            "{:>width$}",
            format!("{arrow}{}", column.label),
            width = widths[c]
        );
        head.push(if c == state.column {
            Span::from(label).fg(accent()).bold()
        } else {
            Span::from(label).bold()
        });
    }
    frame.render_widget(Paragraph::new(Line::from(head)), header);
    let lines: Vec<Line> = state
        .order
        .iter()
        .enumerate()
        .skip(state.offset)
        .take(rows.height as usize)
        .map(|(position, &item)| {
            let mut spans = vec![Span::from(format!(
                " {:name_width$}",
                truncate(&table.names[item], name_width)
            ))];
            for (c, column) in table.columns.iter().enumerate().take(last + 1).skip(first) {
                let text = format!("{:>width$}", show(column, item), width = widths[c]);
                spans.push(if c == state.column {
                    Span::from(text).fg(accent())
                } else {
                    Span::from(text)
                });
            }
            let line = Line::from(spans);
            if position == state.row {
                line.add_modifier(Modifier::REVERSED)
            } else {
                line
            }
        })
        .collect();
    frame.render_widget(Paragraph::new(lines), rows);
}

fn render_bars(frame: &mut Frame, area: Rect, app: &mut App) {
    let [title, rows] = Layout::vertical([Constraint::Length(1), Constraint::Min(1)]).areas(area);
    app.page_rows = rows.height as usize;
    app.keep_row_visible();
    let Some((table, state)) = app.table() else {
        return;
    };
    let Some(column) = table.columns.get(state.column) else {
        return;
    };
    frame.render_widget(Paragraph::new(plot_title(app, table, column)), title);
    let visible: Vec<usize> = state
        .order
        .iter()
        .skip(state.offset)
        .take(rows.height as usize)
        .copied()
        .collect();
    let names: Vec<&str> = visible.iter().map(|&i| table.names[i].as_str()).collect();
    let values: Vec<Option<f64>> = visible.iter().map(|&i| column.values[i]).collect();
    let largest = column
        .values
        .iter()
        .flatten()
        .fold(0.0_f64, |a, &b| a.max(b));
    let highlight = state.row.checked_sub(state.offset);
    let lines = bar_lines_scaled(
        &names,
        &values,
        column.format,
        rows.width as usize,
        highlight,
        largest,
        name_width(table),
    );
    frame.render_widget(Paragraph::new(lines), rows);
}

fn render_histogram(frame: &mut Frame, area: Rect, app: &mut App) {
    let [title, chart, note] = Layout::vertical([
        Constraint::Length(1),
        Constraint::Min(3),
        Constraint::Length(1),
    ])
    .areas(area);
    let Some((table, _)) = app.table() else {
        return;
    };
    let Some(column) = app.column() else {
        return;
    };
    frame.render_widget(Paragraph::new(plot_title(app, table, column)), title);
    let values: Vec<f64> = column.values.iter().flatten().copied().collect();
    if values.is_empty() {
        frame.render_widget(Paragraph::new(" No values to plot").dim(), chart);
        return;
    }
    let max_bins = (chart.width as usize / 4).clamp(1, MAX_BINS);
    let bins = Histogram::of(&values, max_bins);
    let selected = app
        .selected()
        .and_then(|item| column.values[item])
        .map(|value| bins.bin_of(value));
    let bar_width = (chart.width as usize / bins.counts.len())
        .saturating_sub(1)
        .clamp(1, MAX_BAR_WIDTH) as u16;
    let bars: Vec<Bar> = bins
        .counts
        .iter()
        .enumerate()
        .map(|(i, &count)| {
            let style = if Some(i) == selected {
                Style::new().fg(accent()).add_modifier(Modifier::REVERSED)
            } else {
                Style::new().fg(accent())
            };
            Bar::default()
                .value(count)
                .label(Line::from(bins.label(i, column.format)))
                .style(style)
        })
        .collect();
    let widget = BarChart::default()
        .data(BarGroup::default().bars(&bars))
        .bar_width(bar_width)
        .bar_gap(1)
        .value_style(Style::new().add_modifier(Modifier::REVERSED));
    frame.render_widget(widget, chart);
    if let (Some(item), Some(bin)) = (app.selected(), selected) {
        let text = format!(
            " Selected {}: {} (range {})",
            table.names[item],
            show(column, item),
            bins.label(bin, column.format)
        );
        frame.render_widget(Paragraph::new(text).dim(), note);
    }
}

fn render_status(frame: &mut Frame, area: Rect, app: &App) {
    let left = match (&app.message, app.table()) {
        (Some(message), _) => Line::from(format!(" {message}")).fg(accent()),
        (None, Some((table, state))) => match (app.selected(), app.column()) {
            (Some(item), Some(column)) => Line::from(format!(
                " {} · {}: {} · {} of {}",
                table.names[item],
                column.label,
                show(column, item),
                state.row + 1,
                state.order.len()
            )),
            _ => Line::from(format!(" No {}", table.name.to_lowercase())),
        },
        (None, None) => Line::from(" Tab to see the tables").dim(),
    };
    let [left_area, right_area] =
        Layout::horizontal([Constraint::Min(1), Constraint::Length(10)]).areas(area);
    frame.render_widget(Paragraph::new(left), left_area);
    frame.render_widget(
        Paragraph::new(Line::from("? Help ").right_aligned()).dim(),
        right_area,
    );
}

fn render_help(frame: &mut Frame, area: Rect) {
    let width = 70.min(area.width);
    let height = (HELP.len() as u16 + 2).min(area.height);
    let popup = centered(area, width, height);
    let lines: Vec<Line> = HELP
        .iter()
        .map(|(keys, action)| {
            Line::from(vec![
                Span::from(format!(" {:18}", keys)).bold(),
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

/// Name, order, and keys of the open table. [hidden] columns do not fit
/// and scroll into view with the arrow keys.
fn table_title<'a>(app: &App, table: &'a SummaryTable, hidden: usize) -> Line<'a> {
    let state = &app.tables[app.section - 1];
    let order = match state.sort {
        Some((c, true)) => format!("sorted by {}, highest first", table.columns[c].label),
        Some((c, false)) => format!("sorted by {}, lowest first", table.columns[c].label),
        None => String::from("report order"),
    };
    let more = match hidden {
        0 => String::new(),
        1 => String::from(" · 1 more column"),
        n => format!(" · {n} more columns"),
    };
    Line::from(vec![
        Span::from(format!(" {}", table.name)).bold(),
        Span::from(format!(" · {order}{more} · ← → column · s sort · v plot")).dim(),
    ])
}

fn plot_title<'a>(app: &App, table: &'a SummaryTable, column: &'a Column) -> Line<'a> {
    let mut values: Vec<f64> = column.values.iter().flatten().copied().collect();
    values.sort_by(f64::total_cmp);
    let stats = match (values.first(), values.last()) {
        (Some(&min), Some(&max)) => format!(
            " · min {} · median {} · max {}",
            column.format.show(min),
            column.format.show(median(&values)),
            column.format.show(max)
        ),
        _ => String::new(),
    };
    let hint = match app.plot {
        PlotKind::Bars => " · h histogram · v table",
        PlotKind::Histogram => " · h bars · v table",
    };
    Line::from(vec![
        Span::from(format!(
            " {} by {}",
            column.label,
            table.item.to_lowercase()
        ))
        .bold(),
        Span::from(stats),
        Span::from(hint).dim(),
    ])
}

/// Bars scaled to the largest of [values], each line a label, the bar, and
/// the value. [highlight] marks one line.
fn bar_lines<'a>(
    labels: &[&'a str],
    values: &[Option<f64>],
    format: Format,
    width: usize,
    highlight: Option<usize>,
) -> Vec<Line<'a>> {
    let largest = values.iter().flatten().fold(0.0_f64, |a, &b| a.max(b));
    let label_width = labels
        .iter()
        .map(|l| l.chars().count())
        .max()
        .unwrap_or(0)
        .clamp(1, NAME_WIDTH[1]);
    bar_lines_scaled(
        labels,
        values,
        format,
        width,
        highlight,
        largest,
        label_width,
    )
}

fn bar_lines_scaled<'a>(
    labels: &[&'a str],
    values: &[Option<f64>],
    format: Format,
    width: usize,
    highlight: Option<usize>,
    largest: f64,
    label_width: usize,
) -> Vec<Line<'a>> {
    let texts: Vec<String> = values
        .iter()
        .map(|v| v.map_or_else(|| String::from("-"), |v| format.show(v)))
        .collect();
    let value_width = texts.iter().map(|t| t.len()).max().unwrap_or(1);
    let bar_width = width.saturating_sub(label_width + value_width + 4);
    labels
        .iter()
        .zip(values.iter().zip(texts))
        .enumerate()
        .map(|(i, (label, (value, text)))| {
            let bar = value.map_or_else(String::new, |v| bar(v, largest, bar_width));
            let pad = bar_width.saturating_sub(bar.chars().count());
            let line = Line::from(vec![
                Span::from(format!(" {:label_width$} ", truncate(label, label_width))),
                Span::from(bar).fg(accent()),
                Span::from(" ".repeat(pad + 1)),
                Span::from(format!("{:>value_width$}", text)),
            ]);
            if Some(i) == highlight {
                line.add_modifier(Modifier::REVERSED)
            } else {
                line
            }
        })
        .collect()
}

/// A bar of [value] out of [largest] in [width] cells, with eighths at
/// the end.
fn bar(value: f64, largest: f64, width: usize) -> String {
    if largest <= 0.0 || value <= 0.0 || width == 0 {
        return String::new();
    }
    let eighths = ((value / largest).min(1.0) * width as f64 * 8.0).round() as usize;
    let mut text = "█".repeat(eighths / 8);
    text.push_str(PARTIAL[eighths % 8]);
    text
}

fn show(column: &Column, item: usize) -> String {
    column.values[item].map_or_else(|| String::from("-"), |v| column.format.show(v))
}

fn name_width(table: &SummaryTable) -> usize {
    table
        .names
        .iter()
        .map(|n| n.chars().count())
        .max()
        .unwrap_or(0)
        .max(table.item.len())
        .clamp(NAME_WIDTH[0], NAME_WIDTH[1])
}

fn column_width(column: &Column) -> usize {
    let longest = column
        .values
        .iter()
        .map(|v| v.map_or(1, |v| column.format.show(v).len()))
        .max()
        .unwrap_or(1);
    // Room for the sort arrow and a gap before the column.
    longest.max(column.label.chars().count() + 1) + 2
}

/// First value column to draw so the chosen one fits beside the names.
fn first_visible_column(widths: &[usize], used: usize, width: usize, chosen: usize) -> usize {
    let mut first = 0;
    while first < chosen {
        let end: usize = used + widths[first..=chosen].iter().sum::<usize>();
        if end <= width {
            break;
        }
        first += 1;
    }
    first
}

/// Last value column, from [first], that fits whole beside the names.
fn last_visible_column(widths: &[usize], used: usize, width: usize, first: usize) -> usize {
    let mut end = used;
    let mut last = first;
    for (c, &w) in widths.iter().enumerate().skip(first) {
        end += w;
        if end > width {
            break;
        }
        last = c;
    }
    last
}

fn median(sorted: &[f64]) -> f64 {
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 1 {
        sorted[mid]
    } else {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    }
}

fn truncate(text: &str, width: usize) -> String {
    if text.chars().count() <= width {
        return text.to_string();
    }
    let mut out: String = text.chars().take(width.saturating_sub(1)).collect();
    out.push('…');
    out
}

fn centered(area: Rect, width: u16, height: u16) -> Rect {
    let x = area.x + (area.width.saturating_sub(width)) / 2;
    let y = area.y + (area.height.saturating_sub(height)) / 2;
    Rect::new(x, y, width, height)
}

/// Counts of values in ranges of one round width.
struct Histogram {
    start: f64,
    step: f64,
    /// Whether every value is a whole number, so ranges end one below the
    /// next start.
    is_whole: bool,
    counts: Vec<u64>,
}

impl Histogram {
    /// Bins of [values] in at most [max_bins] ranges. Whole numbers with
    /// few distinct values get one bin each.
    fn of(values: &[f64], max_bins: usize) -> Self {
        let low = values.iter().copied().fold(f64::INFINITY, f64::min);
        let high = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let is_whole = values.iter().all(|v| v.fract() == 0.0);
        let mut step = if is_whole && high - low + 1.0 <= max_bins as f64 {
            1.0
        } else {
            round_step((high - low) / max_bins as f64)
        };
        if is_whole && step < 1.0 {
            step = 1.0;
        }
        let start = (low / step).floor() * step;
        let bins = ((high - start) / step).floor() as usize + 1;
        let mut counts = vec![0; bins];
        for &v in values {
            let bin = (((v - start) / step).floor() as usize).min(bins - 1);
            counts[bin] += 1;
        }
        Self {
            start,
            step,
            is_whole,
            counts,
        }
    }

    fn bin_of(&self, value: f64) -> usize {
        (((value - self.start) / self.step).floor().max(0.0) as usize).min(self.counts.len() - 1)
    }

    fn label(&self, bin: usize, format: Format) -> String {
        let from = self.start + bin as f64 * self.step;
        if self.is_whole && self.step == 1.0 {
            return format.show(from);
        }
        let to = if self.is_whole {
            from + self.step - 1.0
        } else {
            from + self.step
        };
        let show = |v: f64| match format {
            Format::Count => format.show(v),
            // Bounds of round steps need no forced decimal.
            Format::Decimal if v.fract() == 0.0 => format!("{v}"),
            Format::Decimal => format!("{v:.1}"),
        };
        format!("{}-{}", show(from), show(to))
    }
}

/// A width of 1, 2, or 5 times a power of ten, at least [raw].
fn round_step(raw: f64) -> f64 {
    if raw <= 0.0 {
        return 1.0;
    }
    let scale = 10_f64.powf(raw.log10().floor());
    [1.0, 2.0, 5.0, 10.0]
        .into_iter()
        .map(|s| s * scale)
        .find(|&s| s >= raw)
        .unwrap_or(10.0 * scale)
}

#[cfg(test)]
mod test {
    use super::*;

    use std::path::Path;

    use crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
    use ratatui::Terminal;
    use ratatui::backend::TestBackend;

    use super::super::data::SummaryData;

    fn screen(app: &mut App) -> String {
        let mut terminal =
            Terminal::new(TestBackend::new(100, 24)).expect("Failed creating terminal");
        terminal
            .draw(|frame| render(frame, app))
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

    fn press(app: &mut App, code: KeyCode) {
        app.handle_key(KeyEvent::new(code, KeyModifiers::NONE));
    }

    #[test]
    fn test_render_alignment() {
        let data = SummaryData::from_path(Path::new(
            "crates/segul-core/tests/files/plot/alignment-report.json",
        ))
        .unwrap();
        let mut app = App::new(data, String::from("alignment-report.json"));
        let text = screen(&mut app);
        assert!(text.contains("alignment-report.json"));
        assert!(text.contains("3 taxa · 4 loci · DNA"));
        assert!(text.contains("2 Loci (4)"));
        assert!(text.contains("gene_1"));
        assert!(text.contains("Informative"));
        assert!(text.contains("gene_1 · Taxa: 1 · 1 of 4"));

        press(&mut app, KeyCode::Char('v'));
        let text = screen(&mut app);
        assert!(text.contains("Taxa by locus"));
        assert!(text.contains("█"));

        press(&mut app, KeyCode::Char('h'));
        let text = screen(&mut app);
        assert!(text.contains("(histogram)"));
        assert!(text.contains("Selected gene_1"));

        press(&mut app, KeyCode::Char('1'));
        let text = screen(&mut app);
        assert!(text.contains("Loci with at least this share of taxa"));
        assert!(text.contains("GC content"));

        app.show_help = true;
        let text = screen(&mut app);
        assert!(text.contains("Keyboard shortcuts"));
    }

    #[test]
    fn test_render_reads() {
        let data = SummaryData::from_path(Path::new(
            "crates/segul-core/tests/files/plot/read-report.json",
        ))
        .unwrap();
        let mut app = App::new(data, String::from("read-report.json"));
        let text = screen(&mut app);
        assert!(text.contains("Read summary"));
        assert!(text.contains("Samples ("));
        // Columns scroll so the chosen one stays in view.
        for _ in 0..12 {
            press(&mut app, KeyCode::Right);
        }
        let text = screen(&mut app);
        assert!(!text.contains(" Reads "));
    }

    #[test]
    fn test_histogram() {
        let bins = Histogram::of(&[2.0, 3.0, 3.0, 4.0], 10);
        assert_eq!(bins.counts, vec![1, 2, 1]);
        assert_eq!(bins.label(0, Format::Count), "2");
        let values: Vec<f64> = (100..300).map(f64::from).collect();
        let bins = Histogram::of(&values, 10);
        assert_eq!(bins.counts.len(), 10);
        assert_eq!(bins.label(0, Format::Count), "100-119");
        assert_eq!(bins.bin_of(299.0), 9);
        let bins = Histogram::of(&[0.5, 1.2, 9.9], 5);
        assert_eq!(bins.label(0, Format::Decimal), "0-2");
        assert_eq!(bins.counts.iter().sum::<u64>(), 3);
    }

    #[test]
    fn test_bar() {
        assert_eq!(bar(10.0, 10.0, 4), "████");
        assert_eq!(bar(5.0, 10.0, 3), "█▌");
        assert_eq!(bar(0.0, 10.0, 4), "");
        assert_eq!(truncate("sample_long", 6), "sampl…");
    }
}
