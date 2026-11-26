//! Draws the read report viewer: toolbar, tabs, the open sample, the
//! sample list, the charts or tables of the open tab, then the status bar.

use ratatui::Frame;
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::style::{Color, Modifier, Style, Stylize};
use ratatui::symbols::Marker;
use ratatui::text::{Line, Span};
use ratatui::widgets::{
    Axis, Block, BorderType, Chart, Clear, Dataset, GraphType, Paragraph, Wrap,
};
use segul_core::color::residue::{ResiduePalette, ResidueScheme};
use segul_core::color::ui::ChartTokens;
use segul_core::color::{Palette, Rgb, Theme, best_text_on};
use segul_core::helper::types::DataType;
use segul_core::stats::fastq::qc::{QcStatus, overall_status, pair_qc_checks};

use super::super::ui::{
    accent, bar, centered, first_visible_column, last_visible_column, scale, selected_bar, truncate,
};
use super::app::{App, Popup, SpreadView, SummaryView, Tab, column_label, keep_visible};
use super::metric::{ReadMetric, ReadSample, compact};
use super::spread::{Points, SampleMatrix, SampleSpread};
use crate::cli::view::align::palette::{is_truecolor, to_color};
use crate::cli::view::thousands;

const MIN_WIDTH: u16 = 50;
const MIN_HEIGHT: u16 = 14;
/// Width of the sample list, and the narrowest screen that shows it.
const LIST_WIDTH: u16 = 32;
const LIST_MIN_SCREEN: u16 = 96;
/// Shortest chart block; shorter stacks show one chart at a time.
const MIN_SECTION_HEIGHT: u16 = 9;
const NAME_WIDTH: [usize; 2] = [6, 28];
const SEPARATOR: &str = "│";
const HEAT: &str = "█";
const DOT: &str = "●";
/// FastQC quality zones: poor below 20, reasonable below 28.
const POOR_QUALITY: f64 = 20.0;
const GOOD_QUALITY: f64 = 28.0;

const HELP: &[(&str, &str)] = &[
    ("Tab / Shift+Tab", "Next or previous tab"),
    ("1 to 8", "Go to a tab"),
    ("Alt+↑↓ or [ ]", "Previous or next sample in the list"),
    ("Up / Down", "Move in a table, or step through samples"),
    ("PgUp / PgDn", "Move a page"),
    ("Left / Right", "Choose a metric or column"),
    ("v", "Switch the summary: bars, table, heatmap"),
    ("s", "Sort by the column"),
    ("c", "Checks of the open sample and their ranges"),
    ("x", "Choose samples to compare"),
    ("b", "Show many samples as a band or heatmap"),
    ("r", "Toggle absolute (bp) vs relative (0-100%) position"),
    ("n / p", "Next or previous chart when one fits"),
    ("l", "Show or hide the sample list"),
    ("/", "Filter samples by name"),
    ("o / O", "Sort the sample list, reverse it"),
    ("y", "Copy the selected sequence"),
    ("? or F1", "This help"),
    ("q or Esc", "Quit"),
];

/// Terminal colors of statuses, series, and bases for the theme.
struct Colors {
    pass: Color,
    warn: Color,
    fail: Color,
    open: Color,
    outlier: Color,
    axis: Color,
    series: Vec<Color>,
    /// A, C, G, T, and N.
    bases: [Color; 5],
}

impl Colors {
    fn new(theme: Theme) -> Self {
        let truecolor = is_truecolor();
        let tokens = ChartTokens::new(&theme);
        let line = |rgb: Rgb| to_color(tokens.line(rgb), truecolor);
        // Deep cells of the other theme read as lines on this background.
        let other = match theme {
            Theme::Dark => Theme::Light,
            Theme::Light => Theme::Dark,
        };
        let residues = ResiduePalette::new(&ResidueScheme::Clustal, &other, &DataType::Dna);
        let base = |code: u8| line(residues.residue(code));
        let status = tokens.status;
        Self {
            pass: line(status.pass),
            warn: line(status.warn),
            fail: line(status.fail),
            open: line(status.open),
            outlier: line(status.outlier),
            axis: to_color(tokens.axis, truecolor),
            series: Palette::OkabeIto.colors(8).into_iter().map(line).collect(),
            bases: [
                base(b'A'),
                base(b'C'),
                base(b'G'),
                base(b'T'),
                to_color(tokens.axis, truecolor),
            ],
        }
    }

    fn status(&self, status: QcStatus) -> Color {
        match status {
            QcStatus::Pass => self.pass,
            QcStatus::Warn => self.warn,
            QcStatus::Fail => self.fail,
        }
    }

    fn series_at(&self, index: usize) -> Color {
        self.series[index % self.series.len()]
    }
}

pub(super) fn render(frame: &mut Frame, app: &mut App) {
    let area = frame.area();
    if area.width < MIN_WIDTH || area.height < MIN_HEIGHT {
        let note = format!("Enlarge the terminal to at least {MIN_WIDTH}x{MIN_HEIGHT}");
        frame.render_widget(Paragraph::new(note).wrap(Wrap { trim: true }), area);
        return;
    }
    let colors = Colors::new(app.theme);
    let [toolbar, tabs, sample, main, status] = Layout::vertical([
        Constraint::Length(1),
        Constraint::Length(1),
        Constraint::Length(1),
        Constraint::Min(5),
        Constraint::Length(1),
    ])
    .areas(area);
    render_toolbar(frame, toolbar, app, &colors);
    render_tabs(frame, tabs, app);
    render_sample_bar(frame, sample, app, &colors);
    let shows_list = app.show_list && app.samples.len() > 1 && area.width >= LIST_MIN_SCREEN;
    let content = if shows_list {
        let [list, content] =
            Layout::horizontal([Constraint::Length(LIST_WIDTH), Constraint::Min(10)]).areas(main);
        render_list(frame, list, app, &colors);
        content
    } else {
        main
    };
    match app.tab {
        Tab::Summary => render_summary(frame, content, app, &colors),
        Tab::Overrepresented => render_sequences(frame, content, app),
        Tab::Cohort => render_cohort(frame, content, app, &colors),
        _ => render_sections(frame, content, app, &colors),
    }
    render_status(frame, status, app);
    match app.popup {
        Some(Popup::Help) => render_help(frame, area),
        Some(Popup::Checks(row)) => render_checks_popup(frame, area, app, row, &colors),
        Some(Popup::Compare(row)) => render_compare_popup(frame, area, app, row, &colors),
        None => (),
    }
}

fn render_toolbar(frame: &mut Frame, area: Rect, app: &App, colors: &Colors) {
    let reads: u64 = app.samples.iter().map(|s| s.metrics().read_count).sum();
    let count = |status: QcStatus| app.samples.iter().filter(|s| s.status == status).count();
    let mut spans = vec![
        Span::from(format!(" {}", app.file_name)).bold(),
        Span::from("  Read summary"),
        Span::from(" · ").dim(),
        Span::from(format!(
            "{} samples · {} reads",
            thousands(app.samples.len() as i128),
            thousands(reads as i128)
        )),
    ];
    for status in [QcStatus::Pass, QcStatus::Warn, QcStatus::Fail] {
        spans.push(Span::from(" · ").dim());
        spans.push(Span::from(DOT).fg(colors.status(status)));
        spans.push(Span::from(format!(
            " {} {}",
            count(status),
            status.label().to_lowercase()
        )));
    }
    frame.render_widget(Paragraph::new(Line::from(spans)), area);
}

fn render_tabs(frame: &mut Frame, area: Rect, app: &App) {
    let mut spans = vec![Span::from(" ")];
    for (index, tab) in Tab::ALL.iter().enumerate() {
        let text = format!(" {} {} ", index + 1, tab.label());
        spans.push(if *tab == app.tab {
            Span::from(text).add_modifier(Modifier::REVERSED | Modifier::BOLD)
        } else {
            Span::from(text)
        });
        spans.push(Span::from(" "));
    }
    frame.render_widget(Paragraph::new(Line::from(spans)), area);
}

fn render_sample_bar(frame: &mut Frame, area: Rect, app: &App, colors: &Colors) {
    let sample = app.open_sample();
    let place = match app.view.iter().position(|&i| i == app.open) {
        Some(p) => format!("{} of {} listed", p + 1, app.view.len()),
        None => String::from("not listed"),
    };
    let line = Line::from(vec![
        Span::from(" "),
        Span::from(DOT).fg(colors.status(sample.status)),
        Span::from(format!(" {}", sample.name())).bold(),
        Span::from(format!("  {}", sample.status.label())).fg(colors.status(sample.status)),
        Span::from(format!(" · {}", sample.profile.label())).dim(),
        Span::from(format!(" · {place}")),
        Span::from(" · Alt+↑↓ or [ ] sample").dim(),
    ]);
    frame.render_widget(Paragraph::new(line), area);
}

/// Every listed sample with its status, reads, and Q30 share.
fn render_list(frame: &mut Frame, area: Rect, app: &mut App, colors: &Colors) {
    let block = Block::bordered()
        .border_type(BorderType::Rounded)
        .title(format!(" Samples ({}) ", app.samples.len()))
        .title_bottom(Line::from(" / filter · o sort ").dim());
    let inner = block.inner(area);
    frame.render_widget(block, area);
    let filter = app.filtering || !app.query.is_empty();
    let [head, rows] = Layout::vertical([
        Constraint::Length(if filter { 2 } else { 1 }),
        Constraint::Min(1),
    ])
    .areas(inner);
    let mut head_lines = vec![Line::from(format!(" {}", app.sort.describe(app.ascending))).dim()];
    if filter {
        let cursor = if app.filtering { "▏" } else { "" };
        head_lines.push(Line::from(vec![
            Span::from(" / ").fg(accent()),
            Span::from(format!("{}{cursor}", app.query)),
        ]));
    }
    frame.render_widget(Paragraph::new(head_lines), head);
    if app.view.is_empty() {
        frame.render_widget(Paragraph::new(" No samples match").dim(), rows);
        return;
    }
    let position = app.view.iter().position(|&i| i == app.open).unwrap_or(0);
    keep_visible(
        &mut app.list_offset,
        position,
        rows.height as usize,
        app.view.len(),
    );
    let name_width = (rows.width as usize).saturating_sub(16);
    let lines: Vec<Line> = app
        .view
        .iter()
        .skip(app.list_offset)
        .take(rows.height as usize)
        .map(|&i| {
            let sample = &app.samples[i];
            let text = format!(
                " {:name_width$} {:>6} {:>5.1}",
                truncate(sample.name(), name_width),
                compact(sample.metrics().read_count as f64),
                sample.metrics().q30_pct
            );
            let line = Line::from(vec![
                Span::from(format!(" {DOT}")).fg(colors.status(sample.status)),
                Span::from(text),
            ]);
            if i == app.open {
                line.add_modifier(Modifier::REVERSED)
            } else {
                line
            }
        })
        .collect();
    frame.render_widget(Paragraph::new(lines), rows);
}

fn render_summary(frame: &mut Frame, area: Rect, app: &mut App, colors: &Colors) {
    let checks = check_lines(app, colors, area.width);
    let [controls, samples, checks_area] = Layout::vertical([
        Constraint::Length(1),
        Constraint::Min(3),
        Constraint::Length(checks.len() as u16),
    ])
    .areas(area);
    let mut spans = vec![Span::from(" ")];
    for (i, view) in [SummaryView::Bars, SummaryView::Table, SummaryView::Heatmap]
        .into_iter()
        .enumerate()
    {
        if i > 0 {
            spans.push(Span::from(" │ ").dim());
        }
        spans.push(if view == app.summary_view {
            Span::from(view.label()).fg(accent()).bold()
        } else {
            Span::from(view.label()).dim()
        });
    }
    spans.push(Span::from("  v view").dim());
    match (app.summary_view, app.metric()) {
        (SummaryView::Bars, Some(metric)) => {
            spans.push(Span::from("   ◀ ").dim());
            spans.push(Span::from(metric.label()).bold());
            spans.push(Span::from(" ▶  ← → metric").dim());
        }
        (_, _) => {
            let order = match app.table_sort {
                Some((c, true)) => format!("sorted by {}, highest first", column_label(c)),
                Some((c, false)) => format!("sorted by {}, lowest first", column_label(c)),
                None => String::from("list order"),
            };
            spans.push(Span::from(format!("   {order} · ← → column · s sort")).dim());
        }
    }
    frame.render_widget(Paragraph::new(Line::from(spans)), controls);
    if app.view.is_empty() {
        frame.render_widget(
            Paragraph::new(" No samples match the filter.").dim(),
            samples,
        );
    } else {
        match app.summary_view {
            SummaryView::Bars => render_summary_bars(frame, samples, app, colors),
            SummaryView::Table => render_summary_table(frame, samples, app, colors, false),
            SummaryView::Heatmap => render_summary_table(frame, samples, app, colors, true),
        }
    }
    frame.render_widget(Paragraph::new(checks), checks_area);
}

fn render_summary_bars(frame: &mut Frame, area: Rect, app: &mut App, colors: &Colors) {
    let Some(metric) = app.metric() else {
        return;
    };
    let rows = app.summary_rows();
    let open = rows.iter().position(|&i| i == app.open).unwrap_or(0);
    app.page_rows = area.height as usize;
    keep_visible(
        &mut app.summary_offset,
        open,
        area.height as usize,
        rows.len(),
    );
    let value = |i: usize| metric.value_of(app.samples[i].metrics()).unwrap_or(0.0);
    let largest = rows.iter().map(|&i| value(i)).fold(0.0, f64::max);
    let name_width = name_width(rows.iter().map(|&i| app.samples[i].name()));
    let texts: Vec<String> = rows.iter().map(|&i| metric.short(value(i))).collect();
    let value_width = texts.iter().map(|t| t.len()).max().unwrap_or(1);
    let bar_width = (area.width as usize).saturating_sub(name_width + value_width + 7);
    let lines: Vec<Line> = rows
        .iter()
        .zip(&texts)
        .enumerate()
        .skip(app.summary_offset)
        .take(area.height as usize)
        .map(|(position, (&i, text))| {
            let sample = &app.samples[i];
            let cells = bar(value(i), largest, bar_width);
            let pad = " ".repeat(bar_width.saturating_sub(cells.chars().count()) + 1);
            let label = Span::from(format!(
                " {:name_width$} ",
                truncate(sample.name(), name_width)
            ));
            let text = Span::from(format!("{text:>value_width$}"));
            let dot = Span::from(format!(" {DOT}")).fg(colors.status(sample.status));
            // Reversing the whole line would hide the bar, so only the
            // label and value are reversed.
            if position == open {
                Line::from(vec![
                    dot,
                    label.add_modifier(Modifier::REVERSED),
                    Span::from(cells).style(selected_bar()),
                    Span::from(pad),
                    text.add_modifier(Modifier::REVERSED),
                ])
            } else {
                Line::from(vec![
                    dot,
                    label,
                    Span::from(cells).fg(accent()),
                    Span::from(pad),
                    text,
                ])
            }
        })
        .collect();
    frame.render_widget(Paragraph::new(lines), area);
}

/// Every metric of every listed sample, colored by status in columns with
/// a check, or as a [heatmap] scaled per column.
fn render_summary_table(
    frame: &mut Frame,
    area: Rect,
    app: &mut App,
    colors: &Colors,
    heatmap: bool,
) {
    let [header, body, legend] = Layout::vertical([
        Constraint::Length(1),
        Constraint::Min(1),
        Constraint::Length(1),
    ])
    .areas(area);
    let rows = app.summary_rows();
    let open = rows.iter().position(|&i| i == app.open).unwrap_or(0);
    app.page_rows = body.height as usize;
    keep_visible(
        &mut app.summary_offset,
        open,
        body.height as usize,
        rows.len(),
    );
    let name_width = name_width(rows.iter().map(|&i| app.samples[i].name()));
    let columns = ReadMetric::TABLE.len() + 1;
    let text_of = |i: usize, c: usize| -> String {
        match c {
            0 => format!("{DOT} {}", app.samples[i].status.label()),
            c => {
                let metric = ReadMetric::TABLE[c - 1];
                metric
                    .value_of(app.samples[i].metrics())
                    .map_or_else(|| String::from("-"), |v| metric.exact(v))
            }
        }
    };
    let widths: Vec<usize> = (0..columns)
        .map(|c| {
            let longest = rows
                .iter()
                .map(|&i| text_of(i, c).chars().count())
                .max()
                .unwrap_or(1);
            (longest + 1).max(column_label(c).chars().count() + 1) + 2
        })
        .collect();
    let used = name_width + 2;
    let first = first_visible_column(&widths, used, area.width as usize, app.column);
    let last = last_visible_column(&widths, used, area.width as usize, first);
    let mut head = vec![Span::from(format!(" {:name_width$} ", "Sample")).bold()];
    for (c, width) in widths.iter().enumerate().take(last + 1).skip(first) {
        let arrow = match app.table_sort {
            Some((s, true)) if s == c => "↓",
            Some((s, false)) if s == c => "↑",
            _ => "",
        };
        let label = format!(
            "{:>width$} ",
            format!("{arrow}{}", column_label(c)),
            width = width - 2
        );
        head.push(Span::from(SEPARATOR).dim());
        head.push(if c == app.column {
            Span::from(label).fg(accent()).bold()
        } else {
            Span::from(label).bold()
        });
    }
    frame.render_widget(Paragraph::new(Line::from(head)), header);
    let ranges: Vec<Option<(f64, f64)>> = (0..columns)
        .map(|c| {
            rows.iter()
                .filter_map(|&i| app.summary_value(i, c))
                .fold(None, |range, v| {
                    Some(range.map_or((v, v), |(lo, hi): (f64, f64)| (lo.min(v), hi.max(v))))
                })
        })
        .collect();
    let truecolor = is_truecolor();
    let lines: Vec<Line> = rows
        .iter()
        .enumerate()
        .skip(app.summary_offset)
        .take(body.height as usize)
        .map(|(position, &i)| {
            let sample = &app.samples[i];
            let selected = position == open;
            let name = Span::from(format!(
                " {:name_width$} ",
                truncate(sample.name(), name_width)
            ));
            let mut spans = vec![if selected {
                name.add_modifier(Modifier::REVERSED)
            } else {
                name
            }];
            for c in first..=last {
                spans.push(Span::from(SEPARATOR).dim());
                let width = widths[c] - 2;
                let text = format!("{:>width$} ", text_of(i, c));
                let value = app.summary_value(i, c);
                let span = match (c, value, ranges[c]) {
                    (0, _, _) => Span::from(text).fg(colors.status(sample.status)),
                    (_, Some(v), Some((lo, hi))) if heatmap => {
                        let color = app.colormap.at(scale(v, lo, hi));
                        Span::from(text)
                            .bg(to_color(color, truecolor))
                            .fg(to_color(best_text_on(color), truecolor))
                    }
                    (_, Some(v), _) => match ReadMetric::TABLE[c - 1].rule() {
                        Some(rule) => Span::from(text).fg(colors.status(rule.status_of(v))),
                        None => Span::from(text),
                    },
                    _ => Span::from(text).dim(),
                };
                spans.push(if selected && !heatmap {
                    span.add_modifier(Modifier::REVERSED)
                } else {
                    span
                });
            }
            Line::from(spans)
        })
        .collect();
    frame.render_widget(Paragraph::new(lines), body);
    let key = if heatmap {
        let mut key = vec![Span::from(" low ").dim()];
        key.extend(
            app.colormap
                .stops(8)
                .into_iter()
                .map(|color| Span::from(HEAT).fg(to_color(color, truecolor))),
        );
        key.push(Span::from(" high · scaled per column").dim());
        key
    } else {
        let mut key = vec![Span::from(" ")];
        for status in [QcStatus::Pass, QcStatus::Warn, QcStatus::Fail] {
            key.push(Span::from(DOT).fg(colors.status(status)));
            key.push(Span::from(format!(" {}  ", status.label())).dim());
        }
        key.push(Span::from("Values colored by their check").dim());
        key
    };
    frame.render_widget(Paragraph::new(Line::from(key)), legend);
}

/// Checks of the open sample under a heading, one colored item each,
/// wrapped between items to fit [width].
fn check_lines(app: &App, colors: &Colors, width: u16) -> Vec<Line<'static>> {
    let sample = app.open_sample();
    let mut lines = vec![Line::from(vec![
        Span::from(format!(" Checks for {}", sample.name())).bold(),
        Span::from(" · c ranges").dim(),
    ])];
    let mut line: Vec<Span> = vec![];
    let mut used = 0;
    for check in &sample.checks {
        let text = format!(" {}: {}  ", check.name, check.value_text());
        let item = text.chars().count() + 2;
        if used > 0 && used + item > width as usize {
            lines.push(Line::from(std::mem::take(&mut line)));
            used = 0;
        }
        line.push(Span::from(format!(" {DOT}")).fg(colors.status(check.status)));
        line.push(Span::from(text));
        used += item;
    }
    lines.push(Line::from(line));
    lines
}

/// Charts of a tab other than the summary, stacked when they fit and one
/// at a time when they do not.
fn render_sections(frame: &mut Frame, area: Rect, app: &mut App, colors: &Colors) {
    let sections = app.sections();
    app.is_stacked = area.height / sections as u16 >= MIN_SECTION_HEIGHT;
    let shown: Vec<usize> = if app.is_stacked {
        (0..sections).collect()
    } else {
        vec![app.section.min(sections - 1)]
    };
    let areas = Layout::vertical(shown.iter().map(|_| Constraint::Fill(1))).split(area);
    for (&section, &section_area) in shown.iter().zip(areas.iter()) {
        render_section(frame, section_area, app, colors, section, sections);
    }
}

fn render_section(
    frame: &mut Frame,
    area: Rect,
    app: &App,
    colors: &Colors,
    section: usize,
    sections: usize,
) {
    let sample = app.open_sample();
    let name = sample.name();
    let report = &sample.report;
    let (title, description, comparison) = match (app.tab, section) {
        (Tab::Quality, 0) => (
            if app.relative_pos && !report.relative_positions.is_empty() {
                format!("Relative-position quality: {name}")
            } else {
                format!("Per-base quality: {name}")
            },
            "Median, mean, and effective quality across read position, with the middle \
            half and 10th to 90th percentile of bases.",
            None,
        ),
        (Tab::Quality, 1) => (
            String::from("Per-read quality"),
            "Reads by their mean quality score.",
            Some(Comparison {
                points: |s| Points::indexed(&to_f64(&s.report.read_quality_distribution)),
                is_count: true,
                x_title: "Mean quality score",
                y_title: "Reads",
            }),
        ),
        (Tab::Quality, _) => (
            format!("Expected errors per read: {name}"),
            "Distribution of expected sequencing errors per read (sum of 10^(-Q/10)).",
            None,
        ),
        (Tab::Composition, 0) => (
            format!("GC content: {name}"),
            "Reads by GC content, with a normal curve of the same mean and deviation. \
            Peaks away from the curve can mean contamination.",
            None,
        ),
        (Tab::Composition, 1) => (
            if app.relative_pos && !report.relative_positions.is_empty() {
                format!("Relative base composition: {name}")
            } else {
                format!("Base composition: {name}")
            },
            "Share of each base at each position. Lines should be flat after the first \
            few bases of most libraries.",
            None,
        ),
        (Tab::Composition, _) => (
            String::from("N content"),
            "Share of bases called as N at each position.",
            Some(Comparison {
                points: |s| position_points(s, |p| share(p.n_count, p.base_count())),
                is_count: false,
                x_title: "Position (bp)",
                y_title: "N (%)",
            }),
        ),
        (Tab::Length, 0) => (
            String::from("Read length"),
            "Reads by length.",
            Some(Comparison {
                points: |s| Points {
                    x: s.report
                        .length_distribution
                        .iter()
                        .map(|l| l.length as f64)
                        .collect(),
                    y: s.report
                        .length_distribution
                        .iter()
                        .map(|l| l.reads as f64)
                        .collect(),
                },
                is_count: true,
                x_title: "Read length (bp)",
                y_title: "Reads",
            }),
        ),
        (Tab::Length, _) => (
            format!("Quality by read length: {name}"),
            "Effective and arithmetic mean quality stratified by read-length range.",
            None,
        ),
        (Tab::Duplication, _) => (
            format!("Duplication: {name}"),
            "Copies of each sequence across every read. Most reads of a diverse library \
            appear once.",
            None,
        ),
        _ => (
            format!("Adapter content: {name}"),
            "Reads with an adapter starting at or before each position. Rising lines \
            suggest the reads need trimming.",
            None,
        ),
    };
    let mut hints = String::new();
    if matches!(
        (app.tab, section),
        (Tab::Quality, 0) | (Tab::Composition, 1)
    ) && !report.relative_positions.is_empty()
    {
        hints.push_str(if app.relative_pos {
            " · r absolute bp"
        } else {
            " · r relative %"
        });
    }
    if comparison.is_some() {
        hints.push_str(if app.is_overlay() {
            " · x samples"
        } else {
            " · b band or heatmap"
        });
    }
    if !app.is_stacked && sections > 1 {
        hints.push_str(&format!(" · {} of {sections}, n p", section + 1));
    }
    let block = Block::bordered()
        .border_type(BorderType::Rounded)
        .title(Line::from(vec![
            Span::from(format!(" {title}")).bold(),
            Span::from(format!("{hints} ")).dim(),
        ]));
    let inner = block.inner(area);
    frame.render_widget(block, area);
    let mut description = description.to_string();
    if let Some(comparison) = &comparison
        && !app.is_overlay()
    {
        if comparison.is_count {
            description.push_str(" Shares of each sample's reads.");
        }
        description.push_str(match app.spread_view {
            SpreadView::Band => {
                " Every listed sample as percentiles, with the open and outlying samples on top."
            }
            SpreadView::Heatmap => " Every listed sample as a row, in list order.",
        });
    }
    // Long descriptions give way to the chart on short screens.
    let description_rows = if inner.height >= 12 { 2 } else { 0 };
    let [text, chart] =
        Layout::vertical([Constraint::Length(description_rows), Constraint::Min(3)]).areas(inner);
    frame.render_widget(
        Paragraph::new(description).dim().wrap(Wrap { trim: true }),
        inset(text),
    );
    if let Some(comparison) = comparison {
        render_comparison(frame, chart, app, colors, &comparison);
        return;
    }
    match (app.tab, section) {
        (Tab::Quality, 0) => render_base_quality(frame, chart, sample, app.relative_pos, colors),
        (Tab::Quality, _) => render_expected_errors(frame, chart, sample, colors),
        (Tab::Length, _) => render_length_quality(frame, chart, sample, colors),
        (Tab::Composition, 0) => {
            let series = vec![
                Series::named(
                    "Observed",
                    Points::indexed(&report.gc_distribution).pairs(),
                    colors.series_at(0),
                ),
                Series::named(
                    "Theoretical",
                    Points::indexed(&report.gc_theoretical).pairs(),
                    colors.series_at(1),
                )
                .marker(Marker::Dot),
            ];
            line_chart(
                frame,
                chart,
                &series,
                "GC content (%)",
                "Reads",
                None,
                colors,
            );
        }
        (Tab::Composition, _) => {
            render_base_composition(frame, chart, sample, app.relative_pos, colors);
        }
        (Tab::Duplication, _) => render_duplication(frame, chart, sample, colors),
        _ => {
            // Only adapters found in the reads; each keeps its color across samples.
            let series: Vec<Series> = report
                .adapters
                .iter()
                .enumerate()
                .filter(|(_, adapter)| adapter.cumulative_pct.iter().any(|&v| v > 0.0))
                .map(|(i, adapter)| {
                    let points = Points {
                        x: report.positions.iter().map(|p| p.start as f64).collect(),
                        y: adapter.cumulative_pct.clone(),
                    };
                    Series::named(&adapter.name, points.pairs(), colors.series_at(i))
                })
                .collect();
            if series.is_empty() {
                placeholder(frame, chart, "No adapters found.");
                return;
            }
            line_chart(
                frame,
                chart,
                &series,
                "Position (bp)",
                "Reads (%)",
                Some(0.0),
                colors,
            );
        }
    }
}

fn render_base_quality(
    frame: &mut Frame,
    area: Rect,
    sample: &ReadSample,
    relative: bool,
    colors: &Colors,
) {
    if relative && !sample.report.relative_positions.is_empty() {
        let rel = &sample.report.relative_positions;
        let line =
            |value: fn(&segul_core::stats::fastq::report::RelativePositionSummary) -> f64| {
                rel.iter()
                    .map(|p| ((p.start_pct + p.end_pct) as f64 / 2.0, value(p)))
                    .collect::<Vec<_>>()
            };
        let first = rel
            .first()
            .map_or(0.0, |p| (p.start_pct + p.end_pct) as f64 / 2.0);
        let last = rel
            .last()
            .map_or(100.0, |p| (p.start_pct + p.end_pct) as f64 / 2.0);
        let zone = |y: f64| vec![(first, y), (last, y)];
        let series = vec![
            Series::named("Poor <20", zone(POOR_QUALITY), colors.fail).marker(Marker::Dot),
            Series::named("Good 28+", zone(GOOD_QUALITY), colors.pass).marker(Marker::Dot),
            Series::named(
                "10th to 90th",
                line(|p| p.p10_qscore as f64),
                colors.outlier,
            ),
            Series::unnamed(line(|p| p.p90_qscore as f64), colors.outlier),
            Series::named("Middle half", line(|p| p.q1_qscore as f64), colors.axis),
            Series::unnamed(line(|p| p.q3_qscore as f64), colors.axis),
            Series::named(
                "Eff Q",
                line(|p| p.effective_mean_qscore),
                colors.series_at(2),
            ),
            Series::named("Mean", line(|p| p.mean_qscore), colors.series_at(1)),
            Series::named(
                "Median",
                line(|p| p.median_qscore as f64),
                colors.series_at(0),
            ),
        ];
        line_chart(
            frame,
            area,
            &series,
            "Relative position (%)",
            "Quality score",
            Some(0.0),
            colors,
        );
        return;
    }
    let positions = &sample.report.positions;
    let line = |value: fn(&segul_core::stats::fastq::report::PositionSummary) -> f64| {
        positions
            .iter()
            .map(|p| (p.start as f64, value(p)))
            .collect::<Vec<_>>()
    };
    let first = positions.first().map_or(0.0, |p| p.start as f64);
    let last = positions.last().map_or(1.0, |p| p.start as f64);
    let zone = |y: f64| vec![(first, y), (last, y)];
    let mut series = vec![
        Series::named("Poor <20", zone(POOR_QUALITY), colors.fail).marker(Marker::Dot),
        Series::named("Good 28+", zone(GOOD_QUALITY), colors.pass).marker(Marker::Dot),
        Series::named(
            "10th to 90th",
            line(|p| p.p10_qscore as f64),
            colors.outlier,
        ),
        Series::unnamed(line(|p| p.p90_qscore as f64), colors.outlier),
        Series::named("Middle half", line(|p| p.q1_qscore as f64), colors.axis),
        Series::unnamed(line(|p| p.q3_qscore as f64), colors.axis),
        Series::named("Mean", line(|p| p.mean_qscore), colors.series_at(1)),
        Series::named(
            "Median",
            line(|p| p.median_qscore as f64),
            colors.series_at(0),
        ),
    ];
    if positions.iter().any(|p| p.effective_mean_qscore > 0.0) {
        series.push(Series::named(
            "Eff Q",
            line(|p| p.effective_mean_qscore),
            colors.series_at(2),
        ));
    }
    if positions.is_empty() {
        placeholder(frame, area, "No per-position statistics in this report.");
        return;
    }
    line_chart(
        frame,
        area,
        &series,
        "Position (bp)",
        "Quality score",
        Some(0.0),
        colors,
    );
}

fn render_base_composition(
    frame: &mut Frame,
    area: Rect,
    sample: &ReadSample,
    relative: bool,
    colors: &Colors,
) {
    if relative && !sample.report.relative_positions.is_empty() {
        type RelShare = fn(&segul_core::stats::fastq::report::RelativePositionSummary) -> f64;
        let bases: [(&str, RelShare); 5] = [
            ("A", |p| {
                share(
                    p.a_count,
                    p.a_count + p.c_count + p.g_count + p.t_count + p.n_count,
                )
            }),
            ("C", |p| {
                share(
                    p.c_count,
                    p.a_count + p.c_count + p.g_count + p.t_count + p.n_count,
                )
            }),
            ("G", |p| {
                share(
                    p.g_count,
                    p.a_count + p.c_count + p.g_count + p.t_count + p.n_count,
                )
            }),
            ("T", |p| {
                share(
                    p.t_count,
                    p.a_count + p.c_count + p.g_count + p.t_count + p.n_count,
                )
            }),
            ("N", |p| {
                share(
                    p.n_count,
                    p.a_count + p.c_count + p.g_count + p.t_count + p.n_count,
                )
            }),
        ];
        let series: Vec<Series> = bases
            .iter()
            .zip(colors.bases)
            .map(|((base, value), color)| {
                let pairs: Vec<(f64, f64)> = sample
                    .report
                    .relative_positions
                    .iter()
                    .map(|p| ((p.start_pct + p.end_pct) as f64 / 2.0, value(p)))
                    .collect();
                Series::named(base, pairs, color)
            })
            .collect();
        line_chart(
            frame,
            area,
            &series,
            "Relative position (%)",
            "Bases (%)",
            Some(0.0),
            colors,
        );
        return;
    }
    type BaseShare = fn(&segul_core::stats::fastq::report::PositionSummary) -> f64;
    let bases: [(&str, BaseShare); 5] = [
        ("A", |p| share(p.a_count, p.base_count())),
        ("C", |p| share(p.c_count, p.base_count())),
        ("G", |p| share(p.g_count, p.base_count())),
        ("T", |p| share(p.t_count, p.base_count())),
        ("N", |p| share(p.n_count, p.base_count())),
    ];
    let series: Vec<Series> = bases
        .iter()
        .zip(colors.bases)
        .map(|((base, value), color)| {
            Series::named(base, position_points(sample, *value).pairs(), color)
        })
        .collect();
    line_chart(
        frame,
        area,
        &series,
        "Position (bp)",
        "Bases (%)",
        Some(0.0),
        colors,
    );
}

fn render_expected_errors(frame: &mut Frame, area: Rect, sample: &ReadSample, colors: &Colors) {
    let buckets = &sample.report.expected_error_distribution;
    if buckets.is_empty() {
        placeholder(frame, area, "No expected-error statistics in this report.");
        return;
    }
    let m = sample.metrics();
    let [summary, rows] = Layout::vertical([Constraint::Length(1), Constraint::Min(1)]).areas(area);
    frame.render_widget(
        Paragraph::new(format!(
            " Mean EE: {:.2}/read ({:.2}/100 bp) · Error-free est: {:.1}% · EE < 1: {:.1}%",
            m.mean_expected_errors,
            m.expected_errors_per_100bp,
            m.expected_error_free_reads_pct,
            m.reads_ee_lt_1_pct
        )),
        summary,
    );
    let label_width = 10;
    let value_width = 16;
    let bar_width = (area.width as usize).saturating_sub(label_width + value_width + 4);
    let max_pct = buckets
        .iter()
        .map(|b| b.reads_pct)
        .fold(0.0f64, f64::max)
        .max(1.0);
    let lines: Vec<Line> = buckets
        .iter()
        .take(rows.height as usize)
        .map(|b| {
            let cells = bar(b.reads_pct, max_pct, bar_width);
            let pad = " ".repeat(bar_width.saturating_sub(cells.chars().count()) + 1);
            Line::from(vec![
                Span::from(format!(" {:label_width$}", b.label)),
                Span::from(cells).fg(colors.series_at(0)),
                Span::from(pad),
                Span::from(format!(
                    "{:>5.1}% ({})",
                    b.reads_pct,
                    compact(b.reads as f64)
                )),
            ])
        })
        .collect();
    frame.render_widget(Paragraph::new(lines), rows);
}

fn render_length_quality(frame: &mut Frame, area: Rect, sample: &ReadSample, colors: &Colors) {
    let bins = &sample.report.length_quality;
    if bins.is_empty() {
        placeholder(frame, area, "No length-quality statistics in this report.");
        return;
    }
    let [header, rows] = Layout::vertical([Constraint::Length(1), Constraint::Min(1)]).areas(area);
    let label_width = 10;
    let stats_width = 34;
    let bar_width = (area.width as usize).saturating_sub(label_width + stats_width + 4);
    frame.render_widget(
        Paragraph::new(Line::from(vec![
            Span::from(format!(" {:label_width$}", "Length")).bold(),
            Span::from(format!("{:bar_width$} ", "Effective Q")).bold(),
            Span::from(format!(
                "{:>6} {:>6} {:>7} {:>7}",
                "Eff Q", "MeanQ", "Reads", "MeanEE"
            ))
            .bold(),
        ])),
        header,
    );
    let max_q = bins
        .iter()
        .map(|b| b.effective_mean_qscore.max(b.mean_qscore))
        .fold(40.0f64, f64::max);
    let lines: Vec<Line> = bins
        .iter()
        .take(rows.height as usize)
        .map(|b| {
            let cells = bar(b.effective_mean_qscore, max_q, bar_width);
            let pad = " ".repeat(bar_width.saturating_sub(cells.chars().count()) + 1);
            let color = if b.effective_mean_qscore < POOR_QUALITY {
                colors.fail
            } else if b.effective_mean_qscore < GOOD_QUALITY {
                colors.warn
            } else {
                colors.pass
            };
            Line::from(vec![
                Span::from(format!(" {:label_width$}", b.label)),
                Span::from(cells).fg(color),
                Span::from(pad),
                Span::from(format!(
                    "{:>6.1} {:>6.1} {:>7} {:>7.2}",
                    b.effective_mean_qscore,
                    b.mean_qscore,
                    compact(b.reads as f64),
                    b.mean_expected_errors
                )),
            ])
        })
        .collect();
    frame.render_widget(Paragraph::new(lines), rows);
}

fn render_cohort(frame: &mut Frame, area: Rect, app: &App, colors: &Colors) {
    let block = Block::bordered()
        .border_type(BorderType::Rounded)
        .title(Line::from(vec![
            Span::from(" Cohort balance & paired-end integrity ").bold(),
            Span::from(format!("· profile: {} ", app.qc_profile.label())).dim(),
        ]));
    let inner = block.inner(area);
    frame.render_widget(block, area);
    if app.cohort.is_none() && app.pairs.is_empty() {
        placeholder(
            frame,
            inner,
            "No cohort or paired-end summary in this report.",
        );
        return;
    }
    let mut lines: Vec<Line> = Vec::new();
    if let Some(cohort) = &app.cohort {
        lines.push(Line::from(Span::from(" Cohort depth balance").bold()));
        lines.push(Line::from(format!(
            " Samples: {} ({} empty) · Total reads: {} · Total bases: {}",
            thousands(cohort.sample_count as i128),
            thousands(cohort.empty_sample_count as i128),
            thousands(cohort.total_reads as i128),
            compact(cohort.total_bases as f64),
        )));
        lines.push(Line::from(format!(
            " Median reads: {} · P10-P90: {}-{} · P90/P10: {:.2}x · CV: {:.2} · Gini: {:.2}",
            compact(cohort.median_reads),
            compact(cohort.p10_reads),
            compact(cohort.p90_reads),
            cohort.depth_p90_p10_ratio,
            cohort.depth_cv,
            cohort.depth_gini,
        )));
        if !cohort.low_yield_samples.is_empty() {
            lines.push(Line::from(vec![
                Span::from(" Low-yield dropouts (<10% median): ")
                    .fg(colors.warn)
                    .bold(),
                Span::from(cohort.low_yield_samples.join(", ")),
            ]));
        }
        lines.push(Line::from(""));
    }
    if !app.pairs.is_empty() {
        lines.push(Line::from(
            Span::from(format!(
                " Paired-end R1/R2 consistency ({})",
                app.pairs.len()
            ))
            .bold(),
        ));
        lines.push(
            Line::from(format!(
                " {:<24} {:>10} {:>10} {:>7} {:>9} {:>9} {:>8}",
                "Pair", "R1 reads", "R2 reads", "Sync", "Q30 drop", "EffQ drop", "GC diff"
            ))
            .bold(),
        );
        for pair in &app.pairs {
            let status = overall_status(&pair_qc_checks(pair));
            let sync = if pair.in_sync { "OK" } else { "DIFF" };
            lines.push(Line::from(vec![
                Span::from(format!(" {DOT}")).fg(colors.status(status)),
                Span::from(format!(
                    " {:<22} {:>10} {:>10} ",
                    truncate(&pair.pair_name, 22),
                    compact(pair.r1_reads as f64),
                    compact(pair.r2_reads as f64),
                )),
                Span::from(format!("{:>7}", sync)).fg(if pair.in_sync {
                    colors.pass
                } else {
                    colors.fail
                }),
                Span::from(format!(
                    " {:>+8.1}% {:>+9.1} {:>+7.1}%",
                    pair.q30_drop_pct, pair.effective_q_drop, pair.gc_diff_pct
                )),
            ]));
        }
        lines.push(Line::from(""));
    }
    if let Some(cohort) = &app.cohort
        && !cohort.outliers.is_empty()
    {
        lines.push(Line::from(
            Span::from(format!(
                " Cohort-relative outliers ({})",
                cohort.outliers.len()
            ))
            .bold(),
        ));
        lines.push(
            Line::from(format!(
                " {:<24} {:<20} {:>10} {:>10} {:>8}",
                "Sample", "Metric", "Value", "Median", "Mod Z"
            ))
            .bold(),
        );
        for outlier in &cohort.outliers {
            lines.push(Line::from(vec![
                Span::from(format!(
                    " {:<24} {:<20} ",
                    truncate(&outlier.sample, 24),
                    outlier.metric
                )),
                Span::from(format!("{:>10.1}", outlier.value)).fg(colors.warn),
                Span::from(format!(
                    " {:>10.1} {:>+8.1}",
                    outlier.cohort_median, outlier.modified_z
                )),
            ]));
        }
    }
    frame.render_widget(Paragraph::new(lines).wrap(Wrap { trim: false }), inner);
}

fn render_duplication(frame: &mut Frame, area: Rect, sample: &ReadSample, colors: &Colors) {
    let Some(duplication) = &sample.report.duplication else {
        placeholder(
            frame,
            area,
            "Run the read summary in Complete mode to estimate duplication.",
        );
        return;
    };
    let method = if duplication.exact {
        String::from("counted exactly")
    } else {
        format!(
            "estimated from {} sampled sequences",
            thousands(duplication.sketch_size as i128)
        )
    };
    let [text, header, rows] = Layout::vertical([
        Constraint::Length(2),
        Constraint::Length(1),
        Constraint::Min(1),
    ])
    .areas(area);
    frame.render_widget(
        Paragraph::new(format!(
            " {:.1}% of {} reads repeat another read, {method}.",
            duplication.duplicate_pct,
            thousands(duplication.reads as i128)
        )),
        text,
    );
    // Copies, then a bar and share for distinct sequences and for reads.
    let label_width = 7;
    let value_width = 7;
    let bar_width = (area.width as usize).saturating_sub(label_width + 2 * value_width + 4) / 2;
    let (distinct, reads) = (colors.series_at(0), colors.series_at(1));
    frame.render_widget(
        Paragraph::new(Line::from(vec![
            Span::from(format!(" {:label_width$}", "Copies")).bold(),
            Span::from(format!(
                "{:width$}",
                "Distinct sequences (%)",
                width = bar_width + value_width + 1
            ))
            .fg(distinct)
            .bold(),
            Span::from(SEPARATOR).dim(),
            Span::from(" Reads (%)").fg(reads).bold(),
        ])),
        header,
    );
    let column = |pct: f64, color: Color| {
        let cells = bar(pct, 100.0, bar_width);
        let pad = bar_width.saturating_sub(cells.chars().count());
        [
            Span::from(cells).fg(color),
            Span::from(format!("{}{:>value_width$.1} ", " ".repeat(pad), pct)),
        ]
    };
    let lines: Vec<Line> = duplication
        .levels
        .iter()
        .take(rows.height as usize)
        .map(|level| {
            let mut spans = vec![Span::from(format!(" {:label_width$}", level.label))];
            spans.extend(column(level.distinct_pct, distinct));
            spans.push(Span::from(SEPARATOR).dim());
            spans.push(Span::from(" "));
            spans.extend(column(level.reads_pct, reads));
            Line::from(spans)
        })
        .collect();
    frame.render_widget(Paragraph::new(lines), rows);
}

/// The overrepresented sequences of the open sample.
fn render_sequences(frame: &mut Frame, area: Rect, app: &mut App) {
    let sample = app.open_sample();
    let block = Block::bordered()
        .border_type(BorderType::Rounded)
        .title(Line::from(vec![
            Span::from(format!(" Overrepresented sequences: {}", sample.name())).bold(),
            Span::from(" · ← → column · s sort · y copy ").dim(),
        ]));
    let inner = block.inner(area);
    frame.render_widget(block, area);
    let [text, header, body] = Layout::vertical([
        Constraint::Length(2),
        Constraint::Length(1),
        Constraint::Min(1),
    ])
    .areas(inner);
    frame.render_widget(
        Paragraph::new(
            "Read starts above 0.1% of all reads. Reads longer than 75 bp are compared \
            by their first 50 bp.",
        )
        .dim()
        .wrap(Wrap { trim: true }),
        inset(text),
    );
    let sequences = &sample.report.overrepresented;
    if sequences.is_empty() {
        frame.render_widget(
            Paragraph::new(" No overrepresented sequences.").dim(),
            header,
        );
        return;
    }
    let rows = app.sequence_rows();
    app.page_rows = body.height as usize;
    app.sequence_row = app.sequence_row.min(rows.len() - 1);
    keep_visible(
        &mut app.sequence_offset,
        app.sequence_row,
        body.height as usize,
        rows.len(),
    );
    let sequences = &app.open_sample().report.overrepresented;
    let labels = ["Sequence", "Count", "Reads (%)", "Possible source"];
    let count_width = 12;
    let percent_width = 18;
    let rest = (inner.width as usize).saturating_sub(count_width + percent_width + 5);
    let source_width = rest.saturating_sub(rest * 2 / 3).max(8);
    let sequence_width = rest.saturating_sub(source_width).max(10);
    let widths = [sequence_width, count_width, percent_width, source_width];
    let cell = |c: usize, text: String| -> String {
        if c == 1 || c == 2 {
            format!("{:>width$} ", text, width = widths[c] - 1)
        } else {
            format!(
                " {:width$}",
                truncate(&text, widths[c] - 1),
                width = widths[c] - 1
            )
        }
    };
    let mut head = vec![Span::from("")];
    for (c, label) in labels.iter().enumerate() {
        let arrow = match app.sequence_sort {
            Some((s, true)) if s == c => "↓",
            Some((s, false)) if s == c => "↑",
            _ => "",
        };
        let text = cell(c, format!("{arrow}{label}"));
        if c > 0 {
            head.push(Span::from(SEPARATOR).dim());
        }
        head.push(if c == app.sequence_column {
            Span::from(text).fg(accent()).bold()
        } else {
            Span::from(text).bold()
        });
    }
    frame.render_widget(Paragraph::new(Line::from(head)), header);
    let lines: Vec<Line> = rows
        .iter()
        .enumerate()
        .skip(app.sequence_offset)
        .take(body.height as usize)
        .map(|(position, &row)| {
            let seq = &sequences[row];
            let percent = if seq.max_error_pct > 0.0 {
                format!("{:.2} (+{:.3})", seq.percent, seq.max_error_pct)
            } else {
                format!("{:.2}", seq.percent)
            };
            let texts = [
                seq.sequence.clone(),
                thousands(seq.count as i128),
                percent,
                seq.source.clone(),
            ];
            let mut spans = vec![];
            for (c, text) in texts.into_iter().enumerate() {
                if c > 0 {
                    spans.push(Span::from(SEPARATOR).dim());
                }
                spans.push(Span::from(cell(c, text)));
            }
            let line = Line::from(spans);
            if position == app.sequence_row {
                line.add_modifier(Modifier::REVERSED)
            } else {
                line
            }
        })
        .collect();
    frame.render_widget(Paragraph::new(lines), body);
}

/// Points of a chart comparing samples.
struct Comparison {
    points: fn(&ReadSample) -> Points,
    /// Histograms of reads, compared as shares in a band or heatmap since
    /// depths differ. X values a sample lacks count as zero.
    is_count: bool,
    x_title: &'static str,
    y_title: &'static str,
}

/// Chosen samples as lines in small reports, and every listed sample as
/// a band or heatmap in larger ones.
fn render_comparison(
    frame: &mut Frame,
    area: Rect,
    app: &App,
    colors: &Colors,
    comparison: &Comparison,
) {
    if app.is_overlay() {
        if app.compared.is_empty() {
            placeholder(frame, area, "No samples chosen. Press x to choose some.");
            return;
        }
        let series: Vec<Series> = app
            .compared
            .iter()
            .map(|&i| {
                let sample = &app.samples[i];
                Series::named(
                    sample.name(),
                    (comparison.points)(sample).pairs(),
                    colors.series_at(i),
                )
            })
            .collect();
        let min_y = (!comparison.is_count).then_some(0.0);
        line_chart(
            frame,
            area,
            &series,
            comparison.x_title,
            comparison.y_title,
            min_y,
            colors,
        );
        return;
    }
    if app.view.is_empty() {
        placeholder(frame, area, "No samples match the filter.");
        return;
    }
    let points: Vec<Points> = app
        .view
        .iter()
        .map(|&i| {
            let points = (comparison.points)(&app.samples[i]);
            if comparison.is_count {
                points.shares()
            } else {
                points
            }
        })
        .collect();
    let y_title = if comparison.is_count {
        format!("{} (%)", comparison.y_title)
    } else {
        comparison.y_title.to_string()
    };
    let matrix = SampleMatrix::of(&points, comparison.is_count);
    match app.spread_view {
        SpreadView::Band => {
            let spread = SampleSpread::of(&matrix);
            let line = |y: &[f64]| -> Vec<(f64, f64)> {
                spread.x.iter().copied().zip(y.iter().copied()).collect()
            };
            let mut series = vec![
                Series::named("10th to 90th", line(&spread.p10), colors.axis),
                Series::unnamed(line(&spread.p90), colors.axis),
                Series::named("Middle half", line(&spread.q1), colors.series_at(4)),
                Series::unnamed(line(&spread.q3), colors.series_at(4)),
                Series::named("Median", line(&spread.median), colors.series_at(0)),
            ];
            for (n, &row) in spread.outliers.iter().enumerate() {
                let pairs = points[row].pairs();
                series.push(if n == 0 {
                    Series::named("Outliers", pairs, colors.outlier)
                } else {
                    Series::unnamed(pairs, colors.outlier)
                });
            }
            if let Some(row) = app.view.iter().position(|&i| i == app.open) {
                series.push(Series::named(
                    app.open_sample().name(),
                    points[row].pairs(),
                    colors.open,
                ));
            }
            line_chart(
                frame,
                area,
                &series,
                comparison.x_title,
                &y_title,
                Some(0.0),
                colors,
            );
        }
        SpreadView::Heatmap => render_spread_heatmap(frame, area, app, &matrix, comparison),
    }
}

/// One row per listed sample, one cell per x range, colored by value.
fn render_spread_heatmap(
    frame: &mut Frame,
    area: Rect,
    app: &App,
    matrix: &SampleMatrix,
    comparison: &Comparison,
) {
    let [body, legend] = Layout::vertical([Constraint::Min(1), Constraint::Length(1)]).areas(area);
    let name_width = name_width(app.view.iter().map(|&i| app.samples[i].name()));
    let cells = (body.width as usize)
        .saturating_sub(name_width + 3)
        .min(matrix.x.len())
        .max(1);
    let largest = matrix.max_value();
    let truecolor = is_truecolor();
    let open = app.view.iter().position(|&i| i == app.open).unwrap_or(0);
    let rows = body.height as usize;
    // Keeps the open sample in view without state of its own.
    let offset = open.saturating_sub(rows.saturating_sub(1));
    let lines: Vec<Line> = app
        .view
        .iter()
        .zip(&matrix.values)
        .enumerate()
        .skip(offset)
        .take(rows)
        .map(|(position, (&i, values))| {
            let name = Span::from(format!(
                " {:name_width$} ",
                truncate(app.samples[i].name(), name_width)
            ));
            let mut spans = vec![if position == open {
                name.add_modifier(Modifier::REVERSED)
            } else {
                name
            }];
            for cell in 0..cells {
                let from = cell * values.len() / cells;
                let to = ((cell + 1) * values.len() / cells).max(from + 1);
                let present: Vec<f64> = values[from..to].iter().flatten().copied().collect();
                spans.push(if present.is_empty() {
                    Span::from("·").dim()
                } else {
                    let mean = present.iter().sum::<f64>() / present.len() as f64;
                    let color = app.colormap.at(scale(mean, 0.0, largest));
                    Span::from(HEAT).fg(to_color(color, truecolor))
                });
            }
            Line::from(spans)
        })
        .collect();
    frame.render_widget(Paragraph::new(lines), body);
    let mut key = vec![Span::from(" 0 ").dim()];
    key.extend(
        app.colormap
            .stops(8)
            .into_iter()
            .map(|color| Span::from(HEAT).fg(to_color(color, truecolor))),
    );
    let range = match (matrix.x.first(), matrix.x.last()) {
        (Some(first), Some(last)) => format!(
            " · {} from {} to {}",
            comparison.x_title,
            axis_number(*first),
            axis_number(*last)
        ),
        _ => String::new(),
    };
    key.push(Span::from(format!(" {}{range}", axis_number(largest))).dim());
    frame.render_widget(Paragraph::new(Line::from(key)), legend);
}

/// A line of a chart, named in the legend or not.
struct Series {
    name: Option<String>,
    points: Vec<(f64, f64)>,
    color: Color,
    marker: Marker,
}

impl Series {
    fn named(name: &str, points: Vec<(f64, f64)>, color: Color) -> Self {
        Self {
            name: Some(name.to_string()),
            points,
            color,
            marker: Marker::Braille,
        }
    }

    fn unnamed(points: Vec<(f64, f64)>, color: Color) -> Self {
        Self {
            name: None,
            points,
            color,
            marker: Marker::Braille,
        }
    }

    fn marker(mut self, marker: Marker) -> Self {
        self.marker = marker;
        self
    }
}

/// Lines of [series] on shared axes starting at [min_y], or at the lowest
/// value when none, under a legend of the named series.
fn line_chart(
    frame: &mut Frame,
    area: Rect,
    series: &[Series],
    x_title: &str,
    y_title: &str,
    min_y: Option<f64>,
    colors: &Colors,
) {
    let all = || series.iter().flat_map(|s| s.points.iter());
    if all().next().is_none() {
        placeholder(frame, area, "No values to plot.");
        return;
    }
    let [legend, area] = Layout::vertical([Constraint::Length(1), Constraint::Min(2)]).areas(area);
    let mut key = vec![Span::from(" ")];
    for s in series {
        if let Some(name) = &s.name {
            let symbol = if s.marker == Marker::Dot {
                "┅"
            } else {
                "━"
            };
            key.push(Span::from(symbol).fg(s.color));
            key.push(Span::from(format!(" {name}  ")));
        }
    }
    frame.render_widget(Paragraph::new(Line::from(key)), legend);
    let bounds = |values: Vec<f64>| {
        let low = values.iter().copied().fold(f64::INFINITY, f64::min);
        let high = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if high > low {
            [low, high]
        } else {
            [low, low + 1.0]
        }
    };
    let x = bounds(all().map(|p| p.0).collect());
    let mut y = bounds(all().map(|p| p.1).collect());
    if let Some(min) = min_y {
        y[0] = min.min(y[0]);
    }
    let labels = |[low, high]: [f64; 2]| {
        let mut middle = (low + high) / 2.0;
        // Whole-number axes keep a whole middle tick.
        if low.fract() == 0.0 && high.fract() == 0.0 && high - low >= 2.0 {
            middle = middle.round();
        }
        [low, middle, high]
            .into_iter()
            .map(axis_number)
            .collect::<Vec<_>>()
    };
    let datasets: Vec<Dataset> = series
        .iter()
        .map(|s| {
            Dataset::default()
                .marker(s.marker)
                .graph_type(GraphType::Line)
                .style(Style::new().fg(s.color))
                .data(&s.points)
        })
        .collect();
    let axis_style = Style::new().fg(colors.axis);
    let chart = Chart::new(datasets)
        .x_axis(
            Axis::default()
                .title(Line::from(x_title).dim())
                .bounds(x)
                .labels(labels(x))
                .style(axis_style),
        )
        .y_axis(
            Axis::default()
                .title(Line::from(y_title).dim())
                .bounds(y)
                .labels(labels(y))
                .style(axis_style),
        )
        .legend_position(None);
    frame.render_widget(chart, area);
}

fn render_status(frame: &mut Frame, area: Rect, app: &App) {
    let left = if app.filtering {
        Line::from(vec![
            Span::from(" Filter: ").fg(accent()),
            Span::from(format!("{}▏", app.query)),
            Span::from("  Enter keep · Esc clear").dim(),
        ])
    } else if let Some(message) = &app.message {
        Line::from(format!(" {message}")).fg(accent())
    } else {
        let hint = match app.tab {
            Tab::Summary if app.summary_view == SummaryView::Bars => {
                " ↑↓ sample · ← → metric · v view · c checks"
            }
            Tab::Summary => " ↑↓ sample · ← → column · s sort · v view · c checks",
            Tab::Overrepresented => " ↑↓ sequence · ← → column · s sort · y copy",
            _ if app.is_overlay() => " ↑↓ sample · x compare · Tab next tab",
            _ => " ↑↓ sample · b band or heatmap · Tab next tab",
        };
        match (app.tab, app.sequence_rows().get(app.sequence_row)) {
            (Tab::Overrepresented, Some(&row)) => Line::from(format!(
                " {}",
                app.open_sample().report.overrepresented[row].sequence
            )),
            _ => Line::from(hint).dim(),
        }
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
    let width = 66.min(area.width);
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

/// Each check of the open sample with the range of every status, the
/// sample's own range marked, and what the chosen check measures.
fn render_checks_popup(frame: &mut Frame, area: Rect, app: &App, row: usize, colors: &Colors) {
    let sample = app.open_sample();
    let width = 92.min(area.width);
    let height = (sample.checks.len() as u16 + 7).min(area.height);
    let popup = centered(area, width, height);
    let block = Block::bordered()
        .border_type(BorderType::Rounded)
        .title(format!(" Checks for {} ", sample.name()))
        .title_bottom(Line::from(" ↑↓ check · Esc close ").right_aligned());
    let inner = block.inner(popup);
    frame.render_widget(Clear, popup);
    frame.render_widget(block, popup);
    let [table, about] = Layout::vertical([
        Constraint::Length(sample.checks.len() as u16 + 1),
        Constraint::Min(1),
    ])
    .areas(inner);
    let mut lines = vec![
        Line::from(format!(
            " {:18}{:>8}  {:20}{:20}{:20}",
            "Check", "Value", "Pass", "Warn", "Fail"
        ))
        .bold(),
    ];
    for (i, check) in sample.checks.iter().enumerate() {
        let name = Span::from(format!(" {:16}", check.name));
        let mut spans = vec![
            Span::from(format!(" {DOT}")).fg(colors.status(check.status)),
            if i == row {
                name.add_modifier(Modifier::REVERSED)
            } else {
                name
            },
            Span::from(format!("{:>8}  ", check.value_text())),
        ];
        for (status, range) in check.rule.ranges() {
            let text = format!("{:20}", truncate(&range, 19));
            spans.push(if status == check.status {
                Span::from(text).fg(colors.status(status)).bold()
            } else {
                Span::from(text).dim()
            });
        }
        lines.push(Line::from(spans));
    }
    frame.render_widget(Paragraph::new(lines), table);
    if let Some(check) = sample.checks.get(row) {
        let text = vec![
            Line::from(""),
            Line::from(format!(
                " {}: {} ({})",
                check.name,
                check.value_text(),
                check.status.label().to_lowercase()
            ))
            .bold(),
            Line::from(format!(" {}", check.rule.description)),
        ];
        frame.render_widget(Paragraph::new(text).wrap(Wrap { trim: false }), about);
    }
}

/// Every sample with a check box, for the samples overlaid on charts.
fn render_compare_popup(frame: &mut Frame, area: Rect, app: &App, row: usize, colors: &Colors) {
    let width = 50.min(area.width);
    let height = (app.samples.len() as u16 + 2).min(area.height);
    let popup = centered(area, width, height);
    let block = Block::bordered()
        .border_type(BorderType::Rounded)
        .title(format!(
            " Compared: {} of {} samples ",
            app.compared.len(),
            app.samples.len()
        ))
        .title_bottom(Line::from(" Space toggle · a all · n none · Esc ").right_aligned());
    let inner = block.inner(popup);
    let offset = row.saturating_sub((inner.height as usize).saturating_sub(1));
    let lines: Vec<Line> = app
        .samples
        .iter()
        .enumerate()
        .skip(offset)
        .take(inner.height as usize)
        .map(|(i, sample)| {
            let check = if app.compared.contains(&i) {
                "[x]"
            } else {
                "[ ]"
            };
            let line = Line::from(vec![
                Span::from(format!(" {check} ")),
                Span::from(DOT).fg(colors.series_at(i)),
                Span::from(format!(" {}", sample.name())),
            ]);
            if i == row {
                line.add_modifier(Modifier::REVERSED)
            } else {
                line
            }
        })
        .collect();
    frame.render_widget(Clear, popup);
    frame.render_widget(Paragraph::new(lines).block(block), popup);
}

/// [area] less a column on each side, for wrapped text.
fn inset(area: Rect) -> Rect {
    Rect {
        x: area.x + 1,
        width: area.width.saturating_sub(2),
        ..area
    }
}

fn placeholder(frame: &mut Frame, area: Rect, text: &str) {
    frame.render_widget(
        Paragraph::new(text).dim().wrap(Wrap { trim: true }),
        inset(area),
    );
}

/// Points of a per-position statistic of [sample].
fn position_points(
    sample: &ReadSample,
    value: fn(&segul_core::stats::fastq::report::PositionSummary) -> f64,
) -> Points {
    let positions = &sample.report.positions;
    Points {
        x: positions.iter().map(|p| p.start as f64).collect(),
        y: positions.iter().map(value).collect(),
    }
}

/// [count] out of [total], in percent.
fn share(count: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        100.0 * count as f64 / total as f64
    }
}

fn to_f64(counts: &[u64]) -> Vec<f64> {
    counts.iter().map(|&c| c as f64).collect()
}

/// A tick label: counts shortened, whole numbers plain, others with a
/// decimal.
fn axis_number(value: f64) -> String {
    if value.abs() >= 1000.0 {
        compact(value)
    } else if value.fract() == 0.0 {
        format!("{value:.0}")
    } else {
        format!("{value:.1}")
    }
}

fn name_width<'a>(names: impl Iterator<Item = &'a str>) -> usize {
    names
        .map(|n| n.chars().count())
        .max()
        .unwrap_or(0)
        .max("Sample".len())
        .clamp(NAME_WIDTH[0], NAME_WIDTH[1])
}

#[cfg(test)]
mod test {
    use super::*;

    use std::path::Path;

    use crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
    use ratatui::Terminal;
    use ratatui::backend::TestBackend;
    use segul_core::color::Colormap;
    use segul_core::stats::fastq::report::ReadReport;

    use super::super::metric::samples_with_reads;

    fn app() -> App {
        let report = ReadReport::from_path(Path::new(
            "crates/segul-core/tests/files/plot/read-report.json",
        ))
        .expect("Failed reading fixture");
        App::new(
            samples_with_reads(report),
            String::from("read-report.json"),
            Theme::Dark,
            Colormap::Viridis,
        )
    }

    fn screen(app: &mut App, width: u16, height: u16) -> String {
        let mut terminal =
            Terminal::new(TestBackend::new(width, height)).expect("Failed creating terminal");
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
    fn test_render_summary() {
        let mut app = app();
        // Samples without reads are left out.
        assert_eq!(app.samples.len(), 2);
        let text = screen(&mut app, 120, 30);
        assert!(text.contains("read-report.json"));
        assert!(text.contains("2 samples · 4 reads"));
        for tab in Tab::ALL {
            assert!(text.contains(tab.label()));
        }
        assert!(text.contains("Samples (2)"));
        assert!(text.contains("Checks for read_1"));
        assert!(text.contains("Q30 bases:"));
        press(&mut app, KeyCode::Char('v'));
        let text = screen(&mut app, 120, 30);
        assert!(text.contains("Status"));
        assert!(text.contains("Mean Q"));
        press(&mut app, KeyCode::Char('v'));
        let text = screen(&mut app, 120, 30);
        assert!(text.contains("scaled per column"));
        press(&mut app, KeyCode::Char('c'));
        let text = screen(&mut app, 120, 30);
        assert!(text.contains("Checks for read_1"));
        assert!(text.contains("80% or higher"));
        press(&mut app, KeyCode::Esc);
        press(&mut app, KeyCode::Char('?'));
        assert!(screen(&mut app, 120, 30).contains("Keyboard shortcuts"));
    }

    #[test]
    fn test_render_tabs() {
        let mut app = app();
        press(&mut app, KeyCode::Char('2'));
        let text = screen(&mut app, 120, 40);
        assert!(text.contains("Per-base quality: read_1"));
        assert!(text.contains("Per-read quality"));
        // Too short for both charts, so one shows at a time.
        let text = screen(&mut app, 120, 20);
        assert!(text.contains("1 of 2, n p"));
        press(&mut app, KeyCode::Char('n'));
        assert!(screen(&mut app, 120, 20).contains("Per-read quality"));
        press(&mut app, KeyCode::Char('3'));
        let text = screen(&mut app, 120, 50);
        assert!(text.contains("GC content: read_1"));
        assert!(text.contains("Base composition"));
        assert!(text.contains("N content"));
        press(&mut app, KeyCode::Char('5'));
        assert!(screen(&mut app, 120, 30).contains("repeat another read"));
        press(&mut app, KeyCode::Char('6'));
        let text = screen(&mut app, 120, 30);
        assert!(text.contains("Adapter content: read_1"));
        assert!(text.contains("No adapters found."));
        assert!(!text.contains("Nextera"));
        press(&mut app, KeyCode::Char('7'));
        assert!(screen(&mut app, 120, 30).contains("No overrepresented sequences."));
        press(&mut app, KeyCode::Char(']'));
        assert!(screen(&mut app, 120, 30).contains("Overrepresented sequences: read_2"));
    }

    #[test]
    fn test_render_detected_adapters_only() {
        let mut app = app();
        let report = &mut app.samples[0].report;
        let positions = report.positions.len();
        assert!(positions > 0);
        report.adapters[1].cumulative_pct = vec![25.0; positions];
        press(&mut app, KeyCode::Char('6'));
        let text = screen(&mut app, 120, 30);
        assert!(text.contains("Nextera"));
        assert!(!text.contains("Illumina Universal"));
        assert!(!text.contains("PolyG"));
        assert!(!text.contains("No adapters found."));
    }

    #[test]
    fn test_render_many_samples() {
        let names: Vec<String> = (1..=12).map(|i| format!("s{i}")).collect();
        let names: Vec<&str> = names.iter().map(String::as_str).collect();
        let mut samples = super::super::metric::test::samples(&names, &[10; 12], &[90.0; 12]);
        for (i, sample) in samples.iter_mut().enumerate() {
            sample.report.read_quality_distribution = vec![0, i as u64, 10];
        }
        samples[11].report.read_quality_distribution = vec![50, 0, 0];
        let mut app = App::new(
            samples,
            String::from("read-report.json"),
            Theme::Light,
            Colormap::Viridis,
        );
        press(&mut app, KeyCode::Char('2'));
        let text = screen(&mut app, 120, 40);
        assert!(text.contains("b band or heatmap"));
        assert!(text.contains("Outliers"));
        assert!(text.contains("Reads (%)"));
        press(&mut app, KeyCode::Char('b'));
        let text = screen(&mut app, 120, 40);
        assert!(text.contains("as a row"));
        assert!(text.contains(HEAT));
    }

    #[test]
    fn test_render_cohort_and_relative_toggle() {
        let r1 = segul_core::stats::fastq::FastqAnalyzer::new(
            std::path::Path::new("crates/segul-core/tests/files/raw/read_1.fastq"),
            segul_core::helper::types::PhredOffset::Auto,
            true,
            1,
        )
        .analyze()
        .unwrap();
        let r2 = segul_core::stats::fastq::FastqAnalyzer::new(
            std::path::Path::new("crates/segul-core/tests/files/raw/read_2.fastq"),
            segul_core::helper::types::PhredOffset::Auto,
            true,
            1,
        )
        .analyze()
        .unwrap();
        let report = segul_core::stats::fastq::report::ReadReport::from_samples(
            vec![r1, r2],
            segul_core::stats::fastq::qc::QcProfile::Amplicon,
        );
        let qc_profile = report.qc_profile;
        let pairs = report.pairs.clone();
        let cohort = report.cohort.clone();
        let samples = super::super::metric::samples_with_reads(report);
        let mut app = App::new(
            samples,
            String::from("read-report.json"),
            Theme::Dark,
            Colormap::Viridis,
        )
        .with_report_meta(qc_profile, pairs, cohort);
        press(&mut app, KeyCode::Char('2'));
        press(&mut app, KeyCode::Char('r'));
        let text = screen(&mut app, 120, 40);
        assert!(text.contains("Relative-position quality: read_1"));
        assert!(text.contains("Expected errors per read: read_1"));
        press(&mut app, KeyCode::Char('4'));
        let text = screen(&mut app, 120, 40);
        assert!(text.contains("Quality by read length: read_1"));
        press(&mut app, KeyCode::Char('8'));
        let text = screen(&mut app, 120, 30);
        assert!(text.contains("Cohort depth balance"));
        assert!(text.contains("Paired-end R1/R2 consistency"));
        assert!(text.contains("Amplicon / Capture"));
    }

    #[test]
    fn test_axis_number() {
        assert_eq!(axis_number(4_100_000.0), "4.1M");
        assert_eq!(axis_number(28.0), "28");
        assert_eq!(axis_number(0.25), "0.2");
    }
}
