//! State and key handling of the read report viewer.

use std::collections::BTreeSet;

use crossterm::event::{KeyCode, KeyEvent, KeyModifiers, MouseEvent, MouseEventKind};
use segul_core::color::{Colormap, Theme};
use segul_core::stats::fastq::qc::QcProfile;
use segul_core::stats::fastq::report::{CohortSummary, PairedReadSummary};

use super::metric::{ReadMetric, ReadSample, SampleSort, sample_view, step_through};
use super::spread::MAX_OVERLAID;

/// Samples compared on one chart by default, before it gets crowded.
const DEFAULT_COMPARED: usize = 8;

/// Tabs of the report, in the order of the GUI.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Tab {
    Summary,
    Quality,
    Composition,
    Length,
    Duplication,
    Adapters,
    Overrepresented,
    Cohort,
}

impl Tab {
    pub(super) const ALL: [Self; 8] = [
        Self::Summary,
        Self::Quality,
        Self::Composition,
        Self::Length,
        Self::Duplication,
        Self::Adapters,
        Self::Overrepresented,
        Self::Cohort,
    ];

    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Summary => "Summary",
            Self::Quality => "Quality",
            Self::Composition => "Composition",
            Self::Length => "Length",
            Self::Duplication => "Duplication",
            Self::Adapters => "Adapters",
            Self::Overrepresented => "Overrepresented",
            Self::Cohort => "Cohort",
        }
    }

    fn position(self) -> usize {
        Self::ALL.iter().position(|&t| t == self).unwrap_or(0)
    }
}

/// How the summary tab shows the samples.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SummaryView {
    /// One metric as bars.
    Bars,
    Table,
    Heatmap,
}

impl SummaryView {
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Bars => "Bars",
            Self::Table => "Table",
            Self::Heatmap => "Heatmap",
        }
    }
}

/// How charts of many samples show them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SpreadView {
    /// Percentiles across samples, with the open sample and outliers.
    Band,
    /// One row per listed sample.
    Heatmap,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Popup {
    Help,
    /// Ranges of every check of the open sample, with one chosen.
    Checks(usize),
    /// Samples to overlay, with the cursor on one.
    Compare(usize),
}

pub(super) struct App {
    pub(super) samples: Vec<ReadSample>,
    pub(super) file_name: String,
    pub(super) qc_profile: QcProfile,
    pub(super) pairs: Vec<PairedReadSummary>,
    pub(super) cohort: Option<CohortSummary>,
    /// Show relative read position (0-100%) instead of absolute bp on position charts.
    pub(super) relative_pos: bool,
    pub(super) theme: Theme,
    pub(super) colormap: Colormap,
    pub(super) tab: Tab,
    /// The open sample, as a position in [`App::samples`].
    pub(super) open: usize,
    /// Listed samples, filtered and ordered, as positions in
    /// [`App::samples`].
    pub(super) view: Vec<usize>,
    pub(super) query: String,
    /// Whether keys type into the filter.
    pub(super) filtering: bool,
    pub(super) sort: SampleSort,
    pub(super) ascending: bool,
    pub(super) show_list: bool,
    pub(super) summary_view: SummaryView,
    /// Metric of the bars, among [`ReadMetric::available_for`].
    pub(super) metric: usize,
    /// Chosen column of the summary table: 0 is the status, then
    /// [`ReadMetric::TABLE`].
    pub(super) column: usize,
    /// Summary column sorted by, and whether highest first.
    pub(super) table_sort: Option<(usize, bool)>,
    pub(super) spread_view: SpreadView,
    /// Samples overlaid on comparison charts.
    pub(super) compared: BTreeSet<usize>,
    /// The shown chart of a stacked tab when only one fits.
    pub(super) section: usize,
    /// Whether the stacked charts of the open tab fit together, set while
    /// drawing.
    pub(super) is_stacked: bool,
    /// Selected overrepresented sequence, in display order.
    pub(super) sequence_row: usize,
    pub(super) sequence_column: usize,
    pub(super) sequence_sort: Option<(usize, bool)>,
    /// First visible row of the summary, sample list, and sequences.
    pub(super) summary_offset: usize,
    pub(super) list_offset: usize,
    pub(super) sequence_offset: usize,
    /// Rows that fit in the main table, set while drawing.
    pub(super) page_rows: usize,
    pub(super) popup: Option<Popup>,
    pub(super) message: Option<String>,
    /// Text to send to the terminal clipboard.
    pub(super) clipboard: Option<String>,
    pub(super) should_quit: bool,
}

impl App {
    pub(super) fn new(
        samples: Vec<ReadSample>,
        file_name: String,
        theme: Theme,
        colormap: Colormap,
    ) -> Self {
        let sort = SampleSort::Name;
        let view = sample_view(&samples, "", sort, true);
        let compared = (0..samples.len().min(DEFAULT_COMPARED)).collect();
        Self {
            // Opens the first sample of the list.
            open: view.first().copied().unwrap_or(0),
            view,
            samples,
            file_name,
            qc_profile: QcProfile::Auto,
            pairs: Vec::new(),
            cohort: None,
            relative_pos: false,
            theme,
            colormap,
            tab: Tab::Summary,
            query: String::new(),
            filtering: false,
            sort,
            ascending: true,
            show_list: true,
            summary_view: SummaryView::Bars,
            metric: 0,
            column: 1,
            table_sort: None,
            spread_view: SpreadView::Band,
            compared,
            section: 0,
            is_stacked: true,
            sequence_row: 0,
            sequence_column: 0,
            sequence_sort: None,
            summary_offset: 0,
            list_offset: 0,
            sequence_offset: 0,
            page_rows: 1,
            popup: None,
            message: None,
            clipboard: None,
            should_quit: false,
        }
    }

    pub(super) fn with_report_meta(
        mut self,
        qc_profile: QcProfile,
        pairs: Vec<PairedReadSummary>,
        cohort: Option<CohortSummary>,
    ) -> Self {
        self.qc_profile = qc_profile;
        self.pairs = pairs;
        self.cohort = cohort;
        self
    }

    /// Charts stacked in the open tab, adapting to optional fields in the report.
    pub(super) fn sections(&self) -> usize {
        let report = &self.open_sample().report;
        match self.tab {
            Tab::Quality => {
                if report.expected_error_distribution.is_empty() {
                    2
                } else {
                    3
                }
            }
            Tab::Composition => 3,
            Tab::Length => {
                if report.length_quality.is_empty() {
                    1
                } else {
                    2
                }
            }
            _ => 1,
        }
    }

    pub(super) fn open_sample(&self) -> &ReadSample {
        &self.samples[self.open]
    }

    pub(super) fn metrics(&self) -> Vec<ReadMetric> {
        ReadMetric::available_for(&self.samples)
    }

    /// The metric of the bars.
    pub(super) fn metric(&self) -> Option<ReadMetric> {
        let metrics = self.metrics();
        metrics
            .get(self.metric.min(metrics.len().saturating_sub(1)))
            .copied()
    }

    /// Whether comparison charts draw one line per chosen sample.
    pub(super) fn is_overlay(&self) -> bool {
        self.samples.len() <= MAX_OVERLAID
    }

    /// Rows of the summary table: the listed samples, sorted by the chosen
    /// column when asked.
    pub(super) fn summary_rows(&self) -> Vec<usize> {
        let mut rows = self.view.clone();
        if let Some((column, descending)) = self.table_sort {
            let key = |i: usize| self.summary_value(i, column);
            // Samples without a value go last either way.
            rows.sort_by(|&a, &b| match (key(a), key(b)) {
                (Some(x), Some(y)) if descending => y.total_cmp(&x),
                (Some(x), Some(y)) => x.total_cmp(&y),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            });
        }
        rows
    }

    /// Value of sample [i] in summary column [column], with the status as
    /// its rank.
    pub(super) fn summary_value(&self, i: usize, column: usize) -> Option<f64> {
        let sample = &self.samples[i];
        match column {
            0 => Some(sample.status as usize as f64),
            c => ReadMetric::TABLE[c - 1].value_of(sample.metrics()),
        }
    }

    /// Positions of overrepresented sequences of the open sample in
    /// display order.
    pub(super) fn sequence_rows(&self) -> Vec<usize> {
        let sequences = &self.open_sample().report.overrepresented;
        let mut rows: Vec<usize> = (0..sequences.len()).collect();
        if let Some((column, descending)) = self.sequence_sort {
            rows.sort_by(|&a, &b| {
                let (a, b) = (&sequences[a], &sequences[b]);
                let ordered = match column {
                    0 => a.sequence.cmp(&b.sequence),
                    1 => a.count.cmp(&b.count),
                    2 => a.percent.total_cmp(&b.percent),
                    _ => a.source.cmp(&b.source),
                };
                if descending {
                    ordered.reverse()
                } else {
                    ordered
                }
            });
        }
        rows
    }

    pub(super) fn handle_key(&mut self, key: KeyEvent) {
        self.message = None;
        if key.code == KeyCode::Char('c') && key.modifiers.contains(KeyModifiers::CONTROL) {
            self.should_quit = true;
            return;
        }
        if let Some(popup) = self.popup {
            self.handle_popup(popup, key);
            return;
        }
        if self.filtering {
            self.handle_filter(key);
            return;
        }
        let alt = key.modifiers.contains(KeyModifiers::ALT);
        let page = self.page_rows.saturating_sub(1).max(1) as isize;
        match key.code {
            KeyCode::Esc | KeyCode::Char('q') => self.should_quit = true,
            KeyCode::F(1) | KeyCode::Char('?') => self.popup = Some(Popup::Help),
            KeyCode::Tab => self.go_to_tab(self.tab.position() + 1),
            KeyCode::BackTab => self.go_to_tab(self.tab.position() + Tab::ALL.len() - 1),
            KeyCode::Char(c @ '1'..='8') => self.go_to_tab(c as usize - '1' as usize),
            KeyCode::Up if alt => self.step_sample(-1),
            KeyCode::Down if alt => self.step_sample(1),
            KeyCode::Char('[') => self.step_sample(-1),
            KeyCode::Char(']') => self.step_sample(1),
            KeyCode::Up => self.move_row(-1),
            KeyCode::Down => self.move_row(1),
            KeyCode::PageUp => self.move_row(-page),
            KeyCode::PageDown => self.move_row(page),
            KeyCode::Home => self.move_row(isize::MIN / 2),
            KeyCode::End => self.move_row(isize::MAX / 2),
            KeyCode::Left => self.move_column(-1),
            KeyCode::Right => self.move_column(1),
            KeyCode::Char('v') => self.cycle_summary_view(),
            KeyCode::Char('s') => self.cycle_sort(),
            KeyCode::Char('c') => self.popup = Some(Popup::Checks(0)),
            KeyCode::Char('x') => self.open_compare(),
            KeyCode::Char('b') => self.toggle_spread(),
            KeyCode::Char('r') => self.toggle_relative_pos(),
            KeyCode::Char('n') => self.move_section(1),
            KeyCode::Char('p') => self.move_section(-1),
            KeyCode::Char('l') => self.show_list = !self.show_list,
            KeyCode::Char('/') => self.filtering = true,
            KeyCode::Char('o') => {
                self.sort = self.sort.next();
                self.refresh_view();
                self.message = Some(format!("Samples by {}", self.sort.describe(self.ascending)));
            }
            KeyCode::Char('O') => {
                self.ascending = !self.ascending;
                self.refresh_view();
                self.message = Some(format!("Samples by {}", self.sort.describe(self.ascending)));
            }
            KeyCode::Char('y') => self.copy_sequence(),
            _ => (),
        }
    }

    pub(super) fn handle_mouse(&mut self, mouse: MouseEvent) {
        if self.popup.is_some() {
            return;
        }
        match mouse.kind {
            MouseEventKind::ScrollDown => self.move_row(3),
            MouseEventKind::ScrollUp => self.move_row(-3),
            _ => (),
        }
    }

    fn handle_popup(&mut self, popup: Popup, key: KeyEvent) {
        let close = matches!(key.code, KeyCode::Esc | KeyCode::Enter | KeyCode::Char('q'));
        match popup {
            Popup::Help => {
                if close || matches!(key.code, KeyCode::F(1) | KeyCode::Char('?')) {
                    self.popup = None;
                }
            }
            Popup::Checks(row) => {
                let last = self.open_sample().checks.len().saturating_sub(1);
                match key.code {
                    KeyCode::Up => self.popup = Some(Popup::Checks(row.saturating_sub(1))),
                    KeyCode::Down => self.popup = Some(Popup::Checks((row + 1).min(last))),
                    KeyCode::Char('c') => self.popup = None,
                    _ if close => self.popup = None,
                    _ => (),
                }
            }
            Popup::Compare(row) => {
                let last = self.samples.len().saturating_sub(1);
                match key.code {
                    KeyCode::Up => self.popup = Some(Popup::Compare(row.saturating_sub(1))),
                    KeyCode::Down => self.popup = Some(Popup::Compare((row + 1).min(last))),
                    KeyCode::Char(' ') => {
                        if !self.compared.remove(&row) {
                            self.compared.insert(row);
                        }
                    }
                    KeyCode::Char('a') => self.compared = (0..self.samples.len()).collect(),
                    KeyCode::Char('n') => self.compared.clear(),
                    KeyCode::Char('x') => self.popup = None,
                    _ if close => self.popup = None,
                    _ => (),
                }
            }
        }
    }

    /// Typing edits the filter. Enter keeps it and Esc clears it.
    fn handle_filter(&mut self, key: KeyEvent) {
        match key.code {
            KeyCode::Enter => self.filtering = false,
            KeyCode::Esc => {
                self.filtering = false;
                self.query.clear();
            }
            KeyCode::Backspace => {
                self.query.pop();
            }
            KeyCode::Char(c) => self.query.push(c),
            _ => return,
        }
        self.refresh_view();
    }

    fn refresh_view(&mut self) {
        self.view = sample_view(&self.samples, &self.query, self.sort, self.ascending);
    }

    fn go_to_tab(&mut self, position: usize) {
        self.tab = Tab::ALL[position % Tab::ALL.len()];
        self.section = 0;
    }

    /// Opens sample [index], starting its sequences from the top.
    pub(super) fn select(&mut self, index: usize) {
        if index != self.open {
            self.open = index;
            self.sequence_row = 0;
            self.sequence_offset = 0;
        }
    }

    fn step_sample(&mut self, delta: isize) {
        if self.view.is_empty() {
            self.message = Some(String::from("No samples match the filter"));
            return;
        }
        let next = step_through(&self.view, self.open, delta);
        self.select(next);
    }

    /// Moves through the rows of the open table. On charts, steps through
    /// the samples.
    fn move_row(&mut self, delta: isize) {
        match self.tab {
            Tab::Summary => {
                let rows = self.summary_rows();
                if rows.is_empty() {
                    return;
                }
                let position = rows.iter().position(|&i| i == self.open);
                let next = match position {
                    Some(p) => clamp_step(p, delta, rows.len()),
                    None if delta > 0 => 0,
                    None => rows.len() - 1,
                };
                self.select(rows[next]);
            }
            Tab::Overrepresented => {
                let rows = self.open_sample().report.overrepresented.len();
                if rows > 0 {
                    self.sequence_row = clamp_step(self.sequence_row, delta, rows);
                }
            }
            _ => self.step_sample(delta.signum()),
        }
    }

    fn move_column(&mut self, delta: isize) {
        match (self.tab, self.summary_view) {
            (Tab::Summary, SummaryView::Bars) => {
                let metrics = self.metrics().len();
                if metrics > 0 {
                    self.metric = clamp_step(self.metric.min(metrics - 1), delta, metrics);
                }
            }
            (Tab::Summary, _) => {
                self.column = clamp_step(self.column, delta, ReadMetric::TABLE.len() + 1);
            }
            (Tab::Overrepresented, _) => {
                self.sequence_column = clamp_step(self.sequence_column, delta, 4);
            }
            _ => self.step_sample(delta.signum()),
        }
    }

    fn cycle_summary_view(&mut self) {
        if self.tab != Tab::Summary {
            self.message = Some(String::from("Views switch on the Summary tab"));
            return;
        }
        self.summary_view = match self.summary_view {
            SummaryView::Bars => SummaryView::Table,
            SummaryView::Table => SummaryView::Heatmap,
            SummaryView::Heatmap => SummaryView::Bars,
        };
    }

    /// Sorts by the chosen column, highest first, then lowest first, then
    /// back to the list order.
    fn cycle_sort(&mut self) {
        let next = |sort: Option<(usize, bool)>, column: usize| match sort {
            Some((c, true)) if c == column => Some((c, false)),
            Some((c, false)) if c == column => None,
            _ => Some((column, true)),
        };
        match self.tab {
            Tab::Summary if self.summary_view != SummaryView::Bars => {
                self.table_sort = next(self.table_sort, self.column);
                self.message = Some(match self.table_sort {
                    Some((c, true)) => format!("Sorted by {}, highest first", column_label(c)),
                    Some((c, false)) => format!("Sorted by {}, lowest first", column_label(c)),
                    None => String::from("List order"),
                });
            }
            Tab::Overrepresented => {
                self.sequence_sort = next(self.sequence_sort, self.sequence_column);
                self.sequence_row = 0;
            }
            _ => self.message = Some(String::from("Sort the sample list with o")),
        }
    }

    fn open_compare(&mut self) {
        if !self.is_overlay() {
            self.message = Some(format!(
                "Reports of more than {MAX_OVERLAID} samples show every sample"
            ));
            return;
        }
        self.popup = Some(Popup::Compare(self.open));
    }

    fn toggle_spread(&mut self) {
        if self.is_overlay() {
            self.message = Some(String::from("Choose samples to compare with x"));
            return;
        }
        self.spread_view = match self.spread_view {
            SpreadView::Band => SpreadView::Heatmap,
            SpreadView::Heatmap => SpreadView::Band,
        };
    }

    fn toggle_relative_pos(&mut self) {
        if self.open_sample().report.relative_positions.is_empty() {
            self.message = Some(String::from("No relative-position data in this report"));
            return;
        }
        self.relative_pos = !self.relative_pos;
        self.message = Some(if self.relative_pos {
            String::from("Showing relative read position (0-100%)")
        } else {
            String::from("Showing absolute read position (bp)")
        });
    }

    fn move_section(&mut self, delta: isize) {
        let sections = self.sections();
        if sections < 2 || self.is_stacked {
            self.message = Some(String::from("Every chart of this tab is shown"));
            return;
        }
        self.section = (self.section as isize + delta).rem_euclid(sections as isize) as usize;
    }

    fn copy_sequence(&mut self) {
        if self.tab != Tab::Overrepresented {
            self.message = Some(String::from("Copy sequences on the Overrepresented tab"));
            return;
        }
        let rows = self.sequence_rows();
        let Some(&row) = rows.get(self.sequence_row) else {
            self.message = Some(String::from("No sequence to copy"));
            return;
        };
        self.clipboard = Some(
            self.open_sample().report.overrepresented[row]
                .sequence
                .clone(),
        );
        self.message = Some(String::from("Sequence copied"));
    }
}

/// Heading of summary column [column].
pub(super) fn column_label(column: usize) -> &'static str {
    match column {
        0 => "Status",
        c => ReadMetric::TABLE[c - 1].column_label(),
    }
}

/// [position] moved by [delta] within `0..len`.
fn clamp_step(position: usize, delta: isize, len: usize) -> usize {
    (position as isize)
        .saturating_add(delta)
        .clamp(0, len as isize - 1) as usize
}

/// Moves [offset] so [position] is among [rows] visible rows of [len].
pub(super) fn keep_visible(offset: &mut usize, position: usize, rows: usize, len: usize) {
    let rows = rows.max(1);
    if position < *offset {
        *offset = position;
    } else if position >= *offset + rows {
        *offset = position + 1 - rows;
    }
    *offset = (*offset).min(len.saturating_sub(rows));
}

#[cfg(test)]
mod test {
    use super::*;

    use super::super::metric::test::samples;

    fn key(code: KeyCode) -> KeyEvent {
        KeyEvent::new(code, KeyModifiers::NONE)
    }

    fn sample_app(count: usize) -> App {
        let names: Vec<String> = (1..=count).map(|i| format!("s{i}")).collect();
        let names: Vec<&str> = names.iter().map(String::as_str).collect();
        let reads: Vec<u64> = (1..=count as u64).collect();
        let q30: Vec<f64> = (0..count).map(|i| 60.0 + 5.0 * i as f64).collect();
        App::new(
            samples(&names, &reads, &q30),
            String::from("read-report.json"),
            Theme::Dark,
            Colormap::Viridis,
        )
    }

    #[test]
    fn test_step_samples() {
        let mut app = sample_app(3);
        assert_eq!(app.open, 0);
        app.handle_key(KeyEvent::new(KeyCode::Down, KeyModifiers::ALT));
        assert_eq!(app.open, 1);
        app.handle_key(key(KeyCode::Char(']')));
        app.handle_key(key(KeyCode::Char(']')));
        assert_eq!(app.open, 0);
        app.handle_key(key(KeyCode::Char('[')));
        assert_eq!(app.open, 2);
    }

    #[test]
    fn test_filter_and_sort() {
        let mut app = sample_app(12);
        app.handle_key(key(KeyCode::Char('/')));
        for c in "s1".chars() {
            app.handle_key(key(KeyCode::Char(c)));
        }
        app.handle_key(key(KeyCode::Enter));
        assert!(!app.filtering);
        // s1, s10, s11, s12 in natural order.
        assert_eq!(app.view, vec![0, 9, 10, 11]);
        app.handle_key(key(KeyCode::Char('O')));
        assert_eq!(app.view, vec![11, 10, 9, 0]);
        app.handle_key(key(KeyCode::Char('o')));
        assert_eq!(app.sort, SampleSort::Status);
        app.handle_key(key(KeyCode::Char('/')));
        app.handle_key(key(KeyCode::Esc));
        assert_eq!(app.view.len(), 12);
    }

    #[test]
    fn test_summary_keys() {
        let mut app = sample_app(3);
        app.page_rows = 10;
        app.handle_key(key(KeyCode::Char('v')));
        assert_eq!(app.summary_view, SummaryView::Table);
        app.handle_key(key(KeyCode::Down));
        assert_eq!(app.open, 1);
        // Reads, highest first.
        app.handle_key(key(KeyCode::Char('s')));
        assert_eq!(app.summary_rows(), vec![2, 1, 0]);
        app.handle_key(key(KeyCode::Down));
        assert_eq!(app.open, 0);
        app.handle_key(key(KeyCode::Char('v')));
        app.handle_key(key(KeyCode::Char('v')));
        assert_eq!(app.summary_view, SummaryView::Bars);
        app.handle_key(key(KeyCode::Right));
        assert_eq!(app.metric(), Some(ReadMetric::Bases));
        app.handle_key(key(KeyCode::Char('c')));
        assert_eq!(app.popup, Some(Popup::Checks(0)));
        app.handle_key(key(KeyCode::Down));
        assert_eq!(app.popup, Some(Popup::Checks(1)));
        app.handle_key(key(KeyCode::Esc));
        assert!(app.popup.is_none() && !app.should_quit);
    }

    #[test]
    fn test_compare_or_spread() {
        let mut app = sample_app(10);
        assert!(app.is_overlay());
        assert_eq!(app.compared.len(), 8);
        app.handle_key(key(KeyCode::Char('x')));
        app.handle_key(key(KeyCode::Char('n')));
        assert!(app.compared.is_empty());
        app.handle_key(key(KeyCode::Char(' ')));
        assert!(app.compared.contains(&0));
        app.handle_key(key(KeyCode::Char('x')));
        app.handle_key(key(KeyCode::Char('b')));
        assert!(app.message.is_some());

        let mut app = sample_app(11);
        assert!(!app.is_overlay());
        app.handle_key(key(KeyCode::Char('b')));
        assert_eq!(app.spread_view, SpreadView::Heatmap);
        app.handle_key(key(KeyCode::Char('x')));
        assert!(app.popup.is_none());
    }

    #[test]
    fn test_keep_visible() {
        let mut offset = 0;
        keep_visible(&mut offset, 5, 3, 10);
        assert_eq!(offset, 3);
        keep_visible(&mut offset, 1, 3, 10);
        assert_eq!(offset, 1);
        keep_visible(&mut offset, 1, 20, 10);
        assert_eq!(offset, 0);
    }
}
