//! State and key handling of the summary viewer.

use crossterm::event::{KeyCode, KeyEvent, KeyModifiers, MouseEvent, MouseEventKind};
use segul_core::color::Colormap;

use super::data::{Column, SummaryData, SummaryMatrix, SummaryTable};

/// How a table section shows its items.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum Mode {
    Table,
    Plot,
}

/// What a plot draws for the chosen column.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum PlotKind {
    /// One bar per item, in table order.
    Bars,
    /// Items counted in ranges of the column's values.
    Histogram,
    /// Every column, colored by value, one row per item.
    Heatmap,
}

/// Selection and order of one table.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct TableState {
    /// Items in display order, as positions in the table.
    pub(super) order: Vec<usize>,
    /// The selected item, as a position in [`TableState::order`].
    pub(super) row: usize,
    /// First visible item, as a position in [`TableState::order`].
    pub(super) offset: usize,
    /// The chosen column, among the value columns.
    pub(super) column: usize,
    /// Column sorted by, and whether highest first.
    pub(super) sort: Option<(usize, bool)>,
}

impl TableState {
    fn new(table: &SummaryTable) -> Self {
        Self {
            order: (0..table.names.len()).collect(),
            row: 0,
            offset: 0,
            column: 0,
            sort: None,
        }
    }
}

/// Cursor and scroll of one matrix.
#[derive(Debug, Clone, Default, PartialEq)]
pub(super) struct MatrixState {
    pub(super) row: usize,
    pub(super) column: usize,
    /// First visible row and column.
    pub(super) row_offset: usize,
    pub(super) column_offset: usize,
}

/// What the open section shows.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum Section {
    Overview,
    Table(usize),
    Matrix(usize),
}

pub(super) struct App {
    pub(super) data: SummaryData,
    pub(super) file_name: String,
    /// 0 is the overview; tables follow, then matrices.
    pub(super) section: usize,
    pub(super) mode: Mode,
    pub(super) plot: PlotKind,
    /// Colors of the heatmap.
    pub(super) colormap: Colormap,
    pub(super) tables: Vec<TableState>,
    pub(super) matrices: Vec<MatrixState>,
    pub(super) show_help: bool,
    pub(super) should_quit: bool,
    pub(super) message: Option<String>,
    /// Items that fit in the table or plot, set while drawing.
    pub(super) page_rows: usize,
    /// Matrix columns that fit, set while drawing.
    pub(super) page_columns: usize,
}

impl App {
    pub(super) fn new(data: SummaryData, file_name: String, colormap: Colormap) -> Self {
        let tables = data.tables.iter().map(TableState::new).collect();
        let matrices = vec![MatrixState::default(); data.matrices.len()];
        Self {
            data,
            file_name,
            // Opens on the first table, the part most often wanted.
            section: 1,
            mode: Mode::Table,
            plot: PlotKind::Bars,
            colormap,
            tables,
            matrices,
            show_help: false,
            should_quit: false,
            message: None,
            page_rows: 1,
            page_columns: 1,
        }
    }

    pub(super) fn sections(&self) -> usize {
        self.data.tables.len() + self.data.matrices.len() + 1
    }

    pub(super) fn open_section(&self) -> Section {
        let tables = self.data.tables.len();
        match self.section {
            0 => Section::Overview,
            s if s <= tables => Section::Table(s - 1),
            s => Section::Matrix(s - 1 - tables),
        }
    }

    /// The open table and its state, or none on other sections.
    pub(super) fn table(&self) -> Option<(&SummaryTable, &TableState)> {
        let Section::Table(i) = self.open_section() else {
            return None;
        };
        Some((self.data.tables.get(i)?, self.tables.get(i)?))
    }

    /// The open matrix and its state, or none on other sections.
    pub(super) fn matrix(&self) -> Option<(&SummaryMatrix, &MatrixState)> {
        let Section::Matrix(i) = self.open_section() else {
            return None;
        };
        Some((self.data.matrices.get(i)?, self.matrices.get(i)?))
    }

    /// The chosen column of the open table.
    pub(super) fn column(&self) -> Option<&Column> {
        let (table, state) = self.table()?;
        table.columns.get(state.column)
    }

    /// The selected item of the open table, as a position in the table.
    pub(super) fn selected(&self) -> Option<usize> {
        let (_, state) = self.table()?;
        state.order.get(state.row).copied()
    }

    pub(super) fn handle_key(&mut self, key: KeyEvent) {
        self.message = None;
        let ctrl = key.modifiers.contains(KeyModifiers::CONTROL);
        if self.show_help {
            if matches!(
                key.code,
                KeyCode::Esc | KeyCode::F(1) | KeyCode::Enter | KeyCode::Char('?' | 'q')
            ) {
                self.show_help = false;
            }
            return;
        }
        let page = self.page_rows.saturating_sub(1).max(1) as isize;
        match key.code {
            KeyCode::Esc | KeyCode::Char('q') => self.should_quit = true,
            KeyCode::Char('c') if ctrl => self.should_quit = true,
            KeyCode::F(1) | KeyCode::Char('?') => self.show_help = true,
            KeyCode::Tab => self.section = (self.section + 1) % self.sections(),
            KeyCode::BackTab => {
                self.section = (self.section + self.sections() - 1) % self.sections()
            }
            KeyCode::Char(c @ '1'..='9') => {
                let section = c as usize - '1' as usize;
                if section < self.sections() {
                    self.section = section;
                }
            }
            KeyCode::Char('v') | KeyCode::F(2) => self.toggle_mode(),
            KeyCode::Char('h') => self.toggle_plot(),
            KeyCode::Char('s') => self.cycle_sort(),
            KeyCode::Left => self.move_column(-1),
            KeyCode::Right => self.move_column(1),
            KeyCode::Up => self.move_row(-1),
            KeyCode::Down => self.move_row(1),
            KeyCode::PageUp => self.move_row(-page),
            KeyCode::PageDown => self.move_row(page),
            KeyCode::Home => self.move_row(isize::MIN / 2),
            KeyCode::End => self.move_row(isize::MAX / 2),
            _ => (),
        }
    }

    pub(super) fn handle_mouse(&mut self, mouse: MouseEvent) {
        match mouse.kind {
            MouseEventKind::ScrollDown => self.move_row(3),
            MouseEventKind::ScrollUp => self.move_row(-3),
            _ => (),
        }
    }

    fn toggle_mode(&mut self) {
        if self.matrix().is_some() {
            self.message = Some(String::from("A matrix has a single view"));
            return;
        }
        if self.table().is_none() {
            self.message = Some(String::from("Open a table with Tab to switch views"));
            return;
        }
        self.mode = match self.mode {
            Mode::Table => Mode::Plot,
            Mode::Plot => Mode::Table,
        };
    }

    fn toggle_plot(&mut self) {
        if self.table().is_none() || self.mode != Mode::Plot {
            self.message = Some(String::from("Press v for the plot first"));
            return;
        }
        self.plot = match self.plot {
            PlotKind::Bars => PlotKind::Histogram,
            PlotKind::Histogram => PlotKind::Heatmap,
            PlotKind::Heatmap => PlotKind::Bars,
        };
    }

    fn move_column(&mut self, delta: isize) {
        if let Section::Matrix(i) = self.open_section() {
            let last = self.data.matrices[i].names.len().saturating_sub(1) as isize;
            let state = &mut self.matrices[i];
            state.column = (state.column as isize).saturating_add(delta).clamp(0, last) as usize;
            self.keep_cell_visible();
            return;
        }
        let Section::Table(i) = self.open_section() else {
            return;
        };
        let columns = self.data.tables[i].columns.len();
        if columns == 0 {
            return;
        }
        let state = &mut self.tables[i];
        state.column = (state.column as isize + delta).clamp(0, columns as isize - 1) as usize;
    }

    fn move_row(&mut self, delta: isize) {
        if let Section::Matrix(i) = self.open_section() {
            let last = self.data.matrices[i].names.len().saturating_sub(1) as isize;
            let state = &mut self.matrices[i];
            state.row = (state.row as isize).saturating_add(delta).clamp(0, last) as usize;
            self.keep_cell_visible();
            return;
        }
        let Section::Table(i) = self.open_section() else {
            return;
        };
        let state = &mut self.tables[i];
        let last = state.order.len().saturating_sub(1) as isize;
        state.row = (state.row as isize).saturating_add(delta).clamp(0, last) as usize;
        self.keep_row_visible();
    }

    /// Sorts by the chosen column, highest first, then lowest first, then
    /// back to the report order. The selected item stays selected.
    fn cycle_sort(&mut self) {
        let Section::Table(i) = self.open_section() else {
            return;
        };
        let table = &self.data.tables[i];
        let state = &mut self.tables[i];
        let selected = state.order.get(state.row).copied();
        state.sort = match state.sort {
            Some((column, true)) if column == state.column => Some((column, false)),
            Some((column, false)) if column == state.column => None,
            _ => Some((state.column, true)),
        };
        state.order = (0..table.names.len()).collect();
        if let Some((column, descending)) = state.sort {
            let values = &table.columns[column].values;
            // Items without a value go last either way.
            state.order.sort_by(|&a, &b| match (values[a], values[b]) {
                (Some(x), Some(y)) if descending => y.total_cmp(&x),
                (Some(x), Some(y)) => x.total_cmp(&y),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            });
        }
        if let Some(selected) = selected {
            state.row = state.order.iter().position(|&x| x == selected).unwrap_or(0);
        }
        self.message = Some(match state.sort {
            Some((column, true)) => {
                format!("Sorted by {}, highest first", table.columns[column].label)
            }
            Some((column, false)) => {
                format!("Sorted by {}, lowest first", table.columns[column].label)
            }
            None => String::from("Report order"),
        });
        self.keep_row_visible();
    }

    /// Scrolls the open table so the selected item is in view.
    pub(super) fn keep_row_visible(&mut self) {
        let Section::Table(i) = self.open_section() else {
            return;
        };
        let rows = self.page_rows.max(1);
        let state = &mut self.tables[i];
        if state.row < state.offset {
            state.offset = state.row;
        } else if state.row >= state.offset + rows {
            state.offset = state.row + 1 - rows;
        }
        let last_page = state.order.len().saturating_sub(rows);
        state.offset = state.offset.min(last_page);
    }

    /// Scrolls the open matrix so the cursor cell is in view.
    pub(super) fn keep_cell_visible(&mut self) {
        let Section::Matrix(i) = self.open_section() else {
            return;
        };
        let items = self.data.matrices[i].names.len();
        let state = &mut self.matrices[i];
        let scroll = |cursor: usize, offset: usize, page: usize| {
            let page = page.max(1);
            let offset = if cursor < offset {
                cursor
            } else if cursor >= offset + page {
                cursor + 1 - page
            } else {
                offset
            };
            offset.min(items.saturating_sub(page))
        };
        state.row_offset = scroll(state.row, state.row_offset, self.page_rows);
        state.column_offset = scroll(state.column, state.column_offset, self.page_columns);
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use super::super::data::Format;

    fn data() -> SummaryData {
        SummaryData {
            kind: "Alignment summary",
            headline: String::from("3 loci"),
            overview: vec![(String::from("Loci"), String::from("3"))],
            charts: Vec::new(),
            tables: vec![SummaryTable {
                name: "Loci",
                item: "Locus",
                names: vec![String::from("a"), String::from("b"), String::from("c")],
                columns: vec![
                    Column {
                        label: "Sites",
                        values: vec![Some(20.0), Some(30.0), Some(10.0)],
                        format: Format::Count,
                    },
                    Column {
                        label: "GC %",
                        values: vec![Some(40.0), None, Some(45.0)],
                        format: Format::Decimal,
                    },
                ],
            }],
            matrices: vec![SummaryMatrix {
                name: "Matrix",
                item: "samples",
                names: vec![String::from("a"), String::from("b"), String::from("c")],
                values: vec![
                    vec![Some(0.0), Some(0.1), None],
                    vec![Some(0.1), Some(0.0), Some(0.2)],
                    vec![None, Some(0.2), Some(0.0)],
                ],
                format: Format::Distance,
            }],
        }
    }

    fn key(code: KeyCode) -> KeyEvent {
        KeyEvent::new(code, KeyModifiers::NONE)
    }

    #[test]
    fn test_sections_and_modes() {
        let mut app = App::new(data(), String::from("report.json"), Colormap::Viridis);
        assert_eq!(app.section, 1);
        app.handle_key(key(KeyCode::Tab));
        assert_eq!(app.section, 2);
        app.handle_key(key(KeyCode::Tab));
        assert_eq!(app.section, 0);
        // The overview has no plot of its own.
        app.handle_key(key(KeyCode::Char('v')));
        assert_eq!(app.mode, Mode::Table);
        assert!(app.message.is_some());
        app.handle_key(key(KeyCode::Char('2')));
        assert_eq!(app.section, 1);
        app.handle_key(key(KeyCode::Char('v')));
        assert_eq!(app.mode, Mode::Plot);
        app.handle_key(key(KeyCode::Char('h')));
        assert_eq!(app.plot, PlotKind::Histogram);
        app.handle_key(key(KeyCode::Char('h')));
        assert_eq!(app.plot, PlotKind::Heatmap);
        app.handle_key(key(KeyCode::Char('h')));
        assert_eq!(app.plot, PlotKind::Bars);
        app.handle_key(key(KeyCode::F(1)));
        assert!(app.show_help);
        app.handle_key(key(KeyCode::Esc));
        assert!(!app.show_help && !app.should_quit);
        app.handle_key(key(KeyCode::Char('q')));
        assert!(app.should_quit);
    }

    #[test]
    fn test_sort_keeps_selection() {
        let mut app = App::new(data(), String::from("report.json"), Colormap::Viridis);
        app.page_rows = 10;
        app.handle_key(key(KeyCode::Down));
        assert_eq!(app.selected(), Some(1));
        app.handle_key(key(KeyCode::Char('s')));
        assert_eq!(app.tables[0].order, vec![1, 0, 2]);
        assert_eq!(app.selected(), Some(1));
        app.handle_key(key(KeyCode::Char('s')));
        assert_eq!(app.tables[0].order, vec![2, 0, 1]);
        app.handle_key(key(KeyCode::Char('s')));
        assert_eq!(app.tables[0].order, vec![0, 1, 2]);
        // Missing values sort last.
        app.handle_key(key(KeyCode::Right));
        assert_eq!(app.column().unwrap().label, "GC %");
        app.handle_key(key(KeyCode::Char('s')));
        assert_eq!(app.tables[0].order, vec![2, 0, 1]);
        app.handle_key(key(KeyCode::Char('s')));
        assert_eq!(app.tables[0].order, vec![0, 2, 1]);
    }

    #[test]
    fn test_matrix_cursor() {
        let mut app = App::new(data(), String::from("report.json"), Colormap::Viridis);
        app.handle_key(key(KeyCode::Char('3')));
        assert_eq!(app.open_section(), Section::Matrix(0));
        assert!(app.table().is_none());
        app.page_rows = 2;
        app.page_columns = 2;
        app.handle_key(key(KeyCode::Right));
        app.handle_key(key(KeyCode::Right));
        app.handle_key(key(KeyCode::Down));
        assert_eq!((app.matrices[0].row, app.matrices[0].column), (1, 2));
        assert_eq!(app.matrices[0].column_offset, 1);
        app.handle_key(key(KeyCode::End));
        assert_eq!(app.matrices[0].row, 2);
        assert_eq!(app.matrices[0].row_offset, 1);
        app.handle_key(key(KeyCode::Char('v')));
        assert_eq!(app.mode, Mode::Table);
        assert!(app.message.is_some());
        // Sorting does not apply to a matrix.
        app.handle_key(key(KeyCode::Char('s')));
        assert_eq!(app.tables[0].sort, None);
    }

    #[test]
    fn test_scroll() {
        let mut app = App::new(data(), String::from("report.json"), Colormap::Viridis);
        app.page_rows = 2;
        app.handle_key(key(KeyCode::End));
        assert_eq!(app.tables[0].row, 2);
        assert_eq!(app.tables[0].offset, 1);
        app.handle_key(key(KeyCode::Home));
        assert_eq!(app.tables[0].offset, 0);
    }
}
