use std::ffi::OsStr;
use std::path::{Path, PathBuf};

use crate::handler::sequence::id::Id;
use crate::helper::files::create_output_fname_for_text;
use crate::helper::logger::AlignSeqLogger;

use super::args::SequenceIdArgs;
use super::{collect_paths, AlignSeqInput, InputCli, OutputCli};

impl InputCli for IdParser<'_> {}
impl OutputCli for IdParser<'_> {}
impl AlignSeqInput for IdParser<'_> {}

pub(in crate::cli) struct IdParser<'a> {
    args: &'a SequenceIdArgs,
    input_dir: Option<PathBuf>,
}

impl<'a> IdParser<'a> {
    pub(in crate::cli) fn new(args: &'a SequenceIdArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn extract(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        let log = AlignSeqLogger::new(
            self.input_dir.as_deref(),
            &input_fmt,
            &datatype,
            files.len(),
        );
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        let output = create_output_fname_for_text(&self.args.output, &self.args.prefix);
        let id = Id::new(&input_fmt, &datatype, &output);
        if self.args.map {
            let task = "Sequence ID Mapping";
            log.log(task);
            let map_fname = self.create_map_fname(&output);
            id.map_id(&files, &map_fname);
        } else {
            let task = "Sequence ID Generation";
            log.log(task);
            id.generate_id(&files);
        }
    }

    fn create_map_fname(&self, output: &Path) -> PathBuf {
        let parent = output.parent().expect("Failed getting parent dir");
        let file_stem = output
            .file_stem()
            .and_then(OsStr::to_str)
            .expect("Failed getting file stem for mapping IDs");
        parent
            .join(format!("{}_map", file_stem))
            .with_extension("csv")
    }
}
