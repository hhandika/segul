use std::ffi::OsStr;
use std::path::{Path, PathBuf};

use crate::handler::sequence::id::Id;

use super::args::SequenceIdArgs;
use super::{collect_paths, AlignSeqInput, AlignSeqPrint, InputCli, OutputCli};

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

    pub(in crate::cli) fn find(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let task_desc = "IDs finding";
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);

        AlignSeqPrint::new(
            &self.input_dir,
            &input_fmt,
            &datatype,
            task_desc,
            files.len(),
        )
        .print();

        let output = self.args.output.with_extension("txt");
        self.check_output_file_exist(&output, self.args.io.force);
        let id = Id::new(&output, &input_fmt, &datatype);
        if self.args.map {
            let map_fname = self.create_map_fname(&output);
            self.check_output_file_exist(&map_fname, self.args.io.force);
            id.map_id(&files, &map_fname);
        } else {
            id.generate_id(&files);
        }
    }

    fn create_map_fname(&self, output: &Path) -> PathBuf {
        let parent = output.parent().expect("Failed getting parent dir");
        let fstem = output
            .file_stem()
            .and_then(OsStr::to_str)
            .expect("Failed getting file stem for mapping IDs");
        parent.join(format!("{}_map", fstem)).with_extension("csv")
    }
}
