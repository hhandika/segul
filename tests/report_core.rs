//! Integration tests for unified `report.json`, checksums, verification,
//! replay, workflow lineage, backward compatibility, and methods export.

use std::fs;
use std::path::PathBuf;

use tempfile::TempDir;

use segul_core::report::checksum::{ChecksumMode, sha256_file};
use segul_core::report::paths::{
    common_base_dir, relative_posix, resolve_portable_path, to_posix_string,
};
use segul_core::report::spec::AlignConcatParams;
use segul_core::report::verify::{FileVerifyStatus, VerifyScope, verify_report};
use segul_core::report::{FrontendInfo, ReportContext, TaskKind, write_task_report};

#[test]
fn test_portable_posix_paths() {
    let dir = TempDir::new().unwrap();
    let a = dir.path().join("sub").join("a.fas");
    let b = dir.path().join("sub").join("nested").join("b.fas");
    let inputs = vec![a.clone(), b.clone()];
    let base = common_base_dir(&inputs);
    assert_eq!(base, dir.path().join("sub"));
    assert_eq!(relative_posix(&a, &base), "a.fas");
    assert_eq!(relative_posix(&b, &base), "nested/b.fas");
    assert!(!to_posix_string(&b).contains('\\'));
    assert_eq!(resolve_portable_path(&base, "nested/b.fas"), b);
}

#[test]
fn test_sha256_and_synthetic_verify_statuses() {
    let tmp = TempDir::new().unwrap();
    let in_dir = tmp.path().join("inputs");
    let out_dir = tmp.path().join("outputs");
    fs::create_dir_all(&in_dir).unwrap();
    fs::create_dir_all(&out_dir).unwrap();

    let in_file = in_dir.join("locus1.fas");
    let out_ok = out_dir.join("summary.csv");
    let out_missing = out_dir.join("removed.csv");
    fs::write(&in_file, b">t1\nACTG\n").unwrap();
    fs::write(&out_ok, b"col\n1\n").unwrap();
    fs::write(&out_missing, b"gone").unwrap();

    let params = AlignConcatParams {
        input_format: "fasta".to_string(),
        datatype: "dna".to_string(),
        output_format: "nexus".to_string(),
        partition_format: "nexus".to_string(),
        sort: false,
        prefix: Some("run".to_string()),
    };

    // 1. ChecksumMode::All verifies Ok, then Mismatch, then Missing.
    let ctx_all = ReportContext::new(FrontendInfo::cli("0.25.0"), ChecksumMode::All);
    let (rpath, report) = write_task_report::<_, serde_json::Value>(
        &ctx_all,
        TaskKind::AlignConcat,
        std::slice::from_ref(&in_file),
        Vec::new(),
        &params,
        &out_dir,
        vec![("table", out_ok.clone()), ("extra", out_missing.clone())],
        None,
        None,
    )
    .unwrap();
    assert_eq!(report.schema_version, 2);
    assert!(report.inputs.primary[0].sha256.is_some());
    assert!(
        report
            .outputs
            .files
            .iter()
            .all(|f| !f.path.ends_with("report.json"))
    );

    let vr_ok = verify_report(&rpath, None, None, VerifyScope::All).unwrap();
    assert!(vr_ok.passed);
    assert_eq!(vr_ok.mismatched + vr_ok.missing, 0);

    // Tamper output & remove missing file.
    fs::write(&out_ok, b"col\n9\n").unwrap();
    fs::remove_file(&out_missing).unwrap();
    let vr_fail = verify_report(&rpath, None, None, VerifyScope::Outputs).unwrap();
    assert!(!vr_fail.passed);
    assert_eq!(vr_fail.mismatched, 1);
    assert_eq!(vr_fail.missing, 1);
    let out_missing_res = vr_fail
        .outputs
        .iter()
        .find(|f| f.path == "removed.csv")
        .unwrap();
    let out_ok_res = vr_fail
        .outputs
        .iter()
        .find(|f| f.path == "summary.csv")
        .unwrap();
    assert_eq!(out_missing_res.status, FileVerifyStatus::Missing);
    assert!(matches!(
        out_ok_res.status,
        FileVerifyStatus::Mismatch { .. }
    ));

    // Scoped input verification should still pass since input wasn't touched.
    let vr_in = verify_report(&rpath, None, None, VerifyScope::Inputs).unwrap();
    assert!(vr_in.passed);
    assert_eq!(
        sha256_file(&in_file).unwrap(),
        report.inputs.primary[0].sha256.as_deref().unwrap()
    );

    // 2. ChecksumMode::None produces Unchecked when size matches, SizeMismatch when size changes.
    let out_unchecked = out_dir.join("plain.txt");
    fs::write(&out_unchecked, b"abcd").unwrap();
    let ctx_none = ReportContext::new(FrontendInfo::python("0.1.0"), ChecksumMode::None);
    let (rpath_none, _) = write_task_report::<_, serde_json::Value>(
        &ctx_none,
        TaskKind::AlignConcat,
        &[in_file],
        Vec::new(),
        &params,
        &out_dir,
        vec![("text", out_unchecked.clone())],
        None,
        Some("unchecked"),
    )
    .unwrap();
    let vr_unchecked = verify_report(&rpath_none, None, None, VerifyScope::All).unwrap();
    assert!(vr_unchecked.passed);
    assert_eq!(vr_unchecked.unchecked, 2);
    assert_eq!(vr_unchecked.outputs[0].status, FileVerifyStatus::Unchecked);

    fs::write(&out_unchecked, b"abcdefgh").unwrap();
    let vr_size = verify_report(&rpath_none, None, None, VerifyScope::Outputs).unwrap();
    assert!(!vr_size.passed);
    assert!(matches!(
        vr_size.outputs[0].status,
        FileVerifyStatus::SizeMismatch {
            expected: 4,
            actual: 8
        }
    ));
}

#[test]
fn test_align_summary_v2_and_v1_compat_and_replay() {
    use segul_core::core::align::summarize::AlignmentSummary;
    use segul_core::helper::finder::SeqFileFinder;
    use segul_core::helper::types::{DataType, InputFmt};
    use segul_core::report::RunReport;
    use segul_core::report::methods::{MethodsFormat, render_methods};
    use segul_core::report::replay::{ReplayOptions, replay_report};
    use segul_core::report::stats::ReportStatistics;
    use segul_core::stats::align_report::AlignReport;
    use segul_core::stats::distance_report::DistanceReport;

    let tmp = TempDir::new().unwrap();
    let out1 = tmp.path().join("run1");
    let out2 = tmp.path().join("run2");
    let files =
        SeqFileFinder::new(PathBuf::from("tests/files/gappy").as_path()).find(&InputFmt::Nexus);
    let ctx = ReportContext::new(FrontendInfo::cli("0.25.0"), ChecksumMode::All)
        .with_input_dir(Some(PathBuf::from("tests/files/gappy")));
    let mut summary =
        AlignmentSummary::new(&InputFmt::Nexus, &out1, 5, &DataType::Dna).with_report_context(ctx);
    let result = summary.summarize_all(&files, None);

    let rpath = result.report;
    assert!(rpath.ends_with("report.json"));
    let v2 = RunReport::from_path(&rpath).unwrap();
    assert_eq!(v2.schema_version, 2);
    assert_eq!(v2.kind, TaskKind::AlignSummary);
    assert_eq!(v2.run.checksum_mode, ChecksumMode::All);
    assert!(v2.inputs.primary.iter().all(|f| f.sha256.is_some()));

    // AlignReport::from_json reads v2 envelope seamlessly.
    let raw = fs::read_to_string(&rpath).unwrap();
    let parsed_v2 = AlignReport::from_json(&raw).unwrap();
    assert_eq!(parsed_v2.loci.len(), v2.inputs.primary.len());

    // Legacy v1 flat JSON (without envelope fields).
    let legacy = serde_json::to_string(&parsed_v2).unwrap();
    let parsed_v1 = AlignReport::from_json(&legacy).unwrap();
    assert_eq!(parsed_v1.loci.len(), parsed_v2.loci.len());

    // Verify checksums pass.
    let verify = verify_report(&rpath, None, None, VerifyScope::All).unwrap();
    assert!(verify.passed);

    // Replay into out2 and verify identical statistics.
    let replayed = replay_report(
        &rpath,
        &ReplayOptions {
            output_dir: Some(out2.clone()),
            input_base: Some(PathBuf::from("tests/files/gappy")),
            frontend: FrontendInfo::gui("0.6.0"),
            checksum_mode: Some(ChecksumMode::All),
            verify_inputs: true,
        },
    )
    .unwrap();
    let ReportStatistics::Align(rep2_stats) = replayed.typed_statistics().unwrap().unwrap() else {
        panic!("Expected Align statistics");
    };
    assert_eq!(rep2_stats.summary, parsed_v2.summary);
    assert_eq!(rep2_stats.loci.len(), parsed_v2.loci.len());
    assert_eq!(v2.parameters, replayed.parameters);

    // Legacy DistanceReport compat check using real legacy fixture.
    let dist_legacy_raw = fs::read_to_string("tests/files/plot/distance-report.json").unwrap();
    let dist_legacy = DistanceReport::from_json(&dist_legacy_raw).unwrap();
    assert!(!dist_legacy.model.is_empty());

    // Methods rendering in both Markdown and PlainText.
    let md = render_methods(std::slice::from_ref(&replayed), MethodsFormat::Markdown);
    let txt = render_methods(&[replayed], MethodsFormat::PlainText);
    assert!(md.contains("## Methods"));
    assert!(txt.contains("Methods"));
}

#[test]
fn test_lineage_workflow_export_and_replay_concat() {
    use segul_core::core::align::concat::AlignmentConcatenation;
    use segul_core::core::align::filter::{AlignmentFiltering, FilteringParameters};
    use segul_core::helper::finder::SeqFileFinder;
    use segul_core::helper::types::{DataType, InputFmt, OutputFmt, PartitionFmt};
    use segul_core::report::methods::{MethodsFormat, render_methods};
    use segul_core::report::replay::{ReplayOptions, replay_workflow};
    use segul_core::report::spec::{AlignFilterCriterion, AlignFilterParams, TaskParameters};
    use segul_core::report::workflow::WorkflowRecipe;

    let tmp = TempDir::new().unwrap();
    let step1_dir = tmp.path().join("step1_filter");
    let step2_dir = tmp.path().join("step2_concat");
    let wf_out_dir = tmp.path().join("workflow_replay");

    let input_dir = PathBuf::from("tests/files/gappy");
    let files = SeqFileFinder::new(&input_dir).find(&InputFmt::Nexus);
    let ctx = ReportContext::new(FrontendInfo::cli("0.25.0"), ChecksumMode::Outputs);

    // Step 1: filter
    let param = FilteringParameters::MinTax(2);
    let step1 =
        AlignmentFiltering::new(&files, &InputFmt::Nexus, &DataType::Dna, &step1_dir, &param)
            .with_report_context(ctx.clone())
            .filter_with_report(
                AlignFilterParams {
                    input_format: "nexus".to_string(),
                    datatype: "dna".to_string(),
                    criterion: AlignFilterCriterion::MinTaxCount { min_taxa: 2 },
                    concat: false,
                    output_format: None,
                    partition_format: None,
                    prefix: None,
                },
                &[],
            )
            .unwrap();
    assert!(!step1.outputs.files.is_empty());

    // Step 2: concat from step1 outputs -> auto-discovers parent report.
    let mut filtered_files = SeqFileFinder::new(&step1_dir).find(&InputFmt::Nexus);
    let prefix_path = step2_dir.join("concat");
    let mut concat = AlignmentConcatenation::new(
        &InputFmt::Nexus,
        &step2_dir,
        &OutputFmt::Nexus,
        &PartitionFmt::Nexus,
        &prefix_path,
    )
    .with_report_context(ctx);
    let step2 = concat
        .concat_with_report(&mut filtered_files, &DataType::Dna, Some("concat"))
        .unwrap();
    assert_eq!(step2.run.parent_reports.len(), 1);
    assert_eq!(step2.run.parent_reports[0].kind, TaskKind::AlignFilter);

    // Export workflow from step2 report chain and replay it.
    let step2_report_path = step2_dir.join("concat_report.json");
    let recipe =
        WorkflowRecipe::from_report_chain(&step2_report_path, Some("filter-and-concat")).unwrap();
    assert_eq!(recipe.steps.len(), 2);

    let wf_path = tmp.path().join("segul-workflow.json");
    recipe.write(&wf_path).unwrap();

    let replayed_reports = replay_workflow(
        &wf_path,
        &ReplayOptions {
            output_dir: Some(wf_out_dir),
            input_base: Some(input_dir),
            frontend: FrontendInfo::python("0.1.0"),
            checksum_mode: Some(ChecksumMode::All),
            verify_inputs: false,
        },
    )
    .unwrap();
    assert_eq!(replayed_reports.len(), 2);
    assert_eq!(replayed_reports[1].kind, TaskKind::AlignConcat);
    assert!(matches!(
        replayed_reports[1].typed_parameters().unwrap(),
        TaskParameters::AlignConcat(_)
    ));

    let methods = render_methods(&replayed_reports, MethodsFormat::Markdown);
    assert!(methods.contains("1. **"));
    assert!(methods.contains("2. **"));
}

#[test]
fn test_contig_report_compat_and_replay() {
    use segul_core::core::contig::summarize::ContigSummaryHandler;
    use segul_core::helper::finder::ContigFileFinder;
    use segul_core::helper::types::ContigFmt;
    use segul_core::report::replay::{ReplayOptions, replay_report};
    use segul_core::report::stats::ReportStatistics;
    use segul_core::stats::contig_report::ContigReport;

    let tmp = TempDir::new().unwrap();
    let out1 = tmp.path().join("c1");
    let out2 = tmp.path().join("c2");
    let files = ContigFileFinder::new(PathBuf::from("tests/files/contigs").as_path())
        .find(&ContigFmt::Fasta);
    let ctx = ReportContext::new(FrontendInfo::cli("0.25.0"), ChecksumMode::Outputs);
    let written = ContigSummaryHandler::new(&files, &ContigFmt::Fasta, &out1, None)
        .with_report_context(ctx)
        .summarize();
    let v2_json = fs::read_to_string(&written.report).unwrap();
    let c1 = ContigReport::from_json(&v2_json).unwrap();
    let c1_legacy = ContigReport::from_json(&serde_json::to_string(&c1).unwrap()).unwrap();
    assert_eq!(c1.summary, c1_legacy.summary);

    let replayed = replay_report(
        &written.report,
        &ReplayOptions {
            output_dir: Some(out2),
            input_base: None,
            frontend: FrontendInfo::cli("0.25.0"),
            checksum_mode: Some(ChecksumMode::Outputs),
            verify_inputs: false,
        },
    )
    .unwrap();
    let ReportStatistics::Contig(c2) = replayed.typed_statistics().unwrap().unwrap() else {
        panic!("Expected Contig statistics");
    };
    assert_eq!(c1.summary, c2.summary);
}

#[test]
fn test_read_report_compat_and_replay() {
    use segul_core::core::read::summarize::GenomicReadSummary;
    use segul_core::helper::types::{SeqReadFmt, SummaryMode};
    use segul_core::report::replay::{ReplayOptions, replay_report};
    use segul_core::report::stats::ReportStatistics;
    use segul_core::stats::fastq::report::ReadReport;

    let tmp = TempDir::new().unwrap();
    let out1 = tmp.path().join("r1");
    let out2 = tmp.path().join("r2");
    let mut files = vec![
        PathBuf::from("tests/files/raw/read_1.fastq"),
        PathBuf::from("tests/files/raw/read_2.fastq"),
    ];
    let ctx = ReportContext::new(FrontendInfo::cli("0.25.0"), ChecksumMode::Outputs);
    let written = GenomicReadSummary::new(
        &mut files,
        &SeqReadFmt::Auto,
        &SummaryMode::Complete,
        &out1,
        None,
    )
    .with_report_context(ctx)
    .try_summarize()
    .unwrap();
    let report_path = written.report.unwrap();
    let v2_json = fs::read_to_string(&report_path).unwrap();
    let r1 = ReadReport::from_json(&v2_json).unwrap();
    let r1_legacy = ReadReport::from_json(&serde_json::to_string(&r1).unwrap()).unwrap();
    assert_eq!(r1.samples.len(), r1_legacy.samples.len());

    let replayed = replay_report(
        &report_path,
        &ReplayOptions {
            output_dir: Some(out2),
            input_base: None,
            frontend: FrontendInfo::python("0.1.0"),
            checksum_mode: Some(ChecksumMode::All),
            verify_inputs: false,
        },
    )
    .unwrap();
    let ReportStatistics::Read(r2) = replayed.typed_statistics().unwrap().unwrap() else {
        panic!("Expected Read statistics");
    };
    assert_eq!(r1.samples.len(), r2.samples.len());
}

#[test]
fn test_maf_summary_replay() {
    use segul_core::core::maf::summary::MafSummary;
    use segul_core::helper::finder::MafFileFinder;
    use segul_core::report::replay::{ReplayOptions, replay_report};

    let tmp = TempDir::new().unwrap();
    let out1 = tmp.path().join("m1");
    let out2 = tmp.path().join("m2");
    let files = MafFileFinder::new(PathBuf::from("tests/files/maf").as_path()).find();
    let ctx = ReportContext::new(FrontendInfo::cli("0.25.0"), ChecksumMode::Outputs);
    let written = MafSummary::new(&files, &out1, None)
        .with_report_context(ctx)
        .try_summarize_with_report()
        .unwrap();
    assert_eq!(written.run_report.kind, TaskKind::MafSummary);

    let replayed = replay_report(
        &written.report,
        &ReplayOptions {
            output_dir: Some(out2),
            input_base: None,
            frontend: FrontendInfo::gui("0.6.0"),
            checksum_mode: Some(ChecksumMode::Outputs),
            verify_inputs: false,
        },
    )
    .unwrap();
    assert_eq!(replayed.statistics, written.run_report.statistics);
}

#[test]
fn test_sequence_id_replay() {
    use segul_core::core::sequence::id::SequenceID;
    use segul_core::helper::finder::SeqFileFinder;
    use segul_core::helper::types::{DataType, InputFmt};
    use segul_core::report::replay::{ReplayOptions, replay_report};

    let tmp = TempDir::new().unwrap();
    let out1 = tmp.path().join("s1");
    let out2 = tmp.path().join("s2");
    let files = SeqFileFinder::new(PathBuf::from("tests/files/alignments").as_path())
        .find(&InputFmt::Nexus);
    let ctx = ReportContext::new(FrontendInfo::cli("0.25.0"), ChecksumMode::All);
    let report1 = SequenceID::new(&files, &InputFmt::Nexus, &DataType::Dna, &out1, None)
        .with_report_context(ctx)
        .run_with_report(false)
        .unwrap();

    let replayed = replay_report(
        &out1.join("report.json"),
        &ReplayOptions {
            output_dir: Some(out2),
            input_base: None,
            frontend: FrontendInfo::cli("0.25.0"),
            checksum_mode: Some(ChecksumMode::All),
            verify_inputs: true,
        },
    )
    .unwrap();
    assert_eq!(report1.outputs.files.len(), replayed.outputs.files.len());
}
