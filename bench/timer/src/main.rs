//! Run a command and print its wall time, CPU time, peak memory and exit
//! code as one JSON line on stdout.
//!
//! Peak memory (`ru_maxrss`) of a child includes the memory of the process
//! that spawned it, because Linux keeps the high-water mark across `exec`.
//! Spawning from this small process keeps that floor at about 2 MB instead
//! of the size of the Python benchmark runner.
//!
//! Usage: segul-bench-timer <command> [args...]
//! The command's stdin and stdout are discarded and its stderr is inherited.

use std::os::unix::process::ExitStatusExt;
use std::process::{Command, ExitCode, Stdio};
use std::time::Instant;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some((program, rest)) = args.split_first() else {
        eprintln!("Usage: segul-bench-timer <command> [args...]");
        return ExitCode::from(2);
    };

    let start = Instant::now();
    let status = Command::new(program)
        .args(rest)
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .status();
    let wall = start.elapsed().as_secs_f64();
    let status = match status {
        Ok(status) => status,
        Err(err) => {
            eprintln!("Failed to run {program}: {err}");
            return ExitCode::from(2);
        }
    };

    // The command is the only child, so the children's usage is its usage.
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    unsafe { libc::getrusage(libc::RUSAGE_CHILDREN, &mut usage) };
    let exit_code = status
        .code()
        .unwrap_or_else(|| -status.signal().unwrap_or(1));

    println!(
        "{{\"wall_s\": {wall:.6}, \"user_s\": {:.6}, \"sys_s\": {:.6}, \"max_rss_mb\": {:.3}, \"exit_code\": {exit_code}}}",
        seconds(usage.ru_utime),
        seconds(usage.ru_stime),
        max_rss_mb(usage.ru_maxrss),
    );
    ExitCode::SUCCESS
}

fn seconds(time: libc::timeval) -> f64 {
    time.tv_sec as f64 + time.tv_usec as f64 / 1e6
}

fn max_rss_mb(max_rss: libc::c_long) -> f64 {
    // macOS reports bytes, Linux reports kilobytes.
    if cfg!(target_os = "macos") {
        max_rss as f64 / (1024.0 * 1024.0)
    } else {
        max_rss as f64 / 1024.0
    }
}
