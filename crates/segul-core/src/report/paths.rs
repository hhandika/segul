//! Portable POSIX path helpers for reproducible reports across platforms.

use std::path::{Component, Path, PathBuf};

/// Normalize a path to POSIX forward slashes using [`Component`]s.
pub fn to_posix_string(path: &Path) -> String {
    let mut parts: Vec<String> = Vec::new();
    let mut has_root = false;
    let mut prefix_str: Option<String> = None;

    for comp in path.components() {
        match comp {
            Component::Prefix(p) => {
                prefix_str = Some(p.as_os_str().to_string_lossy().into_owned());
            }
            Component::RootDir => {
                has_root = true;
            }
            Component::CurDir => {
                if parts.is_empty() && !has_root && prefix_str.is_none() {
                    parts.push(String::from("."));
                }
            }
            Component::ParentDir => parts.push(String::from("..")),
            Component::Normal(s) => parts.push(s.to_string_lossy().into_owned()),
        }
    }

    let mut out = String::new();
    if let Some(prefix) = prefix_str {
        out.push_str(&prefix);
    }
    if has_root {
        out.push('/');
    }
    out.push_str(&parts.join("/"));
    if out.is_empty() {
        String::from(".")
    } else {
        out
    }
}

/// Find the longest common parent directory across `paths`.
///
/// Returns `.` when `paths` is empty or shares no common directory prefix.
pub fn common_base_dir(paths: &[PathBuf]) -> PathBuf {
    let parents: Vec<Vec<Component<'_>>> = paths
        .iter()
        .map(|p| {
            let parent = p.parent().unwrap_or_else(|| Path::new(""));
            parent.components().collect::<Vec<_>>()
        })
        .collect();

    let Some(first) = parents.first() else {
        return PathBuf::from(".");
    };

    let mut common_len = first.len();
    for comps in &parents[1..] {
        let shared = first
            .iter()
            .zip(comps.iter())
            .take_while(|(a, b)| a == b)
            .count();
        common_len = common_len.min(shared);
    }

    if common_len == 0 {
        return PathBuf::from(".");
    }

    let mut base = PathBuf::new();
    for comp in &first[..common_len] {
        base.push(comp.as_os_str());
    }
    if base.as_os_str().is_empty() {
        PathBuf::from(".")
    } else {
        base
    }
}

/// Express `path` relative to `base` using POSIX forward slashes.
///
/// Falls back to normalized `path` when `path` and `base` have incompatible
/// prefixes (e.g., different drive letters on Windows) or mixed relative/absolute
/// roots.
pub fn relative_posix(path: &Path, base: &Path) -> String {
    let abs_path = absolutize_best_effort(path);
    let abs_base = absolutize_best_effort(base);

    let p_comps: Vec<Component<'_>> = clean_components(&abs_path);
    let b_comps: Vec<Component<'_>> = clean_components(&abs_base);

    // Different prefixes or one rooted / one unrooted cannot be made relative.
    let same_prefix = match (p_comps.first(), b_comps.first()) {
        (Some(Component::Prefix(a)), Some(Component::Prefix(b))) => a == b,
        (Some(Component::Prefix(_)), _) | (_, Some(Component::Prefix(_))) => false,
        _ => true,
    };
    if !same_prefix || abs_path.is_absolute() != abs_base.is_absolute() {
        return to_posix_string(path);
    }

    let common = p_comps
        .iter()
        .zip(b_comps.iter())
        .take_while(|(a, b)| a == b)
        .count();

    let mut rel_parts: Vec<String> = Vec::new();
    for comp in &b_comps[common..] {
        if matches!(comp, Component::Normal(_)) {
            rel_parts.push(String::from(".."));
        }
    }
    for comp in &p_comps[common..] {
        match comp {
            Component::Normal(s) => rel_parts.push(s.to_string_lossy().into_owned()),
            Component::ParentDir => rel_parts.push(String::from("..")),
            _ => {}
        }
    }

    if rel_parts.is_empty() {
        String::from(".")
    } else {
        rel_parts.join("/")
    }
}

/// Resolve a portable POSIX path string against `base`.
///
/// Relative segments (including `..`) are resolved onto `base`. Absolute POSIX
/// or Windows-prefixed strings are returned as-is.
pub fn resolve_portable_path(base: &Path, portable: &str) -> PathBuf {
    let raw = Path::new(portable);
    if raw.is_absolute() {
        return raw.components().collect();
    }
    let mut out: PathBuf = if base == Path::new(".") {
        PathBuf::new()
    } else {
        base.components().collect()
    };
    for part in portable.split('/') {
        match part {
            "" | "." => {}
            ".." => {
                if !out.pop() {
                    out.push("..");
                }
            }
            segment => out.push(segment),
        }
    }
    if out.as_os_str().is_empty() {
        PathBuf::from(".")
    } else {
        out
    }
}

/// Resolve a recorded `base_dir` (which may itself be relative) against a
/// fallback directory (typically the report's parent directory or CWD).
pub fn resolve_base_dir(recorded_base: &str, anchor: Option<&Path>) -> PathBuf {
    let base_path = Path::new(recorded_base);
    if base_path.is_absolute() {
        return base_path.components().collect();
    }
    match anchor {
        Some(dir) if !dir.as_os_str().is_empty() => resolve_portable_path(dir, recorded_base),
        _ => resolve_portable_path(Path::new("."), recorded_base),
    }
}

fn absolutize_best_effort(path: &Path) -> PathBuf {
    if path.is_absolute() {
        return normalize_lexical(path);
    }
    match std::env::current_dir() {
        Ok(cwd) => normalize_lexical(&cwd.join(path)),
        Err(_) => normalize_lexical(path),
    }
}

fn normalize_lexical(path: &Path) -> PathBuf {
    let mut out = PathBuf::new();
    for comp in path.components() {
        match comp {
            Component::CurDir => {}
            Component::ParentDir => {
                if !out.pop() {
                    out.push(comp.as_os_str());
                }
            }
            other => out.push(other.as_os_str()),
        }
    }
    if out.as_os_str().is_empty() {
        PathBuf::from(".")
    } else {
        out
    }
}

fn clean_components(path: &Path) -> Vec<Component<'_>> {
    path.components()
        .filter(|c| !matches!(c, Component::CurDir))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_common_base_and_relative_posix() {
        let a = Path::new("tests").join("files").join("a.nex");
        let b = Path::new("tests").join("files").join("sub").join("b.nex");
        let base = common_base_dir(&[a.clone(), b.clone()]);
        assert_eq!(to_posix_string(&base), "tests/files");
        assert_eq!(relative_posix(&a, &base), "a.nex");
        assert_eq!(relative_posix(&b, &base), "sub/b.nex");
        assert_eq!(resolve_portable_path(&base, "sub/b.nex"), b);
    }

    #[test]
    fn test_relative_parent_traversal() {
        let base = Path::new("data").join("alignments");
        let aux = Path::new("data").join("groups.csv");
        assert_eq!(relative_posix(&aux, &base), "../groups.csv");
        assert_eq!(resolve_portable_path(&base, "../groups.csv"), aux);
    }
}
