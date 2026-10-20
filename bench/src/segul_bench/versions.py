"""Build the current SEGUL and download release binaries."""

from __future__ import annotations

import platform
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from .config import CACHE_DIR, REPO_DIR, workspace_version

RELEASE_URL = "https://github.com/hhandika/segul/releases/download/{tag}/{asset}"
CURRENT = "current"


@dataclass
class Version:
    id: str  # a release tag such as "v0.23.2", or "current"
    binary: Path
    source: str  # "release" or "build"
    number: tuple[int, ...]
    git: str | None = None

    def to_json(self) -> dict:
        data = {"id": self.id, "source": self.source, "version": ".".join(map(str, self.number))}
        if self.git:
            data["git"] = self.git
        return data


def parse_version(text: str) -> tuple[int, ...]:
    return tuple(int(part) for part in text.lstrip("v").split("."))


def release_asset() -> str:
    """Asset name of the release archive for this machine, as uploaded by
    .github/workflows/release.yml."""
    system = platform.system()
    machine = platform.machine().lower()
    arch = {"x86_64": "x86_64", "amd64": "x86_64", "aarch64": "arm64", "arm64": "arm64"}.get(machine)
    os_name = {"Linux": "Linux", "Darwin": "macOS"}.get(system)
    if arch is None or os_name is None:
        raise SystemExit(f"Release binaries are not supported on {system} {machine}")
    return f"segul-{os_name}-{arch}.tar.gz"


def fetch_release(tag: str, refresh: bool = False) -> Version:
    """Download a release binary once and reuse it from the cache."""
    binary = CACHE_DIR / "bin" / tag / "segul"
    if refresh or not binary.exists():
        asset = release_asset()
        url = RELEASE_URL.format(tag=tag, asset=asset)
        print(f"Downloading {tag} ({asset})")
        binary.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as tmp:
            archive = Path(tmp) / asset
            try:
                urllib.request.urlretrieve(url, archive)
            except Exception as err:
                raise SystemExit(f"Failed downloading {url}: {err}") from err
            with tarfile.open(archive) as tar:
                member = tar.extractfile("segul")
                if member is None:
                    raise SystemExit(f"No segul binary in {url}")
                with open(binary, "wb") as out:
                    shutil.copyfileobj(member, out)
        binary.chmod(0o755)
    return Version(id=tag, binary=binary, source="release", number=parse_version(tag))


def build_current(skip_build: bool = False) -> Version:
    """Build the working tree and copy the binary,
    so a rebuild during the run cannot change it."""
    built = REPO_DIR / "target" / "release" / "segul"
    if not skip_build:
        print("Building current working tree (cargo build --release)")
        subprocess.run(["cargo", "build", "--release"], cwd=REPO_DIR, check=True)
    if not built.exists():
        raise SystemExit(f"{built} not found. Run without --skip-build.")
    binary = CACHE_DIR / "bin" / CURRENT / "segul"
    binary.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(built, binary)
    git = subprocess.run(
        ["git", "describe", "--tags", "--dirty", "--always"],
        cwd=REPO_DIR,
        capture_output=True,
        text=True,
    ).stdout.strip()
    return Version(
        id=CURRENT,
        binary=binary,
        source="build",
        number=parse_version(workspace_version()),
        git=git or None,
    )


def resolve(ids: list[str], skip_build: bool = False, refresh: bool = False) -> list[Version]:
    """Resolve version ids in the given order."""
    versions = []
    for version_id in ids:
        if version_id == CURRENT:
            versions.append(build_current(skip_build))
        else:
            versions.append(fetch_release(version_id, refresh))
    return versions
