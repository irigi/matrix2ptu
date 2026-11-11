"""Build the portable Windows zip of matrix2ptu (no installer needed).

    python tools/build_windows_zip.py [--out-dir dist] [--version 0.1.0]

The zip contains the official embeddable Python from python.org, matrix2ptu and
its dependencies (converter + validation readers) as Windows wheels, the
example configs and two launchers. The user unzips it anywhere and
double-clicks matrix2ptu_terminal.bat; nothing is installed on the computer.

Runs on any OS: the Windows wheels are fetched with pip's --platform option.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
import urllib.request
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PY_VERSION = "3.12.10"  # last 3.12 release with Windows binaries; must match the cp312 wheels
PY_URL = f"https://www.python.org/ftp/python/{PY_VERSION}/python-{PY_VERSION}-embed-amd64.zip"
PY_SHA256 = "4acbed6dd1c744b0376e3b1cf57ce906f9dc9e95e68824584c8099a63025a3c3"

# Runtime packages in the zip, pinned to the versions in requirements.txt.
# matplotlib (figures in tools/) and pytest are left out.
PACKAGES = ["numpy", "PyYAML", "ptufile", "tifffile", "tttrkit", "xarray", "scipy"]

LAUNCHERS = os.path.join(ROOT, "packaging", "windows")


def pinned_requirements():
    pins = {}
    with open(os.path.join(ROOT, "requirements.txt"), encoding="utf-8") as f:
        for line in f:
            m = re.match(r"^([A-Za-z0-9_.-]+)==(\S+)", line.strip())
            if m:
                pins[m.group(1).lower()] = m.group(2)
    missing = [p for p in PACKAGES if p.lower() not in pins]
    if missing:
        sys.exit(f"not pinned in requirements.txt: {missing}")
    return [f"{p}=={pins[p.lower()]}" for p in PACKAGES]


def project_version():
    with open(os.path.join(ROOT, "pyproject.toml"), encoding="utf-8") as f:
        return re.search(r'^version\s*=\s*"([^"]+)"', f.read(), re.M).group(1)


def fetch_python(cache_dir):
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, os.path.basename(PY_URL))
    if not os.path.exists(path):
        print(f"downloading {PY_URL}")
        urllib.request.urlretrieve(PY_URL, path)
    with open(path, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    if digest != PY_SHA256:
        os.remove(path)
        sys.exit(f"checksum mismatch for {path}: {digest}")
    return path


def pip(*args):
    cmd = [sys.executable, "-m", "pip", "--disable-pip-version-check", *args]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def build(out_dir, version):
    name = f"matrix2ptu-{version}-windows-x64"
    stage = os.path.join(out_dir, name)
    if os.path.exists(stage):
        shutil.rmtree(stage)
    pydir = os.path.join(stage, "python")
    site = os.path.join(pydir, "Lib", "site-packages")
    os.makedirs(site)

    # 1. embeddable Python; enable site-packages via its ._pth file
    with zipfile.ZipFile(fetch_python(os.path.join(out_dir, "cache"))) as z:
        z.extractall(pydir)
    pth = os.path.join(pydir, "python312._pth")
    with open(pth, "w", encoding="ascii", newline="\r\n") as f:
        f.write("python312.zip\n.\nLib\\site-packages\nimport site\n")

    # 2. dependencies as Windows wheels, then matrix2ptu itself (pure Python)
    pip("install", "--target", site, "--platform", "win_amd64", "--python-version", "3.12",
        "--implementation", "cp", "--only-binary=:all:", "--no-compile", *pinned_requirements())
    shutil.rmtree(os.path.join(site, "bin"), ignore_errors=True)  # entry scripts, unused
    had_build = os.path.exists(os.path.join(ROOT, "build"))
    pip("install", "--target", site, "--no-deps", "--no-compile", ROOT)
    shutil.rmtree(os.path.join(site, "bin"), ignore_errors=True)
    if not had_build:
        shutil.rmtree(os.path.join(ROOT, "build"), ignore_errors=True)

    # 3. launchers, example configs, README
    for fn in os.listdir(LAUNCHERS):
        shutil.copy2(os.path.join(LAUNCHERS, fn), stage)
    shutil.copytree(os.path.join(ROOT, "configs"), os.path.join(stage, "configs"))
    shutil.copy2(os.path.join(ROOT, "README.md"), stage)

    # 4. zip, with the top folder inside so unzipping gives one folder
    zpath = os.path.join(out_dir, name + ".zip")
    if os.path.exists(zpath):
        os.remove(zpath)
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for dirpath, _, files in os.walk(stage):
            for fn in sorted(files):
                full = os.path.join(dirpath, fn)
                z.write(full, os.path.relpath(full, out_dir))
    print(f"wrote {zpath} ({os.path.getsize(zpath) / 1e6:.0f} MB)")
    return zpath


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "dist"))
    ap.add_argument("--version", default=None, help="default: version in pyproject.toml")
    a = ap.parse_args()
    build(os.path.abspath(a.out_dir), a.version or project_version())


if __name__ == "__main__":
    main()
