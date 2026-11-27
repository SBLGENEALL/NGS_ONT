"""Dorado re-basecalling jobs (GPU) that feed the batch analysis.

A job runs ``dorado basecaller`` on a POD5 folder, then ``dorado demux`` to
split barcodes, and finally arranges the FASTQ files as
``<output>/fastq_pass/barcodeNN/*.fastq.gz`` (plus ``unclassified``) so the
folder can be selected directly as the FASTQ source in the Batch tab.

Like batch jobs, the script runs detached from Streamlit, writes its own
log and exit code, and can be re-attached or stopped from a later session.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import uuid
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

KITS = (
    "SQK-RBK114-96",
    "SQK-RBK114-24",
    "SQK-NBD114-96",
    "SQK-NBD114-24",
    "",  # no barcoding
)

PID_FILE = "basecall.pid"
EXIT_CODE_FILE = "basecall.exit_code"
STOPPED_FILE = "basecall.stopped"
LOG_FILE = "basecall.log"
SETTINGS_FILE = "basecall_job.json"
_MARKER = "ont-basecall"


class BasecallError(RuntimeError):
    pass


@dataclass(frozen=True)
class BasecallSettings:
    pod5_dir: str
    output_dir: str
    model_path: str
    kit_name: str = "SQK-RBK114-96"
    device: str = "cuda:all"
    recursive: bool = True
    min_qscore: int = 0  # 0 = Dorado default


@dataclass(frozen=True)
class BasecallJob:
    job_dir: Path
    settings: BasecallSettings


def dorado_executable() -> str | None:
    configured = os.environ.get("ONT_DORADO", "").strip()
    if configured and Path(configured).is_file():
        return configured
    return shutil.which("dorado")


def models_root() -> Path:
    return Path(os.environ.get("ONT_DORADO_MODELS", "/data/dorado_models")).expanduser()


def available_models(root: Path | None = None) -> list[Path]:
    """Simplex basecall model folders (dna_*) under the models root."""
    root = models_root() if root is None else root
    if not root.is_dir():
        return []
    models = [
        path for path in root.iterdir()
        if path.is_dir() and path.name.startswith("dna_") and "@" in path.name
        and "_mods_" not in path.name and not re.search(r"_(5mC|6mA|4mC)", path.name)
    ]
    # SUP first, newest version first.
    return sorted(models, key=lambda p: ("_sup@" not in p.name, p.name), reverse=False)


def gpu_summary() -> str:
    smi = shutil.which("nvidia-smi")
    if not smi:
        return ""
    try:
        output = subprocess.run(
            [smi, "--query-gpu=name,memory.total", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=10, check=False,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ""
    return "; ".join(line.strip() for line in output.splitlines() if line.strip())


def count_pod5(folder: Path, recursive: bool = True) -> int:
    pattern = "**/*.pod5" if recursive else "*.pod5"
    try:
        return sum(1 for _ in folder.glob(pattern))
    except OSError:
        return 0


def default_output_dir(pod5_dir: Path, model_path: str) -> Path:
    tier = "sup" if "_sup@" in model_path else ("hac" if "_hac@" in model_path else "dorado")
    stamp = datetime.now().strftime("%Y%m%d_%H%M")
    return pod5_dir.parent / f"basecall_{tier}_{stamp}"


def _script(settings: BasecallSettings, dorado: str, job_dir: Path) -> str:
    q = shlex.quote
    out = Path(settings.output_dir)
    bam = job_dir / "calls.bam"
    basecaller = [dorado, "basecaller", "--device", settings.device]
    if settings.recursive:
        basecaller.append("--recursive")
    if settings.kit_name:
        basecaller += ["--kit-name", settings.kit_name]
    if settings.min_qscore > 0:
        basecaller += ["--min-qscore", str(settings.min_qscore)]
    basecaller += [settings.model_path, settings.pod5_dir]
    demux_dir = job_dir / "demux"
    lines = [
        "set -euo pipefail",
        f'echo "[basecall] $(date "+%F %T") start: {settings.model_path}"',
        f"mkdir -p {q(str(out))}",
        f"{' '.join(q(part) for part in basecaller)} > {q(str(bam))}",
        f'echo "[basecall] $(date "+%F %T") basecalling finished"',
    ]
    if settings.kit_name:
        lines += [
            f"rm -rf {q(str(demux_dir))}",
            f"{q(dorado)} demux --no-classify --emit-fastq --output-dir {q(str(demux_dir))} {q(str(bam))}",
        ]
    else:
        lines += [
            f"mkdir -p {q(str(demux_dir))}",
            f"samtools fastq -T '*' {q(str(bam))} > {q(str(demux_dir / 'reads.fastq'))}",
        ]
    lines += [
        f'echo "[basecall] $(date "+%F %T") arranging FASTQ by barcode"',
        f"python {q(str(Path(__file__).resolve()))} arrange {q(str(demux_dir))} {q(str(out))}",
        f'echo "[basecall] $(date "+%F %T") done: {settings.output_dir}"',
    ]
    return "\n".join(lines) + "\n"


def prepare_basecall_job(run_root: Path, settings: BasecallSettings) -> BasecallJob:
    dorado = dorado_executable()
    if not dorado:
        raise BasecallError("dorado 실행 파일을 찾을 수 없습니다 (PATH 또는 ONT_DORADO).")
    pod5 = Path(settings.pod5_dir)
    if not pod5.is_dir() or count_pod5(pod5, settings.recursive) == 0:
        raise BasecallError("선택한 폴더에 POD5 파일이 없습니다.")
    if not Path(settings.model_path).is_dir():
        raise BasecallError(f"Dorado 모델 폴더가 없습니다: {settings.model_path}")
    if settings.kit_name and settings.kit_name not in KITS:
        raise BasecallError(f"지원하지 않는 kit: {settings.kit_name}")
    output = Path(settings.output_dir)
    if output.exists() and any(output.iterdir()):
        raise BasecallError(f"출력 폴더가 이미 있고 비어 있지 않습니다: {output}")
    try:
        output.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise BasecallError(f"출력 폴더를 만들 수 없습니다: {output} ({exc})") from exc
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    job_dir = run_root.resolve() / "basecall" / f"{stamp}_{uuid.uuid4().hex[:8]}"
    job_dir.mkdir(parents=True)
    (job_dir / SETTINGS_FILE).write_text(json.dumps(asdict(settings), ensure_ascii=False, indent=2), encoding="utf-8")
    (job_dir / "run_basecall.sh").write_text(_script(settings, dorado, job_dir), encoding="utf-8")
    return BasecallJob(job_dir, settings)


def load_basecall_job(job_dir: Path) -> BasecallJob:
    data = json.loads((job_dir / SETTINGS_FILE).read_text(encoding="utf-8"))
    return BasecallJob(job_dir, BasecallSettings(**data))


def start_basecall_job(job: BasecallJob) -> int:
    if basecall_job_state(job.job_dir) == "running":
        raise BasecallError("이 basecalling 작업은 이미 실행 중입니다.")
    wrapper = 'bash "$1"; rc=$?; printf "%s\\n" "$rc" > "$2"; exit "$rc"'
    with (job.job_dir / LOG_FILE).open("wb") as writer:
        process = subprocess.Popen(
            ["bash", "-c", wrapper, _MARKER, str(job.job_dir / "run_basecall.sh"), str(job.job_dir / EXIT_CODE_FILE)],
            cwd=job.job_dir, stdin=subprocess.DEVNULL, stdout=writer,
            stderr=subprocess.STDOUT, start_new_session=True,
        )
    (job.job_dir / PID_FILE).write_text(str(process.pid), encoding="utf-8")
    _CHILDREN[process.pid] = process
    return process.pid


_CHILDREN: dict[int, subprocess.Popen] = {}


def _alive(pid: int) -> bool:
    child = _CHILDREN.get(pid)
    if child is not None:
        return child.poll() is None
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    stat = Path(f"/proc/{pid}/stat")
    if stat.exists():
        try:
            if stat.read_text().rsplit(")", 1)[1].split()[0] in {"Z", "X"}:
                return False
            return _MARKER.encode() in Path(f"/proc/{pid}/cmdline").read_bytes()
        except (OSError, IndexError):
            return True
    return True


def basecall_job_state(job_dir: Path) -> str:
    exit_file = job_dir / EXIT_CODE_FILE
    if exit_file.is_file():
        try:
            return "completed" if int(exit_file.read_text().strip() or 1) == 0 else "failed"
        except ValueError:
            return "failed"
    pid_file = job_dir / PID_FILE
    if not pid_file.is_file():
        return "not_started"
    try:
        pid = int(pid_file.read_text().strip())
    except ValueError:
        return "interrupted"
    if _alive(pid):
        return "running"
    return "stopped" if (job_dir / STOPPED_FILE).is_file() else "interrupted"


def basecall_jobs(run_root: Path) -> list[BasecallJob]:
    root = run_root / "basecall"
    if not root.is_dir():
        return []
    dirs = sorted((p for p in root.iterdir() if (p / SETTINGS_FILE).is_file()),
                  key=lambda p: p.stat().st_mtime, reverse=True)
    return [load_basecall_job(path) for path in dirs]


def stop_basecall_job(job: BasecallJob) -> bool:
    if basecall_job_state(job.job_dir) != "running":
        return False
    pid = int((job.job_dir / PID_FILE).read_text().strip())
    (job.job_dir / STOPPED_FILE).write_text(datetime.now().isoformat(), encoding="utf-8")
    try:
        os.killpg(pid, signal.SIGTERM)
    except ProcessLookupError:
        return False
    return True


def basecall_log_tail(job: BasecallJob, lines: int = 15) -> list[str]:
    path = job.job_dir / LOG_FILE
    try:
        with path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            handle.seek(max(0, handle.tell() - 64 * 1024))
            text = handle.read().decode("utf-8", errors="replace")
    except OSError:
        return []
    # Dorado draws progress bars with carriage returns; keep the latest state.
    cleaned = [segment.split("\r")[-1] for segment in text.split("\n")]
    return [line for line in cleaned if line.strip()][-lines:]


_BARCODE = re.compile(r"(barcode\d+|unclassified)", re.IGNORECASE)


def arrange_fastq(demux_dir: Path, output_dir: Path) -> dict[str, int]:
    """Move Dorado demux FASTQ into output/fastq_pass/<barcode>/ and gzip them.

    Works with both the flat (``KIT_barcode01.fastq``) and the nested
    (``.../fastq_pass/barcode01/*.fastq``) layouts of different Dorado versions.
    """
    counts: dict[str, int] = {}
    files = sorted(
        path for path in demux_dir.rglob("*")
        if path.is_file() and re.search(r"\.(fastq|fq)(\.gz)?$", path.name, re.IGNORECASE)
    )
    pigz = shutil.which("pigz")
    for index, path in enumerate(files, 1):
        match = None
        for part in reversed(path.relative_to(demux_dir).parts):
            match = _BARCODE.search(part)
            if match:
                break
        group = match.group(1).lower() if match else "unclassified"
        target_dir = output_dir / "fastq_pass" / group
        target_dir.mkdir(parents=True, exist_ok=True)
        stem = re.sub(r"\.(fastq|fq)(\.gz)?$", "", path.name, flags=re.IGNORECASE)
        target = target_dir / f"{stem}_{index:04d}.fastq.gz"
        if path.name.lower().endswith(".gz"):
            shutil.move(str(path), target)
        else:
            with target.open("wb") as out:
                command = [pigz, "-c", "-1", str(path)] if pigz else ["gzip", "-c", "-1", str(path)]
                subprocess.run(command, stdout=out, check=True)
            path.unlink()
        counts[group] = counts.get(group, 0) + 1
    summary = output_dir / "basecall_summary.json"
    summary.write_text(json.dumps({"files_per_barcode": counts}, indent=2), encoding="utf-8")
    return counts


if __name__ == "__main__":
    import sys

    if len(sys.argv) == 4 and sys.argv[1] == "arrange":
        result = arrange_fastq(Path(sys.argv[2]), Path(sys.argv[3]))
        print("[basecall] barcodes:", ", ".join(f"{k}={v}" for k, v in sorted(result.items())) or "none")
    else:
        raise SystemExit("usage: basecall.py arrange DEMUX_DIR OUTPUT_DIR")
