"""Streaming FASTQ quality summaries used by the pipeline and UI."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import random
import shutil
import subprocess
from collections import Counter
from pathlib import Path
from typing import TextIO


try:
    import numpy as np
except ImportError:  # Standalone CLI installations can use the same lookup table.
    np = None

_ERROR_TABLE = tuple(10 ** (-max(0, code - 33) / 10) for code in range(256))
_VECTOR_ERRORS = np.array(_ERROR_TABLE) if np is not None else None


def _is_gzip(path: Path) -> bool:
    # Detect by magic bytes: temporary files (mktemp) have no .gz suffix.
    try:
        with path.open("rb") as handle:
            return handle.read(2) == b"\x1f\x8b"
    except OSError:
        return path.name.casefold().endswith(".gz")


def _open_fastq(path: Path) -> TextIO:
    if _is_gzip(path):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _n50(length_counts: Counter[int], total_bases: int) -> int:
    target = (total_bases + 1) // 2
    cumulative = 0
    for length, count in sorted(length_counts.items(), reverse=True):
        cumulative += length * count
        if cumulative >= target:
            return length
    return 0


def _read_quality(quality: bytes) -> float:
    if np is not None:
        codes = np.frombuffer(quality, dtype=np.uint8)
        score = -10 * math.log10(float(_VECTOR_ERRORS[codes].mean()))
        # Preserve the previous summation order at reported Q10/20/30 boundaries.
        if all(abs(score - cutoff) > 1e-8 for cutoff in (10, 20, 30)):
            return score
    return -10 * math.log10(sum(map(_ERROR_TABLE.__getitem__, quality)) / len(quality))


HISTOGRAM_BIN = 1000
HISTOGRAM_MAX = 100_000
_MODEL_KEYS = (b"basecall_model_version_id=", b"model_version_id=")


def basecall_model_from_header(header: bytes | str) -> str:
    """Return the Dorado/Guppy basecall model recorded in a FASTQ header, if any."""
    raw = header.encode() if isinstance(header, str) else header
    for token in raw.split():
        for key in _MODEL_KEYS:
            if token.startswith(key):
                return token[len(key):].decode("utf-8", errors="replace")
    return ""


def basecall_tier(model: str) -> str:
    """'sup', 'hac', 'fast' or '' from a model name such as ..._hac@v5.2.0."""
    lowered = model.lower()
    for tier in ("sup", "hac", "fast"):
        if f"_{tier}@" in lowered or lowered.endswith(f"_{tier}") or f"_{tier}_" in lowered:
            return tier
    return ""


class _Stats:
    """Accumulates read length / mean-Q metrics in bounded memory."""

    def __init__(self, reference_length: int = 0) -> None:
        self.reference_length = reference_length
        self.reads = self.total_bases = 0
        self.quality_sum = 0.0
        self.q10 = self.q20 = self.q30 = 0
        self.length_counts: Counter[int] = Counter()

    def add(self, length: int, quality: float) -> None:
        self.reads += 1
        self.total_bases += length
        self.quality_sum += quality
        self.length_counts[length] += 1
        self.q10 += quality >= 10
        self.q20 += quality >= 20
        self.q30 += quality >= 30

    def result(self) -> dict[str, object]:
        reads = self.reads
        payload: dict[str, object] = {
            "reads": reads,
            "total_bases": self.total_bases,
            "mean_length": self.total_bases / reads if reads else 0.0,
            "n50": _n50(self.length_counts, self.total_bases),
            "mean_read_quality": self.quality_sum / reads if reads else 0.0,
            "q10_read_fraction": self.q10 / reads if reads else 0.0,
            "q20_read_fraction": self.q20 / reads if reads else 0.0,
            "q30_read_fraction": self.q30 / reads if reads else 0.0,
        }
        bins: Counter[int] = Counter()
        for length, count in self.length_counts.items():
            bins[min(length, HISTOGRAM_MAX) // HISTOGRAM_BIN * HISTOGRAM_BIN] += count
        payload["length_histogram"] = {
            "bin_size": HISTOGRAM_BIN,
            "counts": {str(start): bins[start] for start in sorted(bins)},
        }
        ref = self.reference_length
        if ref > 0 and reads:
            low, high = 0.9 * ref, 1.1 * ref
            short = full = longer = 0
            for length, count in self.length_counts.items():
                if length < low:
                    short += count
                elif length <= high:
                    full += count
                else:
                    longer += count
            payload.update(
                {
                    "reference_length": ref,
                    "full_length_read_fraction": full / reads,
                    "short_read_fraction": short / reads,
                    "long_read_fraction": longer / reads,
                    "n50_to_reference": payload["n50"] / ref,
                }
            )
        return payload


def _records(handle):
    reads = 0
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline()
        separator = handle.readline()
        quality = handle.readline()
        seq = sequence.rstrip(b"\r\n")
        qual = quality.rstrip(b"\r\n")
        reads += 1
        if not seq or not separator or len(seq) != len(qual):
            raise ValueError(f"Malformed FASTQ record near read {reads}")
        yield header, sequence, separator, quality, seq, qual


def _open_binary(path: Path):
    opener = gzip.open if _is_gzip(path) else open
    return opener(path, "rb")


class _CompressedWriter:
    """Write gzip output through pigz when available (much faster than zlib)."""

    def __init__(self, path: Path, threads: int = 1) -> None:
        self._handle = path.open("wb")
        self._process = None
        pigz = shutil.which("pigz")
        if pigz:
            self._process = subprocess.Popen(
                [pigz, "-1", "-p", str(max(1, threads)), "-c"],
                stdin=subprocess.PIPE,
                stdout=self._handle,
            )
            self._stream = self._process.stdin
        else:
            self._stream = gzip.GzipFile(fileobj=self._handle, mode="wb", compresslevel=1)

    def writelines(self, parts) -> None:
        self._stream.writelines(parts)

    def close(self) -> None:
        self._stream.close()
        if self._process is not None and self._process.wait() != 0:
            self._handle.close()
            raise OSError("pigz failed while writing filtered FASTQ")
        self._handle.close()


def summarize_fastq(path: Path, reference_length: int = 0) -> dict[str, object]:
    """Return read-level length and mean-Q metrics without loading reads in memory."""
    path = path.expanduser().resolve()
    stats = _Stats(reference_length)
    with _open_binary(path) as handle:
        for _header, _s, _p, _q, seq, qual in _records(handle):
            stats.add(len(seq), _read_quality(qual))
    return stats.result()


def summarize_and_filter(
    merged: Path,
    filtered_out: Path,
    min_length: int,
    min_quality: float,
    reference_length: int = 0,
    threads: int = 1,
) -> tuple[dict[str, object], dict[str, object], str]:
    """One pass over the merged FASTQ: merged metrics, filtered metrics, filtered FASTQ.

    The read-quality definition (mean of per-base error probabilities, converted
    back to Q) is the same one NanoFilt/chopper use for ``-q``.
    """
    merged_stats = _Stats(reference_length)
    filtered_stats = _Stats(reference_length)
    model = ""
    writer = _CompressedWriter(filtered_out, threads)
    try:
        with _open_binary(merged.expanduser().resolve()) as handle:
            for header, sequence, separator, quality, seq, qual in _records(handle):
                if not model and merged_stats.reads < 50:
                    model = basecall_model_from_header(header)
                length = len(seq)
                read_quality = _read_quality(qual)
                merged_stats.add(length, read_quality)
                if length >= min_length and read_quality >= min_quality:
                    filtered_stats.add(length, read_quality)
                    writer.writelines((header, sequence, separator, quality))
    finally:
        writer.close()
    return merged_stats.result(), filtered_stats.result(), model


def first_basecall_model(path: Path, max_reads: int = 50) -> str:
    try:
        with _open_binary(path) as handle:
            for index, record in enumerate(_records(handle)):
                model = basecall_model_from_header(record[0])
                if model or index >= max_reads:
                    return model
    except (OSError, ValueError):
        return ""
    return ""


def diagnostics(merged: dict, filtered: dict | None, model: str) -> dict[str, object]:
    """Plain-language hints about basecalling and library prep."""
    effective = filtered or merged
    tier = basecall_tier(model)
    hints: list[str] = []
    if tier in {"hac", "fast"}:
        hints.append(
            f"Basecall 모델이 {tier.upper()}입니다 ({model}). Mean read Q는 주로 basecall 모델이 "
            "결정하므로, POD5를 Dorado SUP 모델로 다시 basecalling(re-basecalling)하면 보통 Q20 이상으로 올라갑니다."
        )
    ratio = effective.get("n50_to_reference")
    full = effective.get("full_length_read_fraction")
    if isinstance(ratio, float) and ratio < 0.8:
        hints.append(
            "Read N50이 reference 길이의 80% 미만입니다. Whole-plasmid 데이터는 대부분 전장 read여야 "
            "합니다: 제한효소 자리가 여러 개인지, DNA 절단·분해, 짧은 오염 DNA(gDNA/RNA)를 확인하세요. "
            "Rapid Barcoding(transposase)으로 원형 pDNA를 그대로 쓰면 효소 선형화가 필요 없습니다."
        )
    if isinstance(full, float) and full < 0.3:
        hints.append(f"전장 read(reference 길이의 90-110%)가 {full:.0%}뿐입니다.")
    longer = effective.get("long_read_fraction")
    if isinstance(longer, float) and longer > 0.2:
        hints.append(
            f"Read의 {longer:.0%}가 reference보다 깁니다: plasmid multimer, 더 큰 오염 plasmid, "
            "또는 잘못된 reference일 수 있습니다."
        )
    return {"basecall_model": model, "basecall_tier": tier, "hints": hints}


def load_fastq_qc(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def sample_for_mapping(source: Path, output: Path, metrics: dict, target_bases: int) -> dict:
    """Uniform selection of whole records, preserving input order and using bounded memory."""
    reads = int(metrics["reads"])
    bases = int(metrics["total_bases"])
    if not reads or not bases:
        raise ValueError("No usable FASTQ reads for screening")
    # Convert desired input bases to a read count using the actual mean read length.
    count = min(reads, max(1, (target_bases * reads + bases - 1) // bases))
    summary = {"available_reads": reads, "selected_reads": count, "selected_bases": bases,
               "target_bases": target_bases, "sampled": count < reads, "seed": 42}
    if count == reads:
        return summary
    rng = random.Random(42)
    needed = count
    selected_bases = 0
    with _open_fastq(source) as handle, gzip.open(output, "wt", compresslevel=1) as out:
        for remaining in range(reads, 0, -1):
            header, seq, plus, quality = (handle.readline() for _ in range(4))
            if not quality:
                raise ValueError("FASTQ changed since QC")
            if rng.randrange(remaining) < needed:
                out.writelines((header, seq, plus, quality))
                selected_bases += len(seq.rstrip("\r\n"))
                needed -= 1
                if not needed:
                    break
    summary["selected_bases"] = selected_bases
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Write compact FASTQ QC metrics as JSON.")
    parser.add_argument("--merged", required=True, type=Path)
    parser.add_argument("--filtered", type=Path, help="already-filtered FASTQ to summarise")
    parser.add_argument("--filter-out", type=Path, help="filter merged reads in the same pass")
    parser.add_argument("--min-length", type=int, default=0)
    parser.add_argument("--min-quality", type=float, default=0.0)
    parser.add_argument("--reference-length", type=int, default=0)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--target-bases", type=int, default=0)
    parser.add_argument("--sample-out", type=Path)
    args = parser.parse_args()

    if args.filter_out is not None:
        merged, filtered, model = summarize_and_filter(
            args.merged, args.filter_out, args.min_length, args.min_quality,
            args.reference_length, args.threads,
        )
        filtered_source: Path | None = args.filter_out
    else:
        merged = summarize_fastq(args.merged, args.reference_length)
        filtered = None
        filtered_source = None
        if args.filtered and args.filtered.is_file():
            filtered = summarize_fastq(args.filtered, args.reference_length)
            filtered_source = args.filtered
        model = first_basecall_model(args.merged)
    payload: dict[str, object] = {"merged": merged, "filtered": filtered}
    payload["diagnostics"] = diagnostics(merged, filtered, model)
    if args.target_bases > 0:
        if args.sample_out is None:
            parser.error("--target-bases requires --sample-out")
        source = filtered_source if filtered is not None else args.merged
        metrics = filtered or merged
        payload["screening"] = sample_for_mapping(source, args.sample_out, metrics, args.target_bases)
        print("Mapping input:", json.dumps(payload["screening"]))
    for hint in payload["diagnostics"]["hints"]:
        print("QC hint:", hint)
    args.output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
