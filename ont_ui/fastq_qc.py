"""Streaming FASTQ quality summaries used by the pipeline and UI."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import random
from functools import lru_cache
from collections import Counter
from pathlib import Path
from typing import TextIO


def _open_fastq(path: Path) -> TextIO:
    if path.name.casefold().endswith(".gz"):
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


@lru_cache(maxsize=128)
def _error_probability(char: str) -> float:
    return 10 ** (-max(0, ord(char) - 33) / 10)


def summarize_fastq(path: Path) -> dict[str, int | float]:
    """Return read-level length and mean-Q metrics without loading reads in memory."""
    path = path.expanduser().resolve()
    reads = total_bases = 0
    read_quality_sum = 0.0
    q10 = q20 = q30 = 0
    length_counts: Counter[int] = Counter()

    with _open_fastq(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip("\r\n")
            separator = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            if not sequence or not separator or len(sequence) != len(quality):
                raise ValueError(f"Malformed FASTQ record in {path.name} near read {reads + 1}")

            mean_error_probability = sum(map(_error_probability, quality)) / len(quality)
            read_quality = -10 * math.log10(mean_error_probability)
            length = len(sequence)
            reads += 1
            total_bases += length
            read_quality_sum += read_quality
            length_counts[length] += 1
            q10 += read_quality >= 10
            q20 += read_quality >= 20
            q30 += read_quality >= 30

    return {
        "reads": reads,
        "total_bases": total_bases,
        "mean_length": total_bases / reads if reads else 0.0,
        "n50": _n50(length_counts, total_bases),
        "mean_read_quality": read_quality_sum / reads if reads else 0.0,
        "q10_read_fraction": q10 / reads if reads else 0.0,
        "q20_read_fraction": q20 / reads if reads else 0.0,
        "q30_read_fraction": q30 / reads if reads else 0.0,
    }


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
    parser.add_argument("--filtered", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--target-bases", type=int, default=0)
    parser.add_argument("--sample-out", type=Path)
    args = parser.parse_args()

    payload: dict[str, object] = {"merged": summarize_fastq(args.merged)}
    if args.filtered and args.filtered.is_file():
        payload["filtered"] = summarize_fastq(args.filtered)
    else:
        payload["filtered"] = None
    if args.target_bases > 0:
        if args.sample_out is None:
            parser.error("--target-bases requires --sample-out")
        source = args.filtered if payload["filtered"] is not None else args.merged
        metrics = payload["filtered"] or payload["merged"]
        payload["screening"] = sample_for_mapping(source, args.sample_out, metrics, args.target_bases)
        print("Mapping input:", json.dumps(payload["screening"]))
    args.output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
