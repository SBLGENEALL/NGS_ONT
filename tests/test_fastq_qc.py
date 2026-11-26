import gzip
import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ont_ui import fastq_qc
from ont_ui.fastq_qc import load_fastq_qc, sample_for_mapping, summarize_fastq


class FastqQcTests(unittest.TestCase):
    def test_summarize_fastq_reports_read_level_qscore_fractions(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reads.fastq.gz"
            records = (
                "@q40\nAAAA\n+\nIIII\n"
                "@q20\nAAAAAA\n+\n555555\n"
                "@q10\nAA\n+\n++\n"
            )
            path.write_bytes(gzip.compress(records.encode("utf-8")))

            result = summarize_fastq(path)

            self.assertEqual(result["reads"], 3)
            self.assertEqual(result["total_bases"], 12)
            self.assertEqual(result["n50"], 6)
            self.assertAlmostEqual(result["mean_read_quality"], 70 / 3)
            self.assertAlmostEqual(result["q10_read_fraction"], 1.0)
            self.assertAlmostEqual(result["q20_read_fraction"], 2 / 3)
            self.assertAlmostEqual(result["q30_read_fraction"], 1 / 3)

    def test_lookup_quality_preserves_reported_cutoffs(self):
        for quality in (b"+" * 6000, b"5" * 6000, b"?" * 6000, b"!5I" * 2000):
            error_sum = sum(10 ** (-(q - 33) / 10) for q in quality)
            expected = -10 * math.log10(error_sum / len(quality))
            fast = fastq_qc._read_quality(quality)
            with patch.object(fastq_qc, "np", None):
                portable = fastq_qc._read_quality(quality)
            for actual in (fast, portable):
                self.assertAlmostEqual(actual, expected, places=9)
                for threshold in (10, 20, 30):
                    self.assertEqual(actual >= threshold, expected >= threshold)

    def test_screening_selects_reproducibly_and_skips_low_depth(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input.fastq"
            source.write_text("".join(f"@{i}\nACGT\n+\nIIII\n" for i in range(100)))
            metrics = summarize_fastq(source)
            output = Path(tmp) / "sample.fastq.gz"
            self.assertFalse(sample_for_mapping(source, output, metrics, 500)["sampled"])
            self.assertFalse(output.exists())
            summary = sample_for_mapping(source, output, metrics, 80)
            self.assertEqual(summary["selected_reads"], 20)
            self.assertEqual(summary["selected_bases"], 80)
            first = gzip.decompress(output.read_bytes())
            self.assertEqual(summarize_fastq(output)["reads"], 20)
            sample_for_mapping(source, output, metrics, 80)
            self.assertEqual(first, gzip.decompress(output.read_bytes()))

    def test_load_fastq_qc_returns_empty_dict_for_invalid_or_missing_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "qc.json"
            self.assertEqual(load_fastq_qc(path), {})
            path.write_text("not-json", encoding="utf-8")
            self.assertEqual(load_fastq_qc(path), {})
            path.write_text(json.dumps({"merged": {"reads": 2}}), encoding="utf-8")
            self.assertEqual(load_fastq_qc(path)["merged"]["reads"], 2)


if __name__ == "__main__":
    unittest.main()
