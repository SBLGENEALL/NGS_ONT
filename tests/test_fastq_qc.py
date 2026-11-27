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


class FastqDiagnosticsTests(unittest.TestCase):
    def test_basecall_model_and_tier_from_header(self):
        header = b"@read1 runid=x ch=1 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.2.0\n"
        model = fastq_qc.basecall_model_from_header(header)
        self.assertEqual(model, "dna_r10.4.1_e8.2_400bps_hac@v5.2.0")
        self.assertEqual(fastq_qc.basecall_tier(model), "hac")
        self.assertEqual(fastq_qc.basecall_tier("dna_r10.4.1_e8.2_400bps_sup@v5.2.0"), "sup")
        self.assertEqual(fastq_qc.basecall_model_from_header(b"@plain\n"), "")

    def test_single_pass_filter_matches_two_pass_and_reports_full_length(self):
        with tempfile.TemporaryDirectory() as tmp:
            merged = Path(tmp) / "merged.fastq.gz"
            model = "basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.2.0"
            records = [
                f"@full1 {model}\n{'A' * 1000}\n+\n{'5' * 1000}\n",   # Q20, full length
                f"@full2 {model}\n{'C' * 1050}\n+\n{'+' * 1050}\n",   # Q10, full length
                f"@short {model}\n{'G' * 300}\n+\n{'I' * 300}\n",     # too short
                f"@lowq {model}\n{'T' * 1000}\n+\n{'%' * 1000}\n",    # Q4, fails -q 10
            ]
            merged.write_bytes(gzip.compress("".join(records).encode()))
            out = Path(tmp) / "filtered.fastq.gz"
            m, f, found = fastq_qc.summarize_and_filter(merged, out, 500, 10, reference_length=1000)
            self.assertEqual(found, "dna_r10.4.1_e8.2_400bps_hac@v5.2.0")
            self.assertEqual(m["reads"], 4)
            self.assertEqual(f["reads"], 2)
            self.assertEqual(summarize_fastq(out, 1000), f)
            self.assertEqual(summarize_fastq(merged, 1000), m)
            self.assertAlmostEqual(f["full_length_read_fraction"], 1.0)
            self.assertAlmostEqual(m["short_read_fraction"], 0.25)
            self.assertEqual(f["length_histogram"]["counts"], {"1000": 2})
            hints = fastq_qc.diagnostics(m, f, found)["hints"]
            self.assertTrue(any("SUP" in hint for hint in hints))

    def test_fragmented_library_hint(self):
        stats = {"n50": 4500, "n50_to_reference": 0.3, "full_length_read_fraction": 0.1}
        hints = fastq_qc.diagnostics(stats, None, "dna_r10.4.1_e8.2_400bps_sup@v5.2.0")["hints"]
        self.assertTrue(any("N50" in hint for hint in hints))
        self.assertFalse(any("re-basecalling" in hint for hint in hints))
