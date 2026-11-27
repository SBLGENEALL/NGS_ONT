import gzip
import os
import stat
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from ont_ui import basecall


FAKE_DORADO = """#!/usr/bin/env bash
set -euo pipefail
cmd="$1"; shift
if [[ "$cmd" == basecaller ]]; then
  echo "basecalling args: $*" >&2
  printf 'BAM'
elif [[ "$cmd" == demux ]]; then
  out=""
  while [[ $# -gt 0 ]]; do
    if [[ "$1" == --output-dir ]]; then out="$2"; shift; fi
    shift
  done
  mkdir -p "$out/run1/fastq_pass/barcode01" "$out/run1/fastq_pass/unclassified"
  printf '@r1\\nACGT\\n+\\nIIII\\n' > "$out/run1/fastq_pass/barcode01/a.fastq"
  printf '@r2\\nACGT\\n+\\nIIII\\n' > "$out/SQK-RBK114-96_barcode02.fastq"
  printf '@r3\\nACGT\\n+\\nIIII\\n' > "$out/run1/fastq_pass/unclassified/u.fastq"
fi
"""


class BasecallTests(unittest.TestCase):
    def test_arrange_handles_flat_and_nested_layouts(self):
        with tempfile.TemporaryDirectory() as tmp:
            demux = Path(tmp) / "demux"
            (demux / "x" / "fastq_pass" / "barcode05").mkdir(parents=True)
            (demux / "x" / "fastq_pass" / "barcode05" / "r.fastq").write_text("@r\nA\n+\nI\n")
            (demux / "KIT_barcode07.fastq").write_text("@r\nA\n+\nI\n")
            out = Path(tmp) / "out"
            counts = basecall.arrange_fastq(demux, out)
            self.assertEqual(counts, {"barcode05": 1, "barcode07": 1})
            files = sorted(p.relative_to(out).parts[:2] for p in out.rglob("*.fastq.gz"))
            self.assertEqual(files, [("fastq_pass", "barcode05"), ("fastq_pass", "barcode07")])
            first = next((out / "fastq_pass" / "barcode05").glob("*.fastq.gz"))
            self.assertEqual(gzip.decompress(first.read_bytes()), b"@r\nA\n+\nI\n")

    def test_end_to_end_with_fake_dorado(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dorado = root / "dorado"
            dorado.write_text(FAKE_DORADO)
            dorado.chmod(dorado.stat().st_mode | stat.S_IEXEC)
            pod5 = root / "data" / "run" / "pod5_pass"
            pod5.mkdir(parents=True)
            (pod5 / "a.pod5").write_bytes(b"x")
            model = root / "models" / "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
            model.mkdir(parents=True)
            with patch.dict(os.environ, {"ONT_DORADO": str(dorado), "ONT_DORADO_MODELS": str(root / "models")}):
                self.assertEqual([p.name for p in basecall.available_models()], [model.name])
                out = basecall.default_output_dir(pod5, str(model))
                self.assertTrue(out.name.startswith("basecall_sup_"))
                job = basecall.prepare_basecall_job(
                    root / "runs", basecall.BasecallSettings(str(pod5), str(out), str(model))
                )
                basecall.start_basecall_job(job)
                deadline = time.monotonic() + 20
                while basecall.basecall_job_state(job.job_dir) == "running" and time.monotonic() < deadline:
                    time.sleep(0.05)
                log = "\n".join(basecall.basecall_log_tail(job, 50))
                self.assertEqual(basecall.basecall_job_state(job.job_dir), "completed", log)
                self.assertIn("--kit-name SQK-RBK114-96", log)
                groups = sorted(p.name for p in (out / "fastq_pass").iterdir())
                self.assertEqual(groups, ["barcode01", "barcode02", "unclassified"])
                self.assertEqual([j.job_dir for j in basecall.basecall_jobs(root / "runs")], [job.job_dir])
                with self.assertRaises(basecall.BasecallError):  # output folder now non-empty
                    basecall.prepare_basecall_job(
                        root / "runs", basecall.BasecallSettings(str(pod5), str(out), str(model))
                    )


if __name__ == "__main__":
    unittest.main()
