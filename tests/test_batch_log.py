import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

from ont_ui.batch import BatchExecutionError, BatchJob, run_batch_job


class BatchLogTests(unittest.TestCase):
    def test_complete_lines_and_failure_tail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            job = BatchJob(root, root / "cfg", root / "manifest", "test", 1)
            (root / "run_pipeline.sh").write_text(
                "printf partial\nsleep 0.25\nprintf ' line\\nerror detail' >&2\nexit 7\n"
            )
            lines = []
            with self.assertRaisesRegex(BatchExecutionError, "code 7") as error:
                run_batch_job(root, job, lines.append)
            self.assertEqual(lines, ["partial line", "error detail"])
            self.assertIn("error detail", str(error.exception))
            self.assertEqual((root / "pipeline.log").read_text(), "partial line\nerror detail")

    def test_pipeline_finishes_after_monitor_exits(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "run_pipeline.sh").write_text(
                'echo START\nwhile [ ! -f "$1" ]; do sleep 0.05; done\n'
                "head -c 262144 /dev/zero | tr '\\\\0' x\necho FINISHED\n"
            )
            monitor = subprocess.Popen([
                sys.executable, "-c",
                "from pathlib import Path; from ont_ui.batch import BatchJob, run_batch_job; "
                "import sys; r=Path(sys.argv[1]); "
                "run_batch_job(r, BatchJob(r,r/'ready',r/'manifest','test',1))", tmp
            ])
            log = root / "pipeline.log"
            def wait_for(token):
                deadline = time.monotonic() + 8
                while time.monotonic() < deadline:
                    if log.exists() and token in log.read_bytes():
                        return True
                    time.sleep(0.03)
                return False
            try:
                self.assertTrue(wait_for(b"START"))
                monitor.terminate()
                monitor.wait(timeout=5)
            finally:
                (root / "ready").touch()
                if monitor.poll() is None:
                    monitor.kill()
                    monitor.wait()
            self.assertTrue(wait_for(b"FINISHED"))
            self.assertGreater(log.stat().st_size, 262144)
