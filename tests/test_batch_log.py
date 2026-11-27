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


class BatchJobStateTests(unittest.TestCase):
    def _job(self, root: Path) -> BatchJob:
        import json
        (root / "batch_manifest.json").write_text(json.dumps({"experiment_name": "exp", "samples": [{}, {}]}))
        return BatchJob(root, root / "batch_config.yaml", root / "batch_manifest.json", "exp", 2)

    def test_state_reattach_and_stop(self):
        from ont_ui.batch import (active_batch_jobs, batch_job_state, follow_batch_job,
                                  latest_batch_job, load_batch_job, start_batch_job, stop_batch_job)
        with tempfile.TemporaryDirectory() as tmp:
            run_root = Path(tmp)
            job_dir = run_root / "job1"
            job_dir.mkdir()
            job = self._job(job_dir)
            repo = run_root / "repo"
            repo.mkdir()
            (repo / "run_pipeline.sh").write_text('echo "x Done: a"\nwhile [ ! -f "$1.go" ]; do sleep 0.05; done\necho "x Done: b"\n')
            self.assertEqual(batch_job_state(job_dir), "not_started")
            start_batch_job(repo, job)
            self.assertEqual(batch_job_state(job_dir), "running")
            self.assertEqual([j.job_dir for j in active_batch_jobs(run_root)], [job_dir])
            with self.assertRaises(BatchExecutionError):
                start_batch_job(repo, job)
            (job_dir / "batch_config.yaml.go").touch()
            reloaded = load_batch_job(job_dir)
            lines = []
            follow_batch_job(reloaded, lines.append)
            self.assertEqual(sum(" Done: " in line for line in lines), 2)
            self.assertEqual(batch_job_state(job_dir), "completed")
            self.assertEqual(latest_batch_job(run_root).job_dir, job_dir)

            (repo / "run_pipeline.sh").write_text("sleep 30\n")
            start_batch_job(repo, job)
            self.assertTrue(stop_batch_job(job))
            deadline = time.monotonic() + 5
            while batch_job_state(job_dir) == "running" and time.monotonic() < deadline:
                time.sleep(0.05)
            self.assertIn(batch_job_state(job_dir), {"stopped", "failed"})
