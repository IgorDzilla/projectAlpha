"""Проверки чтения вывода solver-cli и построения графиков без окна."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import plot_toro as plot


class PlotToroTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="plot-toro-test-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.case_dir = self.root / "toro_1"
        self.case_dir.mkdir()
        self.columns = ["time", "cell", "x", "density", "momentum", "energy", "velocity", "pressure"]
        self.rows = []
        for t in (0.0, 0.07, 0.2):
            for cell in range(2):
                rho, u, p = 1 + 0.2 * cell + 0.05 * t, 0.1 * cell, 1 + 0.2 * t
                self.rows.append([t, cell, 0.25 + 0.5 * cell, rho, rho*u,
                                  p/(1.4-1) + 0.5*rho*u*u, u, p])
        self.meta = dict(model="euler_1d", status="complete", cells="2", snapshots="3",
                         gamma="1.4", x_min="0", x_max="1", t_start="0", t_end="0.2", final_time="0.2")
        self.write()

    def write(self, columns=None, rows=None):
        (self.case_dir / "run.txt").write_text(
            "".join(f"{k}={v}\n" for k, v in self.meta.items()), encoding="utf-8")
        with (self.case_dir / "solution.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(columns if columns is not None else self.columns)
            writer.writerows(rows if rows is not None else self.rows)

    def test_current_format_and_shuffled_rows(self):
        expected = plot.read_case(self.case_dir)
        self.write(rows=list(reversed(self.rows)))
        actual = plot.read_case(self.case_dir)
        np.testing.assert_array_equal(actual.times, [0, 0.07, 0.2])
        np.testing.assert_array_equal(actual.values["density"], expected.values["density"])
        self.assertEqual(actual.values["density"].shape, (3, 2))

    def test_derive_primitive_fields(self):
        expected = plot.read_case(self.case_dir)
        self.write(columns=self.columns[:6], rows=[row[:6] for row in self.rows])
        actual = plot.read_case(self.case_dir)
        for key in ("velocity", "pressure", "internal_energy", "sound_speed", "entropy"):
            np.testing.assert_allclose(actual.values[key], expected.values[key], rtol=1e-12)

    def test_partial_fields_and_render(self):
        self.write(columns=["time", "cell", "x", "pressure"],
                   rows=[row[:3] + [row[-1]] for row in self.rows])
        case = plot.read_case(self.case_dir)
        self.assertEqual(set(case.values), {"pressure"})
        out = self.root / "plots"
        plot.render([case], out)
        self.assertTrue((out / "toro_1/profiles.png").is_file())
        self.assertTrue((out / "toro_1/xt.png").is_file())
        self.assertFalse((out / "toro_1/pu.png").exists())
        self.assertIn("пропущена", (out / "REPORT.md").read_text())

    def test_single_cell_single_time(self):
        self.meta.update(cells="1", snapshots="1", t_end="0", final_time="0")
        row = self.rows[0].copy()
        row[2] = 0.5
        self.write(rows=[row])
        case = plot.read_case(self.case_dir)
        out = self.root / "single-plots"
        plot.render([case], out)
        self.assertEqual(len(list((out / "toro_1").glob("*.png"))), 3)

    def test_nonuniform_times(self):
        edges = plot.time_edges(np.array([1.0, 1.1, 1.7, 2.0]))
        np.testing.assert_allclose(edges, [1.0, 1.05, 1.4, 1.85, 2.0])

    def test_bad_layers(self):
        for change in (lambda r: r.pop(), lambda r: r[1].__setitem__(1, 0),
                       lambda r: r[0].__setitem__(1, 0.5),
                       lambda r: r[0].__setitem__(0, 0.001),
                       lambda r: r[0].__setitem__(2, 0.4)):
            rows = [r.copy() for r in self.rows]
            change(rows)
            self.write(rows=rows)
            with self.assertRaises(ValueError):
                plot.read_case(self.case_dir)

    def test_bad_physics_and_inconsistent_fields(self):
        for column, value in ((3, -1), (3, float("nan")), (7, 0), (7, 100), (6, 10)):
            rows = [r.copy() for r in self.rows]
            rows[0][column] = value
            self.write(rows=rows)
            with self.assertRaises(ValueError):
                plot.read_case(self.case_dir)

    def test_failed_run_cli_creates_no_output(self):
        self.meta["status"] = "failed"
        self.write()
        output = self.root / "output"
        process = subprocess.run([sys.executable, plot.__file__, str(self.root),
                                  "--output-dir", str(output)], capture_output=True, text=True)
        self.assertEqual(process.returncode, 1)
        self.assertIn("status=complete", process.stderr)
        self.assertNotIn("Traceback", process.stderr)
        self.assertFalse(output.exists())

    def test_discovery(self):
        self.assertEqual(plot.discover_cases(self.root), [self.case_dir])
        self.assertEqual(plot.discover_cases(self.case_dir), [self.case_dir])
        with self.assertRaises(ValueError):
            plot.discover_cases(self.root / "missing")

    def exact_fixture(self):
        meta = dict(test="1", passed="1", gamma="1.4", x_min="0", x_max="1", x0="0.5",
                    t_start="0", tout="0.2", p_star="0.30313", u_star="0.92745",
                    rho_star_left="0.42632", rho_star_right="0.26557", guess="two-shocks",
                    initial_pressure="0.29", iterations="3", contact="0.92745",
                    left_type="rarefaction", right_type="shock", left_head="-1.1832",
                    left_tail="-0.0703", right_head="1.7522", right_tail="1.7522",
                    left_rho="1", left_u="0", left_p="1", right_rho="0.125", right_u="0", right_p="0.1")
        def write(name, header, rows):
            with (self.case_dir / name).open("w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(header)
                writer.writerows(rows)
        write("summary.csv", ["key", "value"], meta.items())
        values = [(0, 1, 0, 1), (.5, .42632, .92745, .30313), (1, .125, 0, .1)]
        profile = [[x, r, u, p, p/(.4*r), np.sqrt(1.4*p/r), p/r**1.4] for x, r, u, p in values]
        write("profile.csv", ["x", "density", "velocity", "pressure", "internal_energy", "sound_speed", "entropy"], profile)
        write("density_xt.csv", ["time", "cell", "x", "density"],
              [[t, i, row[0], row[1]] for t in (0, .2) for i, row in enumerate(values)])
        write("wave_curves.csv", ["p", "u_left", "u_right"], [[.1, 1.6, 0], [.30313, .92745, .92745], [1, 0, 2.3]])
        write("errors.csv", plot.ERROR_COLUMNS,
              [[k, meta[k], meta[k], 0, 0, 1e-5, 1] for k in plot.QUANTITIES])

    def test_exact_profiles_curves_and_report(self):
        self.exact_fixture()
        case = plot.read_case(self.case_dir)
        self.assertTrue(case.exact)
        self.assertEqual(case.xt_rho.shape, (2, 3))
        np.testing.assert_allclose(case.curves["u_left"][1], float(case.meta["u_star"]))
        out = self.root / "exact-plots"
        plot.render([case], out)
        self.assertEqual(len(list((out / "toro_1").glob("*.pdf"))), 3)
        self.assertTrue((out / "toro_report.pdf").is_file())
        report = (out / "REPORT.md").read_text()
        for text in ("GUESSP", "ρ*L", "ρ*R", "Голова", "Контакт"):
            self.assertIn(text, report)

    def test_failed_exact_case_is_not_reported_as_passed(self):
        self.exact_fixture()
        path = self.case_dir / "summary.csv"
        path.write_text(path.read_text().replace("passed,1", "passed,0"))
        with self.assertRaises(ValueError):
            plot.read_case(self.case_dir)
        self.exact_fixture()
        path = self.case_dir / "profile.csv"
        path.write_text(path.read_text().replace("0,1,0,1,", "0,-1,0,1,"))
        with self.assertRaises(ValueError):
            plot.read_case(self.case_dir)


if __name__ == "__main__":
    unittest.main()
