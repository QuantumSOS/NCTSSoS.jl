import unittest

from run import CASES, run_repeats, stock_diagnostics, validate_bound


class RepeatAndRecord(unittest.TestCase):
    def test_published_sdp_difference_is_not_an_acceptance_threshold(self):
        case = next(c for c in CASES if c.id == "square-16-rdm")
        self.assertEqual(case.paper_sdp_per_site, -0.674580)
        row = dict(method="reference", outcome="ok", termination_status="OPTIMAL",
                   stock_pfeas=1e-9, stock_dfeas=1e-9, solver_objective_relative_gap=0.0,
                   reported_objective_per_site=-0.70, primal_objective_per_site=-0.70)
        validate_bound(row, case, None)
        self.assertEqual(row["outcome"], "ok")
        self.assertAlmostEqual(row["energy_minus_published_sdp_per_site"], -0.025420)
        self.assertNotIn("signed_error_per_site", row)

    def test_stock_last_iteration_is_diagnostic_even_when_rejected(self):
        iteration = ("ITE PFEAS DFEAS GFEAS PRSTATUS POBJ DOBJ MU TIME\n"
                     "26 1e-9 1e-9 1e-17 1 -0.5 -0.5 1e-12 9.28\n"
                     "27 3.1e-8 9.4e-7 1.6e-17 1 -0.499999999 -0.5 1e-12 9.69\n"
                     "Optimizer terminated. Time: 10.14\n")
        for status, last, outcome in (
            ("SLOW_PROGRESS", iteration, "solver_rejected"),
            (None, iteration.replace("3.1e-8 9.4e-7", "2e-9 3e-9"), "ok"),
            (None, iteration.replace("-0.499999999 -0.5", "-0.499999999 NaN"), "solver_rejected"),
        ):
            with self.subTest(status=status, outcome=outcome):
                log = last + (f"termination status: {status}\nsolution status: FEASIBLE_POINT\n" if status else "")
                row = dict(method="reference", outcome="ok", reported_objective_per_site=-0.499999999)
                stock_diagnostics(row, log + "optimum = -0.499999999\n")
                validate_bound(row, CASES[0], -0.5)
                self.assertEqual(row["outcome"], outcome)
                self.assertEqual(row["primal_objective_per_site"], -0.499999999)
                self.assertEqual("energy_per_site" in row, outcome == "ok")
                if status:
                    self.assertEqual(row["termination_status"], status)
                    self.assertEqual(row["stock_pfeas"], 3.1e-8)
                    self.assertEqual(row["stock_dfeas"], 9.4e-7)
                    self.assertIn("solver_objective_relative_gap", row)
                elif outcome != "ok":
                    self.assertNotIn("dual_objective_per_site", row)
                    self.assertNotIn("solver_objective_relative_gap", row)

    def test_fresh_repeats_rotate_after_rejection_and_stop_after_resource_failure(self):
        methods = ["reference", "graph_distance", "local_products"]
        audit = dict(comparison_lane="exploratory", selected_rows=7,
                     retained_translation_orders=[4], retained_point_order=48,
                     row_span_match=False, verified_symmetry_match=True)
        for resource_stop in (None, "timeout", "oom"):
            with self.subTest(resource_stop=resource_stop):
                records = []

                def invoke(method, label, _case):
                    repeat = int(label.split("-")[2])
                    slow = repeat == 1 and method != "local_products"
                    return dict(method=method, label=label, prototype_row_counts=[7], translation_orders=[4],
                                point_order=48, stock_pfeas=1e-9, stock_dfeas=1e-9,
                                primal_status="FEASIBLE_POINT", dual_status="FEASIBLE_POINT",
                                termination_status="SLOW_PROGRESS" if slow else "OPTIMAL",
                                reported_objective_per_site=-0.5, primal_objective_per_site=-0.5,
                                dual_objective_per_site=-0.5, solver_objective_relative_gap=0.0,
                                outcome=resource_stop if resource_stop and repeat == 2 and method == "graph_distance" else "ok")

                self.assertTrue(run_repeats(CASES[0], -0.5, methods,
                                {m: audit for m in methods[1:]}, 3, invoke, records.append))
                attempts = [row for row in records if "method" in row]
                pairs = [row for row in records if "candidate" in row]
                self.assertEqual([r["method"] for r in attempts], methods + methods[1:] + methods[:1] +
                                 [methods[2], methods[0]] + ([] if resource_stop else [methods[1]]))
                self.assertEqual([r["replicate"] for r in pairs],
                                 [2, 3] if resource_stop else [2, 2, 3, 3])
                self.assertTrue(all(("energy_per_site" in r) == (r["outcome"] == "ok") for r in attempts))


if __name__ == "__main__":
    unittest.main()
