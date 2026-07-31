import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts" / "mock_benchmark"))
import metrics


class MockMetricsTest(unittest.TestCase):
    def setUp(self):
        self.truth = pd.DataFrame([
            {"Assembly": "a", "ani_level": 98.0, "read_count": 10, "family": "f__A", "genus": "g__A", "species": "s__A", "subspecies": "t__A"},
            {"Assembly": "b", "ani_level": 85.0, "read_count": 20, "family": "f__B", "genus": "g__B", "species": "s__B", "subspecies": "t__B"},
        ])
        self.prediction = pd.DataFrame([
            {"count": 10, "family": "f__A", "genus": "g__A", "species": "s__A", "subspecies": "t__A"},
            {"count": 20, "family": "f__B", "genus": "g__B", "species": "s__B", "subspecies": "t__B"},
        ])

    def test_ani_caps_include_subspecies_and_genus(self):
        self.assertEqual(metrics.DEFAULT_ANI_TIERS[98.0], "subspecies")
        self.assertEqual(metrics.DEFAULT_ANI_TIERS[95.0], "species")
        self.assertEqual(metrics.DEFAULT_ANI_TIERS[85.0], "genus")

    def test_capped_vectors_roll_to_each_truth_family_cap(self):
        truth, prediction, _ = metrics.capped_vectors(self.truth, self.prediction, metrics.DEFAULT_ANI_TIERS)
        self.assertEqual(truth.to_dict(), {"subspecies|t__A": 10.0, "genus|g__B": 20.0})
        self.assertEqual(prediction.to_dict(), {"subspecies|t__A": 10.0, "genus|g__B": 20.0})

    def test_capped_metrics_preserve_fixed_denominator_recovery(self):
        rows, capped_detection, _ = metrics.compute_community_and_detection(self.truth, self.prediction, "esviritu", 30)
        capped = next(row for row in rows if row["truth_mode"] == "capped")
        self.assertEqual(capped_detection, {"TP": 2, "FP": 0, "FN": 0, "precision": 1.0, "recall": 1.0, "F1": 1.0})
        self.assertAlmostEqual(capped["truth_total_rpm"], 1_000_000)
        self.assertAlmostEqual(capped["pred_total_rpm"], 1_000_000)
        self.assertEqual(capped["bray_curtis_rpm"], 0.0)

    def test_capped_scoring_rejects_shared_truth_families(self):
        duplicated = pd.concat([self.truth, self.truth.iloc[[0]]], ignore_index=True)
        with self.assertRaisesRegex(ValueError, "one truth genome per family"):
            metrics.capped_vectors(duplicated, self.prediction, metrics.DEFAULT_ANI_TIERS)


if __name__ == "__main__":
    unittest.main()
