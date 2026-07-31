import unittest

from scripts.control_read_units import (
    base_tool,
    negative_count_to_individual_reads,
    positive_count_to_read_pairs,
)


class ControlReadUnitTests(unittest.TestCase):
    def test_base_tool_handles_parameterized_and_display_names(self):
        self.assertEqual(base_tool("esviritu_sr"), "esviritu")
        self.assertEqual(base_tool("EsViritu (sense)"), "esviritu")
        self.assertEqual(base_tool("kraken2_conf0.5_mhg2"), "kraken2")

    def test_positive_esviritu_individual_reads_become_pairs(self):
        self.assertEqual(positive_count_to_read_pairs("esviritu_sr", 20), 10)

    def test_positive_paired_tool_counts_remain_pairs(self):
        for tool in ("kraken2_conf0_mhg2", "centrifuger_mhl100_ms0", "metabuli_prec1"):
            self.assertEqual(positive_count_to_read_pairs(tool, 10), 10)

    def test_positive_sylph_has_no_read_pair_count(self):
        self.assertIsNone(positive_count_to_read_pairs("sylph_mnk10", 1))

    def test_negative_paired_short_read_counts_become_individual_reads(self):
        for tool in ("Kraken2", "Centrifuger", "Metabuli"):
            self.assertEqual(
                negative_count_to_individual_reads(tool, 10, "illumina_short"), 20
            )

    def test_negative_esviritu_is_already_individual_reads(self):
        self.assertEqual(
            negative_count_to_individual_reads("EsViritu (sr)", 10, "illumina_short"),
            10,
        )

    def test_negative_single_read_counts_are_not_doubled(self):
        self.assertEqual(
            negative_count_to_individual_reads("Centrifuger", 10, "ONT_long"), 10
        )

    def test_negative_sylph_has_no_read_count(self):
        self.assertIsNone(
            negative_count_to_individual_reads("Sylph", 3, "illumina_short")
        )


if __name__ == "__main__":
    unittest.main()
