import csv
import os
import unittest

from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.parser.Parser import LDAParser, LipidMapsParser, LipidParser

TESTFILES = os.path.join("pygoslin", "data", "goslin", "testfiles")


class LDATest(unittest.TestCase):
    """Sphingolipid names written by LipidDataAnalyzer (LDA.g4)."""

    @classmethod
    def setUpClass(cls):
        cls.parser = LDAParser()
        cls.lipid_parser = LipidParser()

    def test_lda_names(self):
        with open(os.path.join(TESTFILES, "lda-test.csv"), encoding="utf-8") as infile:
            rows = list(csv.reader(infile))[1:]
        self.assertGreater(len(rows), 0)
        for name, goslin_name, level, formula in rows:
            lipid = self.parser.parse(name)
            self.assertEqual(goslin_name, lipid.get_lipid_string(), name)
            self.assertEqual(LipidLevel[level], lipid.lipid.info.level, name)
            self.assertEqual(formula, lipid.get_sum_formula(), name)
            # the multi grammar parser reads the name the same way
            lipid = self.lipid_parser.parse(name)
            self.assertEqual(goslin_name, lipid.get_lipid_string(), name)
            self.assertEqual(formula, lipid.get_sum_formula(), name)

    def test_invalid_lda_names(self):
        with open(os.path.join(TESTFILES, "lda-invalid.txt"), encoding="utf-8") as infile:
            names = [line.strip() for line in infile if line.strip()]
        self.assertGreater(len(names), 0)
        for name in names:
            self.assertIsNone(self.parser.parse(name, raise_error=False), name)

    def test_lipid_maps_m_lcb_with_head_group_at_c1(self):
        # one hydroxy group on a long chain base that carries the head group at C1
        parser = LipidMapsParser()
        self.assertEqual("C47H85N2O5P", parser.parse("SM m42:6").get_sum_formula())
        self.assertEqual("C35H73N2O5P", parser.parse("SM(m18:0/12:0)").get_sum_formula())
        self.assertEqual("C20H42NO4P", parser.parse("S1P m20:1").get_sum_formula())
        # ceramides keep the 1-deoxy (3-OH) reading
        self.assertEqual("Cer 18:0;3OH/12:0", parser.parse("Cer(m18:0/12:0)").get_lipid_string())


if __name__ == "__main__":
    unittest.main()
