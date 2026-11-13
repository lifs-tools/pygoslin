import csv
import os
import unittest

from pygoslin.domain.LipidExceptions import LipidException
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.parser.Parser import ClassLevelParser, LipidParser

TESTFILES = os.path.join("pygoslin", "data", "goslin", "testfiles")


class ClassLevelTest(unittest.TestCase):
    """Names that are a lipid class or a lipid category only."""

    @classmethod
    def setUpClass(cls):
        cls.lipid_parser = LipidParser()
        cls.class_parser = ClassLevelParser()

    def test_class_level_names(self):
        with open(os.path.join(TESTFILES, "class-level-test.csv"), encoding="utf-8") as infile:
            rows = list(csv.reader(infile))[1:]
        self.assertGreater(len(rows), 0)
        for name, level, lipid_class, category in rows:
            for parser in (self.lipid_parser, self.class_parser):
                lipid = parser.parse(name)
                self.assertEqual(LipidLevel[level], lipid.lipid.info.level, name)
                self.assertEqual(lipid_class, lipid.get_lipid_string(LipidLevel.CLASS), name)
                self.assertEqual(category, lipid.get_lipid_string(LipidLevel.CATEGORY), name)
                expected = lipid_class if level == "CLASS" else category
                self.assertEqual(expected, lipid.get_lipid_string(), name)

    def test_invalid_class_level_names(self):
        with open(os.path.join(TESTFILES, "class-level-invalid.txt"), encoding="utf-8") as infile:
            names = [line.strip() for line in infile if line.strip()]
        self.assertGreater(len(names), 0)
        for name in names:
            self.assertIsNone(self.class_parser.parse(name, raise_error=False), name)
            with self.assertRaises(LipidException, msg=name):
                self.lipid_parser.parse(name)

    def test_species_names_are_not_class_level(self):
        lipid = self.lipid_parser.parse("PC 34:1")
        self.assertEqual(LipidLevel.SPECIES, lipid.lipid.info.level)

    def test_no_sum_formula_at_class_level(self):
        lipid = self.lipid_parser.parse("PC")
        with self.assertRaises(LipidException):
            lipid.get_sum_formula()


if __name__ == "__main__":
    unittest.main()
