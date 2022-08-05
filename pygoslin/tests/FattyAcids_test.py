"""
MIT License

Copyright (c) 2020 Dominik Kopczynski   -   dominik.kopczynski {at} univie.ac.at
                   Nils Hoffmann  -  nils.hoffmann {at} cebitec.uni-bielefeld.de

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


import unittest
import os
import csv


import pygoslin
from pygoslin.parser.Parser import *
from pygoslin.parser.ParserCommon import *
from pygoslin.domain.LipidLevel import *
from pygoslin.domain.Element import compute_sum_formula

class LipidMapsTest(unittest.TestCase):
    
    def test_parser(self):
        lipid_data = []
        file_name = os.path.join("pygoslin", "data", "goslin", "testfiles", "fatty-acids-test.csv")
        with open(file_name, mode = "rt", encoding= "utf-8") as infile:
            lipidreader = csv.reader(infile, delimiter=',', quotechar='"')
            for row in lipidreader:
                lipid_data.append(row)
        
        
        parser = FattyAcidParser()
        shorthand_parser = ShorthandParser()
        formula_parser = SumFormulaParser()
        for lipid_row in lipid_data:
            
            lmid, lipid_name, formula, expected_lipid_name = lipid_row
            formula = compute_sum_formula(formula_parser.parse(formula))
            lipid = parser.parse(lipid_name)
            
            
            self.assertEqual(expected_lipid_name, lipid.get_lipid_string(), "%s: %s != %s (computed, fa)" % (lmid, expected_lipid_name, lipid.get_lipid_string()))
            
            lipid_formula = lipid.get_sum_formula()
            
            self.assertEqual(formula, lipid_formula, "formula %s: %s != %s (computed, fa)" % (lmid, formula, lipid_formula))
                
            if lipid_name.lower().find("cyano") >= 0: continue
            
            lipid2 = shorthand_parser.parse(lipid.get_lipid_string())
            lipid_formula = lipid2.get_sum_formula()
            
            self.assertEqual(formula, lipid_formula, "lipid %s: %s != %s (computed, shorthand)" % (lmid, formula, lipid_formula))
                
            lipid2 = shorthand_parser.parse(lipid.get_lipid_string(LipidLevel.MOLECULAR_SPECIES))
            lipid_formula = lipid2.get_sum_formula()
            
            self.assertEqual(formula, lipid_formula, "molecular lipid '%s': %s != %s (computed, shorthand)" % (lmid, formula, lipid_formula))
                
            
            lipid2 = shorthand_parser.parse(lipid.get_lipid_string(LipidLevel.SPECIES))
            lipid_formula = lipid2.get_sum_formula()
            
            self.assertEqual(formula, lipid_formula, "species lipid '%s': %s != %s (computed, shorthand)" % (lmid, formula, lipid_formula))
                    
