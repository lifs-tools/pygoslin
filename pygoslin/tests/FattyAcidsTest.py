"""
MIT License

Copyright (c) 2020 Dominik Kopczynski   -   dominik.kopczynski {at} isas.de
                   Nils Hoffmann  -  nils.hoffmann {at} isas.de

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

"""
try:
    import pyximport
    pyximport.install(setup_args = {"script_args" : ["--force"]}, language_level = 3)
except:
    print("Warning: cython module is not installed, parsing performance will be lower since pure python code will be applied.")
"""

import pygoslin
from pygoslin.parser.Parser import *
from pygoslin.parser.ParserCommon import *
from pygoslin.domain.LipidLevel import *
from pygoslin.domain.Element import *
from pygoslin.domain.LipidExceptions import *

class LipidMapsTest(unittest.TestCase):
    
    def test_parser(self):
        lipidnames = []
        file_name = os.path.join("pygoslin", "tests", "fatty-acids-test.csv")
        with open(file_name, mode = "rt") as infile:
            lipidreader = csv.reader(infile, delimiter=',', quotechar='"')
            for row in lipidreader:
                lipidnames.append(row)
        
        
        lipid_parser = FattyAcidParser()
        shorthand_parser = ShorthandParser()
        formula_parser = SumFormulaParser()
    
        not_implemented, failed, failed_sum = 0, 0, 0
        with open("failed.txt", "wt") as output:
            for i, lipid_name in enumerate(lipidnames):
                name = lipid_name[3]
                if i and i % 100 == 0: print(i)
                
                if name.find("yn") >= 0 or name.find("furan") >= 0 or name[-3:] == "ane" or name[-3:] == "one" or name.find("phosphate") >= 0 or name.find("pyran") >= 0 or name[-5:] == "olide" or name[-4:] == "-one":
                    not_implemented += 1
                    continue                
                
                try:
                    lipid = lipid_parser.parse(name)
                except Exception as e:
                    failed += 1
                    output.write("'%s','%s',''\n" % (lipid_name[0], name.replace("'", "\\'")))
                    continue
                
                lipid_formula = lipid.get_sum_formula()
                formula = compute_sum_formula(formula_parser.parse(lipid_name[2]))
                
                if formula != lipid_formula:
                    print("%i, %s: %s != %s" % (i, lipid_name, formula, lipid_formula))
                    failed_sum += 1
                    exit()
                    
                if name.lower().find("cyano") >= 0:
                    continue
                
                lipid2 = shorthand_parser.parse(lipid.get_lipid_string())
                lipid_formula = lipid2.get_sum_formula()
                
                if formula != lipid_formula:
                    print("current, %i, %s: %s != %s / %s" % (i, lipid_name, formula, lipid_formula, lipid.get_lipid_string()))
                    failed_sum += 1
                    exit()
                    
                
                lipid2 = shorthand_parser.parse(lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES))
                lipid_formula = lipid2.get_sum_formula()
                
                if formula != lipid_formula:
                    print("molecular subspecies, %i, %s: %s != %s" % (i, lipid_name, formula, lipid_formula))
                    failed_sum += 1
                    exit()
                
                lipid2 = shorthand_parser.parse(lipid.get_lipid_string(LipidLevel.SPECIES))
                lipid_formula = lipid2.get_sum_formula()
                
                if formula != lipid_formula:
                    print("species, %i, %s: %s != %s" % (i, lipid_name, formula, lipid_formula))
                    failed_sum += 1
                    exit()
                
                
        print("In the test, %i of %i lipids can not be described by nomenclature" % (not_implemented, len(lipidnames)))
        print("In the test, %i of %i lipids failed" % (failed, len(lipidnames) - not_implemented))
