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

class LipidMapsTest(unittest.TestCase):
    
    def test_parser(self):
        lipidnames = []
        file_name = os.path.join("pygoslin", "tests", "fatty-acids-test.csv")
        with open(file_name, mode = "rt") as infile:
            lipidreader = csv.reader(infile, delimiter=',', quotechar='"')
            for row in lipidreader:
                lipidnames.append(row)
        
        
        lipid_parser = FattyAcidParser()
        formula_parser = SumFormulaParser()
    
        failed, failed_sum = 0, 0
        for i, lipid_name in enumerate(lipidnames):
            #if i and i % 100 == 0: print(i)
            try:
                lipid = lipid_parser.parse(lipid_name[3])
                lipid_formula = lipid.get_sum_formula()
                
                formula = compute_sum_formula(formula_parser.parse(lipid_name[2]))
                
                if formula != lipid_formula:
                    print("%s: %s != %s" % (lipid_name, formula, lipid_formula))
                    failed_sum += 1
                    exit()
                
                
                
            except Exception as e:
                failed += 1
                #if lipid_name[3].find("yn") < 0 and lipid_name[3].find("furan") < 0: 
                #print("'%s','%s',''" % (lipid_name[0], lipid_name[3].replace("'", "\\'")))
                
                
        print("In the test, %i of %i lipids failed" % (failed, len(lipidnames)))
        print("In the test, %i of %i lipid formulas failed" % (failed_sum, len(lipidnames) - failed))
