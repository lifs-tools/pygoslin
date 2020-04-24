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
from pygoslin.parser.Parser import *
from pygoslin.domain.LipidExceptions import *


try:
    import pyximport
    pyximport.install(setup_args = {"script_args" : ["--force"]}, language_level = 3)
except:
    print("Warning: cython module is not installed, parsing performance will be lower since pure python code will be applied.")

swiss_lipids_parser = SwissLipidsParser()
lipids_maps_parser = LipidMapsParser()

class TestFormulas(unittest.TestCase):

    def test_formulas_swiss_lipids(self):
        global swiss_lipids_parser

        with open('pygoslin/tests/formulas-swiss-lipids.csv', newline='') as csvfile:
            lipidreader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            
            for row in lipidreader:
                
                try:
                    lipid = swiss_lipids_parser.parse(row[0])
                except LipidException as e:
                    print(row[0], e)
                    assert False
                    
                try:
                    formula = lipid.get_sum_formula()
                    if len(formula) == 0: continue
                    if formula != row[1]:
                        print("assert for %s: %s != %s" % (row[0], formula, row[1]))
                        assert False
                except Exception as e:
                    print(row[0], e)
                    assert False


    def test_formulas_lipid_maps(self):
        global lipids_maps_parser

        with open('pygoslin/tests/formulas-lipid-maps.csv', newline='') as csvfile:
            lipidreader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            
            for row in lipidreader:
                
                try:
                    lipid = lipids_maps_parser.parse(row[0])
                except LipidException as e:
                    print(row[0], e)
                    assert False
                    
                try:
                    formula = lipid.get_sum_formula()
                    if len(formula) == 0: continue
                    if formula != row[1]:
                        print("assert for %s: %s != %s" % (row[0], formula, row[1]))
                        assert False
                except Exception as e:
                    print(row[0], e)
                    assert False
                    
