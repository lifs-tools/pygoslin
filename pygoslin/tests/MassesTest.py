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
from pygoslin.domain.LipidLevel import *
from pygoslin.domain.Element import *

try:
    import pyximport
    pyximport.install(setup_args = {"script_args" : ["--force"]}, language_level = 3)
except:
    print("Warning: cython module is not installed, parsing performance will be lower since pure python code will be applied.")

parser = GoslinParser()

class TestFormulas(unittest.TestCase):

    def test_masses(self):
        global parser

        with open('pygoslin/tests/lipid-masses.csv', newline='') as csvfile:
            lipidreader = csv.reader(csvfile, delimiter=',', quotechar='\"')
            for i, row in enumerate(lipidreader):
                if i == 0: continue
                
                lipid_name = row[1] + row[3]
                
                try:
                    lipid = parser.parse(lipid_name)
                except LipidException as e:
                    print(row[0], e)
                    assert False
                 
                assert lipid.get_lipid_string(LipidLevel.CLASS) == row[0]
                assert compute_sum_formula(lipid.lipid.get_elements()) == row[2]
                assert abs(lipid.get_mass() - float(row[4])) < 0.001
                assert lipid.adduct.get_charge() == int(row[5])


