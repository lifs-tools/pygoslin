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


try:
    import pyximport
    pyximport.install(setup_args = {"script_args" : ["--force"]}, language_level = 3)
except:
    print("Warning: cython module is not installed, parsing performance will be lower since pure python code will be applied.")


import pygoslin
from pygoslin.parser.Parser import *
from pygoslin.domain.LipidLevel import *

class LipidMapsTest(unittest.TestCase):
    
    def test_parser(self):
        lipidnames = []
        file_name = os.path.join("pygoslin", "data", "goslin", "testfiles", "lipid-maps-test.csv")
        with open(file_name, mode = "rt") as infile:
            lipidreader = csv.reader(infile, delimiter=',', quotechar='"')
            for row in lipidreader:
                lipidnames.append(row)
        
        
        lipid_parser = LipidMapsParser()
        
        for i, lipid_name in enumerate(lipidnames):
            try:
                lipidname, correct_lipid_name = lipid_name
                lipid = lipid_parser.parse(lipidname)
                self.assertTrue(lipid != None)
                if correct_lipid_name != "Unsupported lipid" and len(correct_lipid_name) > 0:
                    self.assertTrue(correct_lipid_name == lipid.get_lipid_string(), "wrong: %s != %s (reference)" % (lipid.get_lipid_string(), correct_lipid_name))
                
            except Exception as e:
                if len(lipid_name[1]) > 0:
                    print("hier: '%s' -> %i" % (lipid_name[0], i))
                    exit()
            
                
