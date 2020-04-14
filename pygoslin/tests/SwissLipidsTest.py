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


try:
    import pyximport
    pyximport.install(setup_args = {"script_args" : ["--force"]}, language_level = 3)
except:
    print("Warning: cython module is not installed, parsing performance will be lower since pure python code will be applied.")

import pygoslin
from pygoslin.parser.Parser import SwissLipidsParser
from pygoslin.domain.LipidLevel import LipidLevel

class SwissLipidsParserTest(unittest.TestCase):
    PARSER_QUOTE = '\''
    
    
    
    def test_swiss_lipids_parser(self):
        lipidnames = []
        with open("pygoslin/tests/swiss-lipids-test.csv", mode = "rt") as infile:
            for line in infile:
                line = line.strip().strip(" ")
                if len(line) > 0: lipidnames.append(line)
        
        
        lipid_parser = SwissLipidsParser()
        
        for i, lipid_name in enumerate(lipidnames):
            try:
                lipid = lipid_parser.parse(lipid_name)
                lipid_class = lipid.get_lipid_string(LipidLevel.CLASS)
                assert lipid_class not in {"Undefined lipid class", "Undefined", "UNDEFINED"}
            except Exception as e:
                print("hier: '%s' -> %i" % (lipid_name, i))
                print(e)
                assert (False)
