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


try:
    import pyximport
    pyximport.install(setup_args = {"script_args" : ["--force"]}, language_level = 3)
except:
    print("Warning: cython module is not installed, parsing performance will be lower since pure python code will be applied.")

import pygoslin
from pygoslin.parser.Parser import GoslinParser
from pygoslin.domain.LipidLevel import *

class GoslinTest(unittest.TestCase):
    PARSER_QUOTE = '\''
    
    
    def test_sphingolipids(self):
        parser = GoslinParser()
        
        l = parser.parse("Cer 18:1(8Z);2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURE_DEFINED), "Cer 18:1(8);(OH)2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SN_POSITION), "Cer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), "Cer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "Cer 42:1;O2")
        self.assertEqual(l.get_sum_formula(), "C42H83NO3")
        
        l = parser.parse("HexCer 18:1(5Z);2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURE_DEFINED), "HexCer 18:1(5);OH/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SN_POSITION), "HexCer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), "HexCer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "HexCer 42:1;O2")
        self.assertEqual(l.get_sum_formula(), "C48H93NO8")
        
        l = parser.parse("LSM 17:1(4E);2")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURE_DEFINED), "LSM 17:1(4);OH")
        self.assertEqual(l.get_lipid_string(LipidLevel.SN_POSITION), "LSM 17:1;O2")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), "LSM 17:1;O2")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "LSM 17:1;O2")
        self.assertEqual(l.get_sum_formula(), "C22H47N2O5P")
        
        l = parser.parse("LCB 18:1(4E);2")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURE_DEFINED), "SPB 18:1(4);(OH)2")
        self.assertEqual(l.get_lipid_string(LipidLevel.SN_POSITION), "SPB 18:1;O2")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), "SPB 18:1;O2")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "SPB 18:1;O2")
        self.assertEqual(l.get_sum_formula(), "C18H37NO2")

        l = parser.parse("EPC 14:1(4E);2/20:1(11Z)")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURE_DEFINED), "EPC 14:1(4);OH/20:1(11)")
        self.assertEqual(l.get_lipid_string(LipidLevel.SN_POSITION), "EPC 14:1;O2/20:1")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), "EPC 14:1;O2/20:1")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "EPC 34:2;O2")
        self.assertEqual(l.get_sum_formula(), "C36H71N2O6P")
        
        l = parser.parse("MIPC 18:0;3/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURE_DEFINED), "MIPC 18:0;(OH)2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SN_POSITION), "MIPC 18:0;O3/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), "MIPC 18:0;O3/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "MIPC 42:0;O3")
        self.assertEqual(l.get_sum_formula(), "C54H106NO17P")
        
        l = parser.parse("EPC 16:2(4E,6E);2/22:1(13Z);1")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURE_DEFINED), "EPC 16:2(4,6);OH/22:1(13);OH")
        self.assertEqual(l.get_lipid_string(LipidLevel.SN_POSITION), "EPC 16:2;O2/22:1;O")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), "EPC 16:2;O2/22:1;O")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "EPC 38:3;O3")
        self.assertEqual(l.get_sum_formula(), "C40H77N2O7P")
    
    
        l = parser.parse("13-oxoODE")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), "FA 18:3;O")
        self.assertEqual(l.get_sum_formula(), "C18H30O3")
        
        
    
    
    def test_parser(self):
        lipidnames = []
        file_name = os.path.join("pygoslin", "data", "goslin", "testfiles", "goslin-test.csv")
        with open(file_name, mode = "rt", encoding= "utf-8") as infile:
            for line in infile:
                line = line.strip().strip(" ")
                if len(line) > 0: lipidnames.append(line)
        
        
        lipid_parser = GoslinParser()
        
        for i, lipid_name in enumerate(lipidnames):
            try:
                lipid = lipid_parser.parse(lipid_name)
                lipid_class = lipid.get_lipid_string(LipidLevel.CLASS)
                assert lipid_class not in {"Undefined lipid class", "Undefined", "UNDEFINED"}
                
            except Exception as e:
                print("hier: '%s' -> %i" % (lipid_name, i))
                print(e)
                assert (False)
