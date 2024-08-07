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


try:
    import pyximport
    pyximport.install(setup_args = {"script_args" : ["--force"]}, language_level = 3)
except:
    print("Warning: cython module is not installed, parsing performance will be lower since pure python code will be applied.")


import pygoslin
from pygoslin.parser.Parser import *
from pygoslin.domain.LipidLevel import *

class LipidMapsTest(unittest.TestCase):
    
    
    def test_sphingolipids(self):
        parser = LipidMapsParser()
        
        l = parser.parse("Cer(d18:1(8Z)/24:0)")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURE_DEFINED), "Cer 18:1(8);(OH)2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SN_POSITION), "Cer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), "Cer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "Cer 42:1;O2")
        self.assertEqual(l.get_sum_formula(), "C42H83NO3")
        
        l = parser.parse("GalCer(d18:1(5Z)/24:0)")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURE_DEFINED), "GalCer 18:1(5);OH/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SN_POSITION), "GalCer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), "GalCer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "GalCer 42:1;O2")
        self.assertEqual(l.get_sum_formula(), "C48H93NO8")
        
        l = parser.parse("LysoSM(d17:1(4E))")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURE_DEFINED), "LSM 17:1(4);OH")
        self.assertEqual(l.get_lipid_string(LipidLevel.SN_POSITION), "LSM 17:1;O2")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), "LSM 17:1;O2")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "LSM 17:1;O2")
        self.assertEqual(l.get_sum_formula(), "C22H47N2O5P")

        l = parser.parse("PE-Cer(d14:1(4E)/20:1(11Z))")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURE_DEFINED), "EPC 14:1(4);OH/20:1(11)")
        self.assertEqual(l.get_lipid_string(LipidLevel.SN_POSITION), "EPC 14:1;O2/20:1")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), "EPC 14:1;O2/20:1")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "EPC 34:2;O2")
        self.assertEqual(l.get_sum_formula(), "C36H71N2O6P")
        
        l = parser.parse("MIPC(t18:0/24:0)")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURE_DEFINED), "MIPC 18:0;(OH)2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SN_POSITION), "MIPC 18:0;O3/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), "MIPC 18:0;O3/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "MIPC 42:0;O3")
        self.assertEqual(l.get_sum_formula(), "C54H106NO17P")
        
        l = parser.parse("PE-Cer(d16:2(4E,6E)/22:1(13Z)(2OH))")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURE_DEFINED), "EPC 16:2(4,6);OH/22:1(13);OH")
        self.assertEqual(l.get_lipid_string(LipidLevel.SN_POSITION), "EPC 16:2;O2/22:1;O")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), "EPC 16:2;O2/22:1;O")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "EPC 38:3;O3")
        self.assertEqual(l.get_sum_formula(), "C40H77N2O7P")

    
    def test_parser(self):
        lipidnames = []
        file_name = os.path.join("pygoslin", "data", "goslin", "testfiles", "lipid-maps-test.csv")
        with open(file_name, mode = "rt", encoding= "utf-8") as infile:
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
            
                
