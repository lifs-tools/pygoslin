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

class ParserTest(unittest.TestCase):
    
    def test_parser(self):
        lipidnames = []
        with open("pygoslin/tests/lipid-maps-test.csv", mode = "rt") as infile:
            lipidreader = csv.reader(infile, delimiter=',', quotechar='"')
            for row in lipidreader:
                lipidnames.append(row)
        
        
        lipid_parser = LipidMapsParser()
        
        for i, lipid_name in enumerate(lipidnames):
            lipid = lipid_parser.parse(lipid_name[0])
            if lipid == None and len(lipid_name[1]) > 0:
                print("hier: '%s' -> %i" % (lipid_name[0], i))
                exit()
                
            if lipid == None: continue
            
            lipid_class = lipid.get_lipid_string(LipidLevel.CLASS)
            assert lipid_class not in {"Undefined lipid class", "Undefined", "UNDEFINED"}
            
                