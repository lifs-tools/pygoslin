import unittest
import os


try:
    import pyximport
    pyximport.install(setup_args = {"script_args" : ["--force"]}, language_level = 3)
except:
    print("Warning: cython module is not installed, parsing performance will be lower since pure python code will be applied.")

import pygoslin
from pygoslin.parser.Parser import GoslinParser

class ParserTest(unittest.TestCase):
    PARSER_QUOTE = '\''
    
    
    
    def test_parser(self):
        lipidnames = []
        with open("pygoslin/tests/goslin-test.csv", mode = "rt") as infile:
            for line in infile:
                line = line.strip().strip(" ")
                if len(line) > 0: lipidnames.append(line)
        
        
        lipid_parser = GoslinParser()
        
        for i, lipid_name in enumerate(lipidnames):
            lipid = lipid_parser.parse(lipid_name)
            if lipid == None: print("hier: '%s' -> %i" % (lipid_name, i))
            assert lipid != None