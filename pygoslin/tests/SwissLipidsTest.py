import unittest
import os


import pyximport
pyximport.install(setup_args = {"script_args" : ["--force"]}, language_level = 3)

import pygoslin
from pygoslin.parser.Parser import SwissLipidsParser
from pygoslin.parser.GoslinParserEventHandler import GoslinParserEventHandler
from pygoslin.parser.GoslinFragmentParserEventHandler import GoslinFragmentParserEventHandler
from pygoslin.parser.LipidMapsParserEventHandler import LipidMapsParserEventHandler
from pygoslin.domain.LipidLevel import LipidLevel
from random import randint

class ParserTest(unittest.TestCase):
    PARSER_QUOTE = '\''
    
    
    def test_parser(self):
        lipidnames = []
        with open("pygoslin/tests/swiss-lipids-test.csv", mode = "rt") as infile:
            for line in infile:
                line = line.strip().strip(" ")
                if len(line) > 0: lipidnames.append(line)
        
        
        lipid_parser = SwissLipidsParser()
        
        for i, lipid_name in enumerate(lipidnames):
            lipid = lipid_parser.parse(lipid_name)
            if lipid == None: print("hier: '%s' -> %i" % (lipid_name, i))
            assert lipid != None