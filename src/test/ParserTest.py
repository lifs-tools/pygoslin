import unittest
import os

from parser.Parser import Parser, Bitfield
from parser.GoslinParserEventHandler import GoslinParserEventHandler
from parser.GoslinFragmentParserEventHandler import GoslinFragmentParserEventHandler
from parser.LipidMapsParserEventHandler import LipidMapsParserEventHandler
from random import randint

class ParserTest(unittest.TestCase):
    PARSER_QUOTE = '\''
    
    
    def test_bitfield(self):
        n = randint(5000, 20000)
        m = randint(500, 2000)
        s = set([randint(1, n - 2) for x in range(m)])
        b = Bitfield(n)
        for i in s:
            b.set(i)
        
        i = 0
        for bb, st in zip(b.get_bit_positions(), sorted(s)):
            assert bb == st
            i += 1
    
        assert i == len(s)
        
    
    
    
    def test_tree_node(self):
        goslin_parser_event_handler = GoslinParserEventHandler()
        goslin_parser = Parser(goslin_parser_event_handler, os.path.abspath(".") + "/test/Goslin.g4", ParserTest.PARSER_QUOTE)
        
        
        lipid_name = "PE 16:1"
        goslin_parser.parse(lipid_name)
        assert goslin_parser.word_in_grammar
        goslin_parser.raise_events()
        
        assert lipid_name == goslin_parser.parse_tree.get_text()
        
    
    def test_parser_read(self):
        lipidnames = []
        with open(os.path.abspath(".") + "/test/lipidnames.txt", mode = "rt") as infile:
            for line in infile:
                line = line.strip().strip(" ")
                if len(line) < 2: continue
                if line[0] == "#": continue
                lipidnames.append(Parser.split_string(line, ",", "\"")[0].strip("\""))
        
        goslin_parser_event_handler = GoslinParserEventHandler()
        goslin_parser = Parser(goslin_parser_event_handler, os.path.abspath(".") + "/test/Goslin.g4", ParserTest.PARSER_QUOTE)
        
        
        goslin_fragment_parser_event_handler = GoslinFragmentParserEventHandler()
        fragment_parser = Parser(goslin_fragment_parser_event_handler, os.path.abspath(".") + "/test/Goslin-Fragments.g4", ParserTest.PARSER_QUOTE)
        
        
        lipid_maps_parser_event_handler = LipidMapsParserEventHandler()
        lipid_maps_parser = Parser(lipid_maps_parser_event_handler, os.path.abspath(".") + "/test/LipidMaps.g4", ParserTest.PARSER_QUOTE)
        
        parsers = [[goslin_parser, goslin_parser_event_handler],
                   [fragment_parser, goslin_fragment_parser_event_handler],
                   [lipid_maps_parser, lipid_maps_parser_event_handler]]
        
        
        for lipidname in lipidnames:
            print("check: %s" % lipidname)
            found = False
            for pp in parsers:
                prs, handler = pp
                prs.parse(lipidname)
                if prs.word_in_grammar:
                    found = True
                    break
                    prs.raise_events()
            assert found