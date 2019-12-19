import unittest
import os

from pygoslin.parser.Parser import Parser, Bitfield
from pygoslin.parser.GoslinParserEventHandler import GoslinParserEventHandler
from pygoslin.parser.GoslinFragmentParserEventHandler import GoslinFragmentParserEventHandler
from pygoslin.parser.LipidMapsParserEventHandler import LipidMapsParserEventHandler
from pygoslin.domain.LipidLevel import LipidLevel
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
        goslin_parser = Parser(goslin_parser_event_handler, "data/goslin/Goslin.g4", ParserTest.PARSER_QUOTE)
        
        lipid_name = "PE 16:1-12:0"
        goslin_parser.parse(lipid_name)
        assert goslin_parser.word_in_grammar
        goslin_parser.raise_events()
        
        assert lipid_name == goslin_parser.parse_tree.get_text()
        
        
        
    
    
    @unittest.expectedFailure
    def test_lipid_fragment_fail(self):
    
        goslin_parser_event_handler = GoslinParserEventHandler()
        goslin_parser = Parser(goslin_parser_event_handler, "data/goslin/Goslin.g4", ParserTest.PARSER_QUOTE)
        
        
        lipid_name = "PE 16:1-12:0 - -(H20)"
        goslin_parser.parse(lipid_name)
        assert goslin_parser.word_in_grammar
        goslin_parser.raise_events()
        
        
        
        
        
        
    def test_lipid_fragment_success(self):
    
        goslin_fragment_parser_event_handler = GoslinFragmentParserEventHandler()
        goslin_fragment_parser = Parser(goslin_fragment_parser_event_handler, "data/goslin/GoslinFragments.g4", ParserTest.PARSER_QUOTE)
        
        
        lipid_name = "PE 16:1-12:0 - -(H20)"
        goslin_fragment_parser.parse(lipid_name)
        assert goslin_fragment_parser.word_in_grammar
        goslin_fragment_parser.raise_events()
        
        assert goslin_fragment_parser_event_handler.lipid.fragment != None
        assert goslin_fragment_parser_event_handler.lipid.fragment.name == "-(H20)"
        
        
        
        
        
        
    def test_lipid_names(self):
        goslin_parser_event_handler = GoslinParserEventHandler()
        goslin_parser = Parser(goslin_parser_event_handler, "data/goslin/Goslin.g4", ParserTest.PARSER_QUOTE)
        
        ## glycerophospholipid
        lipid_name = "PE 16:1/12:0"
        goslin_parser.parse(lipid_name)
        assert goslin_parser.word_in_grammar
        goslin_parser.raise_events()
        
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "PE 16:1/12:0"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "PE 16:1_12:0"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "PE 28:1"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "PE"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "GP"
        
        ## sphingolipid
        lipid_name = "Cer 16:1;2/12:0"
        goslin_parser.parse(lipid_name)
        assert goslin_parser.word_in_grammar
        goslin_parser.raise_events()
        
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "Cer 16:1;2/12:0"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "Cer 16:1;2_12:0"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "Cer 28:1;2"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "Cer"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "SP"
        
        ## glycerolipid
        lipid_name = "TAG 16:1/12:0/20:2"
        goslin_parser.parse(lipid_name)
        assert goslin_parser.word_in_grammar
        goslin_parser.raise_events()
        
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "TAG 16:1/12:0/20:2"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "TAG 16:1_12:0_20:2"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "TAG 48:3"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "TG"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "GL"
        
        ## sterol
        lipid_name = "ChE 16:1"
        goslin_parser.parse(lipid_name)
        assert goslin_parser.word_in_grammar
        goslin_parser.raise_events()
        
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "ChE 16:1"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "ChE 16:1"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "ChE 16:1"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "ChE"
        assert goslin_parser_event_handler.lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "ST"
        
        
        
        
        
    
    def test_parser_read(self):
        lipidnames = []
        with open("pygoslin/tests/lipidnames.txt", mode = "rt") as infile:
            for line in infile:
                line = line.strip().strip(" ")
                if len(line) < 2: continue
                if line[0] == "#": continue
                lipidnames.append(Parser.split_string(line, ",", "\"")[0].strip("\""))
        
        goslin_parser_event_handler = GoslinParserEventHandler()
        goslin_parser = Parser(goslin_parser_event_handler, "data/goslin/Goslin.g4", ParserTest.PARSER_QUOTE)
        
        
        goslin_fragment_parser_event_handler = GoslinFragmentParserEventHandler()
        fragment_parser = Parser(goslin_fragment_parser_event_handler, "data/goslin/GoslinFragments.g4", ParserTest.PARSER_QUOTE)
        
        
        lipid_maps_parser_event_handler = LipidMapsParserEventHandler()
        lipid_maps_parser = Parser(lipid_maps_parser_event_handler, "data/goslin/LipidMaps.g4", ParserTest.PARSER_QUOTE)
        
        parsers = [[goslin_parser, goslin_parser_event_handler],
                   [fragment_parser, goslin_fragment_parser_event_handler],
                   [lipid_maps_parser, lipid_maps_parser_event_handler]]
        
        
        for lipidname in lipidnames:
            found = False
            for pp in parsers:
                prs, handler = pp
                prs.parse(lipidname)
                if prs.word_in_grammar:
                    found = True
                    break
                    prs.raise_events()
            assert found