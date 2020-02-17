import unittest
import os


import pyximport
pyximport.install(setup_args = {"script_args" : ["--force"]}, language_level=3)

import pygoslin
from pygoslin.parser.Parser import LipidParser
from pygoslin.parser.Parser import Parser, GoslinParser, LipidParser, LipidMapsParser
from pygoslin.parser.GoslinParserEventHandler import GoslinParserEventHandler
from pygoslin.parser.GoslinFragmentParserEventHandler import GoslinFragmentParserEventHandler
from pygoslin.parser.LipidMapsParserEventHandler import LipidMapsParserEventHandler
from pygoslin.domain.LipidLevel import LipidLevel
from random import randint

class ParserTest(unittest.TestCase):
    PARSER_QUOTE = '\''
    
        
    
    
    def test_lipid_parser(self):
        lipid_parser = LipidParser()
        
        lipid_name = "PE 16:1-12:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE 16:1_12:0"
        
        lipid_name = "PA 16:1-12:0 - fragment"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PA 16:1_12:0"
        assert lipid.get_lipid_fragment_string() == "PA 16:1_12:0 - fragment"
        
        lipid_name = "PE O-16:1p/12:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE O-16:1p/12:0"
        
        lipid_name = "PAT16 16:1/12:0/14:1/8:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PAT16 16:1/12:0/14:1/8:0"
        
        lipid_name = "SLBPA 16:1/12:0/14:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "SLBPA 16:1/12:0/14:1"
        
        lipid_name = "MLCL 16:1/12:0/14:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "MLCL 16:1/12:0/14:1"
        
        
        lipid_name = "DLCL 14:1/8:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "DLCL 14:1/8:0"
        
        
        lipid_name = "PE O-12:1p/10:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE O-12:1p/10:0"
        
        
        lipid_name = "PE 12:1/10:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE 12:1/10:0"
        
        
        lipid_name = "PIP[3'] 17:0/20:4(5Z,8Z,11Z,14Z)"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PIP[3'] 17:0/20:4"
        
        lipid_name = "AC2SGL 12:0-14:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "AC2SGL 12:0_14:1"
        
        lipid_name = "NAPE 16:1/12:0/14:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "NAPE 16:1/12:0/14:1"
        
        lipid_name = "PE-NMe 12:1/10:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE-NMe 12:1/10:0"
        
        
        
    def test_lipid_maps(self):
        lipid_maps_parser = LipidMapsParser()
        
        for lipid_name_input, lipid_name_output in [["PA(16:1/12:0)", "PA 16:1/12:0"],
                           ["PC(O-14:0/0:0)", "LPC O-14:0a"],
                           ["SQMG(16:1(11Z)/0:0)", "SQMG 16:1"],
                           ["TG(13:0/22:3(10Z,13Z,16Z)/22:5(7Z,10Z,13Z,16Z,19Z))[iso6]", "TAG 13:0/22:3/22:5"],
                           ["13R-HODE", "13R-HODE"],
                           ["CL(1'-[20:0/20:0],3'-[20:4(5Z,8Z,11Z,14Z)/18:2(9Z,12Z)])", "CL 20:0/20:0/20:4/18:2"],
                           ["PA(P-20:0/18:3(6Z,9Z,12Z))", "PA 20:0p/18:3"],
                           ["M(IP)2C(t18:0/20:0(2OH))", "M(IP)2C 18:0;3/20:0;1"],
                           ["Cer(d16:2(4E,6E)/22:0(2OH))", "Cer 16:2;2/22:0;1"],
                           ["MG(18:1(11E)/0:0/0:0)[rac]", "MAG 18:1"],
                           ["PAT18(24:1(2E)(2Me,4Me[S],6Me[S])/25:1(2E)(2Me,4Me[S],6Me[S])/26:1(2E)(2Me,4Me[S],6Me[S])/24:1(2E)(2Me,4Me[S],6Me[S]))", "PAT18 24:1/25:1/26:1/24:1"],
                           ["(3'-sulfo)Galbeta-Cer(d18:1/20:0)", "SHexCer 18:1;2/20:0"],
                           ["GlcCer(d15:2(4E,6E)/22:0(2OH))", "HexCer 15:2;2/22:0;1"]
                          ]:
            lipid = lipid_maps_parser.parse(lipid_name_input)

            assert lipid != None
            assert lipid.get_lipid_string() == lipid_name_output
        


    def test_lyso(self):
        lipid_parser = LipidParser()
        
        lipid_name = "LPA 16:1a"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "LPA 16:1a"
        
        lipid_name = "LPC O-16:1a"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "LPC O-16:1a"
        
        lipid_name = "LPE O-16:1p"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "LPE O-16:1p"
        
        lipid_name = "LPE O-16:1p/12:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid == None
        
        
        
        
        
        
        
    def test_species_level(self):
        goslin_parser = GoslinParser()
        
        lipid_name = "Cer 28:1;2"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "Cer 28:1;2"
        
        lipid_name = "DAG 38:1"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "DAG 38:1"
        
        
        
        
    def test_mediators(self):
        goslin_parser = GoslinParser()
        
        for lipid_name in ["10-HDoHE","11-HDoHE","11-HETE","11,12-DHET","11(12)-EET", "12-HEPE","12-HETE","12-HHTrE","12-OxoETE","12(13)-EpOME","13-HODE","13-HOTrE","14,15-DHET","14(15)-EET","14(15)-EpETE","15-HEPE","15-HETE","15d-PGJ2","16-HDoHE","16-HETE","18-HEPE","5-HEPE","5-HETE","5-HpETE","5-OxoETE","5,12-DiHETE","5,6-DiHETE","5,6,15-LXA4","5(6)-EET","8-HDoHE","8-HETE","8,9-DHET","8(9)-EET","9-HEPE","9-HETE","9-HODE","9-HOTrE","9(10)-EpOME","AA","alpha-LA","DHA","EPA","Linoleic acid","LTB4","LTC4","LTD4","Maresin 1","Palmitic acid","PGB2","PGD2","PGE2","PGF2alpha","PGI2","Resolvin D1","Resolvin D2","Resolvin D3","Resolvin D5","tetranor-12-HETE","TXB1","TXB2","TXB3"]:
            
            lipid = goslin_parser.parse(lipid_name)
            assert lipid != None
            n = lipid.lipid.get_lipid_string()
            assert lipid.get_lipid_string() == lipid_name
    
    @unittest.expectedFailure
    def test_lipid_fragment_fail(self):
        goslin_parser = Parser(GoslinParserEventHandler(), "pygoslin/data/goslin/Goslin.g4", ParserTest.PARSER_QUOTE)
        
        
        lipid_name = "PE 16:1-12:0 - -(H20)"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        
        
        
        
        
        
    def test_lipid_fragment_success(self):
        goslin_fragment_parser = Parser(GoslinFragmentParserEventHandler(), "pygoslin/data/goslin/GoslinFragments.g4", ParserTest.PARSER_QUOTE)
        
        
        lipid_name = "PE 16:1-12:0 - -(H20)"
        lipid = goslin_fragment_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.fragment != None
        assert lipid.fragment.name == "-(H20)"
        
        
        
        
        
        
    def test_lipid_names(self):
        goslin_parser = Parser(GoslinParserEventHandler(), "pygoslin/data/goslin/Goslin.g4", ParserTest.PARSER_QUOTE)
        
        ## glycerophospholipid
        lipid_name = "PE 16:1/12:0"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "PE 16:1/12:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "PE 16:1_12:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "PE 28:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "PE"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "GP"
        
        ## sphingolipid
        lipid_name = "Cer 16:1;2/12:0"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "Cer 16:1;2/12:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "Cer 16:1;2_12:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "Cer 28:1;2"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "Cer"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "SP"
        
        ## glycerolipid
        lipid_name = "TAG 16:1/12:0/20:2"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "TAG 16:1/12:0/20:2"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "TAG 16:1_12:0_20:2"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "TAG 48:3"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "TAG"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "GL"
        
        ## sterol
        lipid_name = "ChE 16:1"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "ChE 16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "ChE 16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "ChE 16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "ChE"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "ST"
        
        
    def test_adduct(self):
        from pygoslin.parser.Parser import Parser
        from pygoslin.parser.GoslinParserEventHandler import GoslinParserEventHandler
        goslin_parser_event_handler = GoslinParserEventHandler()
        goslin_parser = Parser(goslin_parser_event_handler, "pygoslin/data/goslin/Goslin.g4")

        lipid_name = "PE 16:1/12:0[M+H]1+"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE 16:1/12:0[M+H]1+"
        
        
    
    def test_parser_read(self):
        lipidnames = []
        with open("pygoslin/tests/lipidnames.txt", mode = "rt") as infile:
            for line in infile:
                line = line.strip().strip(" ")
                if len(line) < 2: continue
                if line[0] == "#": continue
                lipidnames.append(Parser.split_string(line, ",", "\"")[0].strip("\""))
        
        
        lipid_parser = LipidParser()
        
        for lipid_name in lipidnames:
            lipid = lipid_parser.parse(lipid_name)
            assert lipid != None