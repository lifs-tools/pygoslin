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


from pygoslin.parser.Parser import *
from pygoslin.parser.GoslinParserEventHandler import GoslinParserEventHandler
from pygoslin.parser.GoslinFragmentParserEventHandler import GoslinFragmentParserEventHandler
from pygoslin.parser.LipidMapsParserEventHandler import LipidMapsParserEventHandler
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.LipidMolecularSubspecies import LipidMolecularSubspecies
from pygoslin.domain.FattyAcid import FattyAcid
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from random import randint

lipid_parser = LipidParser()
swiss_lipids_parser = SwissLipidsParser()
goslin_parser = GoslinParser()
goslin_fragment_parser = GoslinFragmentParser()
lipid_maps_parser = LipidMapsParser()
hmdb_parser = HmdbParser()

class ParserTest(unittest.TestCase):
    PARSER_QUOTE = '\''
    
    def test_lipid_parser(self):
        global lipid_parser
        
        lipid_name = "PE 16:1-12:0[M+H]1+"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE 16:1-12:0[M+H]1+"
        
        lipid_name = "PA 16:1-12:0 - fragment"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PA 16:1-12:0"
        assert lipid.get_lipid_fragment_string() == "PA 16:1-12:0 - fragment"
        
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
        assert lipid.get_lipid_string() == "PIP[3'] 17:0/20:4(5Z,8Z,11Z,14Z)"
        
        lipid_name = "AC2SGL 12:0-14:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "AC2SGL 12:0-14:1"
        
        lipid_name = "NAPE 16:1(6Z)/12:0/14:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "NAPE 16:1(6Z)/12:0/14:1"
        
        lipid_name = "PE-NMe 12:1(6Z)/10:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE-NMe 12:1(6Z)/10:0"
        
        lipid_name = "PIMIP 12:0-14:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PIMIP 12:0-14:1"
        
        lipid_name = "LCDPDAG 24:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "LCDPDAG 24:1"
        
        lipid_name = "LPIMIP 10:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "LPIMIP 10:0"
        
        
        
    def test_lipid_maps(self):
        lipid_maps_parser = LipidMapsParser()
        
        for lipid_name_input, lipid_name_output in [["PA(16:1/12:0)", "PA 16:1/12:0"],
                           ["PC(O-14:0/0:0)", "LPC O-14:0a"],
                           ["SQMG(16:1(11Z)/0:0)", "SQMG 16:1(11Z)"],
                           ["TG(13:0/22:3(10Z,13Z,16Z)/22:5(7Z,10Z,13Z,16Z,19Z))[iso6]", "TAG 13:0/22:3(10Z,13Z,16Z)/22:5(7Z,10Z,13Z,16Z,19Z)"],
                           ["13R-HODE", "13R-HODE"],
                           ["CL(1'-[20:0/20:0],3'-[20:4(5Z,8Z,11Z,14Z)/18:2(9Z,12Z)])", "CL 20:0/20:0/20:4(5Z,8Z,11Z,14Z)/18:2(9Z,12Z)"],
                           ["PA(P-20:0/18:3(6Z,9Z,12Z))", "PA 20:1p/18:3(6Z,9Z,12Z)"],
                           ["M(IP)2C(t18:0/20:0(2OH))", "M(IP)2C 18:0;3/20:0;1"],
                           ["Cer(d16:2(4E,6E)/22:0(2OH))", "Cer 16:2(4E,6E);2/22:0;1"],
                           ["MG(18:1(11E)/0:0/0:0)[rac]", "MAG 18:1(11E)"],
                           ["PAT18(24:1(2E)(2Me,4Me[S],6Me[S])/25:1(2E)(2Me,4Me[S],6Me[S])/26:1(2E)(2Me,4Me[S],6Me[S])/24:1(2E)(2Me,4Me[S],6Me[S]))", "PAT18 24:1(2E)/25:1(2E)/26:1(2E)/24:1(2E)"],
                           ["(3'-sulfo)Galbeta-Cer(d18:1/20:0)", "SHexCer 18:1;2/20:0"],
                           ["GlcCer(d15:2(4E,6E)/22:0(2OH))", "HexCer 15:2(4E,6E);2/22:0;1"]
                          ]:
            lipid = lipid_maps_parser.parse(lipid_name_input)
            assert lipid_maps_parser.word_in_grammar
            assert lipid != None
            assert lipid.get_lipid_string() == lipid_name_output
        
        
    @unittest.expectedFailure
    def test_LP(self):
        global lipid_parser
        lipid = lipid_parser.parse("LP 19:1p")
        
        
    def test_hydroxyls(self):
        global goslin_parser, swiss_lipids_parser, lipid_maps_parser, hmdb_parser
        
        lipid = swiss_lipids_parser.parse("Cer(d18:1(4E)/24:0-2OH)")
        assert lipid != None
        assert lipid.get_lipid_string() == "Cer 18:1(4E);2/24:0;1"
        assert lipid.get_sum_formula() == "C42H83NO4"
        assert abs(lipid.get_mass() - 665.632209) < 1e-3
        
        
        lipid = swiss_lipids_parser.parse("Cer(d18:1(4E)/24:0(2OH))")
        assert lipid != None
        assert lipid.get_lipid_string() == "Cer 18:1(4E);2/24:0;1"
        assert lipid.get_sum_formula() == "C42H83NO4"
        assert abs(lipid.get_mass() - 665.632209) < 1e-3
        
        
        lipid = lipid_maps_parser.parse("Cer(d18:1(4E)/24:0(2OH))")
        assert lipid != None
        assert lipid.get_lipid_string() == "Cer 18:1(4E);2/24:0;1"
        assert lipid.get_sum_formula() == "C42H83NO4"
        assert abs(lipid.get_mass() - 665.632209) < 1e-3
        
        
        lipid = goslin_parser.parse("Cer 18:1(4E);2/24:0;1")
        assert lipid != None
        assert lipid.get_lipid_string() == "Cer 18:1(4E);2/24:0;1"
        assert lipid.get_sum_formula() == "C42H83NO4"
        assert abs(lipid.get_mass() - 665.632209) < 1e-3
        
        lipid = hmdb_parser.parse("SM(d18:1/16:1(9Z)(OH))")
        assert lipid != None
        assert lipid.get_lipid_string() == "SM 18:1;2/16:1(9Z);1"
        assert lipid.get_sum_formula() == "C39H77N2O7P"
        assert abs(lipid.get_mass() - 716.546841) < 1e-3


    def test_lyso(self):
        global lipid_parser
        
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
            
            
            
    @unittest.expectedFailure
    def test_lpe_fail(self):
        global lipid_parser
        lipid_name = "LPE O-16:1p/12:0"
        lipid = lipid_parser.parse(lipid_name)
        
        
    @unittest.expectedFailure
    def test_lipid_parser_fail(self):
        global lipid_parser
        lipid_name = "fail"
        lipid = lipid_parser.parse(lipid_name)
        
        
        
    def test_species_level(self):
        global goslin_parser
        
        lipid_name = "PG 22:1(5Z)/12:0"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string(LipidLevel.ISOMERIC_SUBSPECIES) == "PG 22:1(5Z)/12:0"
        assert lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "PG 22:1(5Z)/12:0"
        assert lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "PG 22:1(5Z)-12:0"
        assert lipid.get_lipid_string(LipidLevel.SPECIES) == "PG 34:1"
        assert lipid.get_lipid_string(LipidLevel.CLASS) == "PG"
        assert lipid.get_lipid_string(LipidLevel.CATEGORY) == "GP"
        
        lipid_name = "Cer 28:1;2"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "Cer 28:1;2"
        
        lipid_name = "DAG 38:1"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "DAG 38:1"
        
        
        
        
    def test_info_level(self):
        global swiss_lipids_parser
        global lipid_maps_parser
        
        lipid_name = "PG(22:1(5Z)/12:0)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string(LipidLevel.ISOMERIC_SUBSPECIES) == "PG 22:1(5Z)/12:0"
        assert lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "PG 22:1(5Z)/12:0"
        assert lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "PG 22:1(5Z)-12:0"
        assert lipid.get_lipid_string(LipidLevel.SPECIES) == "PG 34:1"
        assert lipid.get_lipid_string(LipidLevel.CLASS) == "PG"
        assert lipid.get_lipid_string(LipidLevel.CATEGORY) == "GP"
        
        lipid_name = "PG(22:1(5Z)/12:0)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.lipid.info.level == LipidLevel.ISOMERIC_SUBSPECIES
        assert lipid.get_lipid_string() == "PG 22:1(5Z)/12:0"
        
        lipid_name = "PG(22:1/12:0)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid
        assert lipid.lipid.info.level == LipidLevel.STRUCTURAL_SUBSPECIES
        assert lipid.get_lipid_string() == "PG 22:1/12:0"
        
        lipid_name = "PG(22:1_12:0)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid
        assert lipid.lipid.info.level == LipidLevel.MOLECULAR_SUBSPECIES
        assert lipid.get_lipid_string() == "PG 22:1-12:0"
        
        lipid_name = "LPG(O-22:1)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid
        assert lipid.lipid.info.level == LipidLevel.SPECIES
        assert lipid.get_lipid_string() == "LPG 22:1a"
        
        
        lipid_name = "PG(22:1(5Z)/12:0)"
        lipid = lipid_maps_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string(LipidLevel.ISOMERIC_SUBSPECIES) == "PG 22:1(5Z)/12:0"
        assert lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "PG 22:1(5Z)/12:0"
        assert lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "PG 22:1(5Z)-12:0"
        assert lipid.get_lipid_string(LipidLevel.SPECIES) == "PG 34:1"
        assert lipid.get_lipid_string(LipidLevel.CLASS) == "PG"
        assert lipid.get_lipid_string(LipidLevel.CATEGORY) == "GP"
        
        
        
    def test_mediators(self):
        global goslin_parser
        
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
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "PE 16:1-12:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "PE 28:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "PE"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "GP"
        
        ## sphingolipid
        lipid_name = "Cer 16:1;2/12:0"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "Cer 16:1;2/12:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "Cer 16:1;2-12:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "Cer 28:1;2"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "Cer"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "SP"
        
        ## glycerolipid
        lipid_name = "TAG 16:1/12:0/20:2"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "TAG 16:1/12:0/20:2"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "TAG 16:1-12:0-20:2"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "TAG 48:3"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "TAG"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "GL"
        
        
        
        
        ## sterol
        lipid_name = "ChE 16:1"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "SE 27:1/16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "SE 27:1/16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "SE 27:1/16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "SE 27:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "ST"
        
        ## sterol
        lipid_name = "SE 27:1/16:1"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "SE 27:1/16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "SE 27:1/16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "SE 27:1/16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "SE 27:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "ST"
        
        
        
        
        ## PC
        lipid_name = "PC O-18:1a/16:0"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES) == "PC O-18:1a/16:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES) == "PC O-18:1a-16:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "PC O-34:1a"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "PC"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "GP"
        
        
    def test_adduct(self):
        goslin_parser_event_handler = GoslinParserEventHandler()
        goslin_parser = Parser(goslin_parser_event_handler, "pygoslin/data/goslin/Goslin.g4")

        lipid_name = "PE 16:1/12:0[M+H]1+"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE 16:1/12:0[M+H]1+"
        
        
        
        
      
    def test_swiss_lipids(self):
        global swiss_lipids_parser
        
        lipid_name = "TG(O-16:0/18:3(6Z,9Z,12Z)/18:1(11Z))"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "TAG 16:0a/18:3(6Z,9Z,12Z)/18:1(11Z)"
        
        lipid_name = "PIP2[4,5](21:0/24:4(9Z,12Z,15Z,18Z))"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PIP2[4',5'] 21:0/24:4(9Z,12Z,15Z,18Z)"
        
        lipid_name = "GalGb3Cer(d18:0/14:0)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "GalGb3Cer 18:0;2/14:0"
        
        lipid_name = "PI(34:5(19Z,22Z,25Z,28Z,31Z)/18:1(6Z))"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PI 34:5(19Z,22Z,25Z,28Z,31Z)/18:1(6Z)"
        
        lipid_name = "PIP(22:6/20:4)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PIP 22:6/20:4"
        
        lipid_name = "NAPE (12:0/30:4(15Z,18Z,21Z,24Z)/12:0)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "NAPE 12:0/30:4(15Z,18Z,21Z,24Z)/12:0"
      
      
      
      
    
    def test_parser_read(self):
        global lipid_parser
        lipidnames = []
        with open("pygoslin/tests/lipidnames.txt", mode = "rt") as infile:
            for line in infile:
                line = line.strip().strip(" ")
                if len(line) < 2: continue
                if line[0] == "#": continue
                lipidnames.append(Parser.split_string(line, ",", "\"")[0].strip("\""))
        
        for lipid_name in lipidnames:
            lipid = lipid_parser.parse(lipid_name)
            assert lipid != None
