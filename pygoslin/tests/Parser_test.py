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


from pygoslin.parser.Parser import *
from pygoslin.parser.GoslinParserEventHandler import GoslinParserEventHandler
from pygoslin.parser.LipidMapsParserEventHandler import LipidMapsParserEventHandler
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.LipidMolecularSpecies import LipidMolecularSpecies
from pygoslin.domain.FattyAcid import FattyAcid
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from random import randint

lipid_parser = LipidParser()
swiss_lipids_parser = SwissLipidsParser()
goslin_parser = GoslinParser()
lipid_maps_parser = LipidMapsParser()
hmdb_parser = HmdbParser()

class ParserTest(unittest.TestCase):
    PARSER_QUOTE = '\''
    
    def test_lipid_parser(self):
        global lipid_parser
        
        lipid_name = "PE 16:1-12:0[M+H]1+"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE 16:1_12:0[M+H]1+"
        
        lipid_name = "PE O-16:1p/12:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE P-16:0/12:0"
        
        lipid_name = "PAT16 16:1/12:0/14:1/8:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PAT16 16:1/12:0/14:1/8:0"
        
        lipid_name = "SLBPA 16:1/12:0/14:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        
        assert lipid.get_lipid_string() == "SLBPA 16:1_12:0_14:1", "w: %s" % lipid.get_lipid_string()
        
        lipid_name = "MLCL 16:1/12:0/14:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "LCL 16:1_12:0_14:1"
        
        
        lipid_name = "DLCL 14:1/8:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "DLCL 14:1_8:0"
        
        
        lipid_name = "PE O-12:1p/10:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE P-12:0/10:0"
        
        
        lipid_name = "PE 12:1/10:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE 12:1/10:0"
        
        
        lipid_name = "PIP[3'] 17:0/20:4(5Z,8Z,11Z,14Z)"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PIP(3') 17:0/20:4(5Z,8Z,11Z,14Z)"
        
        lipid_name = "AC2SGL 12:0-14:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "AC2SGL 12:0_14:1"
        
        lipid_name = "NAPE 16:1(6Z)/12:0/14:1(5Z)"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE-N(FA 14:1(5Z)) 16:1(6Z)/12:0"
        
        lipid_name = "PE-NMe 12:1(6Z)/10:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE-NMe 12:1(6Z)/10:0"
        
        lipid_name = "PIMIP 12:0-14:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PIMIP 12:0_14:1"
        
        lipid_name = "LCDPDAG 24:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "LCDPDAG 24:1"
        
        lipid_name = "LPIMIP 10:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "LPIMIP 10:0"
        
        
        
    def test_extended_class(self):
        global lipid_parser
        
        lipid_name = "PE O-16:1a-12:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_extended_class() == "PE-O"
        
        lipid_name = "LPE O-16:2p"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_extended_class() == "LPE-P", "wrong: %s" % lipid.get_extended_class()
        
        lipid_name = "PC O-16:1p/12:0"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_extended_class() == "PC-P"
        
        lipid_name = "PC O-16:1p"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_extended_class() == "PC-P"
        
        lipid_name = "LPC O-16:1a"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_extended_class() == "LPC-O"
        
        lipid_name = "LPC O-16:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_extended_class() == "LPC-O"
        
        lipid_name = "LPC 16:1"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_extended_class() == "LPC"
        
        lipid_name = "PC 16:1/12:2"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_extended_class() == "PC"
        
        
        
    def test_oh_count(self):
        lipid_parser = LipidParser()
        lipid = lipid_parser.parse("Cer 36:1;2");
        ohCount = lipid.lipid.info.get_total_functional_group_count("OH")
        assert 2 == ohCount
        asdCount = lipid.lipid.info.get_total_functional_group_count("ASD")
        assert 0 == asdCount
        lipid = lipid_parser.parse("Cer d36:1");
        ohCount = lipid.lipid.info.get_total_functional_group_count("OH")
        assert 2 == ohCount
        lipid = lipid_parser.parse("Cer 18:1;2/18:0");
        ohCount = lipid.lipid.info.get_total_functional_group_count("OH")
        assert 2 == ohCount
        lipid = lipid_parser.parse("Cer d18:1/18:0");
        ohCount = lipid.lipid.info.get_total_functional_group_count("OH")
        assert 2 == ohCount
        lipid = lipid_parser.parse("Cer 18:1;(OH)2/18:0");
        ohCount = lipid.lipid.info.get_total_functional_group_count("OH")
        assert 2 == ohCount
        
        
    def test_lipid_maps(self):
        lipid_maps_parser = LipidMapsParser()
        
        for lipid_name_input, lipid_name_output in [["PA(16:1/12:0)", "PA 16:1/12:0"],
                           ["PC(O-14:0/0:0)", "LPC O-14:0/0:0"],
                           ["SQMG(16:1(11Z)/0:0)", "SQMG 16:1(11Z)/0:0"],
                           ["TG(13:0/22:3(10Z,13Z,16Z)/22:5(7Z,10Z,13Z,16Z,19Z))[iso6]", "TG 13:0/22:3(10Z,13Z,16Z)/22:5(7Z,10Z,13Z,16Z,19Z)"],
                           ["13R-HODE", "13R-HODE"],
                           ["CL(1'-[20:0/20:0],3'-[20:4(5Z,8Z,11Z,14Z)/18:2(9Z,12Z)])", "CL 20:0/20:0/20:4(5Z,8Z,11Z,14Z)/18:2(9Z,12Z)"],
                           ["PA(P-20:0/18:3(6Z,9Z,12Z))", "PA P-20:0/18:3(6Z,9Z,12Z)"],
                           ["M(IP)2C(t18:0/20:0(2OH))", "M(IP)2C(1) 18:0;3OH,4OH/20:0;2OH"],
                           ["Cer(d16:2(4E,6E)/22:0(2OH))", "Cer 16:2(4E,6E);1OH,3OH/22:0;2OH"],
                           ["MG(18:1(11E)/0:0/0:0)[rac]", "MG 18:1(11E)/0:0/0:0"],
                           ["PAT18(24:1(2E)(2Me,4Me[S],6Me[S])/25:1(2E)(2Me,4Me[S],6Me[S])/26:1(2E)(2Me,4Me[S],6Me[S])/24:1(2E)(2Me,4Me[S],6Me[S]))", "PAT18 24:1(2E);2Me,4Me,6Me/25:1(2E);2Me,4Me,6Me/26:1(2E);2Me,4Me,6Me/24:1(2E);2Me,4Me,6Me"],
                           ["(3'-sulfo)Galbeta-Cer(d18:1/20:0)", "SHexCer 18:1;O2/20:0"],
                           ["GlcCer(d15:2(4E,6E)/22:0(2OH))", "GlcCer(1) 15:2(4E,6E);3OH/22:0;2OH"]
                          ]:
            lipid = lipid_maps_parser.parse(lipid_name_input)
            assert lipid_maps_parser.word_in_grammar
            assert lipid != None
            assert lipid.get_lipid_string() == lipid_name_output, "wrong: %s" % lipid.get_lipid_string()
        
        
    @unittest.expectedFailure
    def test_LP(self):
        global lipid_parser
        lipid = lipid_parser.parse("LP 19:1p")
        
        
    def test_hydroxyls(self):
        global goslin_parser, swiss_lipids_parser, lipid_maps_parser, hmdb_parser
        
        lipid = swiss_lipids_parser.parse("Cer(d18:1(4E)/24:0-2OH)")
        assert lipid != None
        assert lipid.get_lipid_string() == "Cer 18:1(4);(OH)2/24:0;OH"
        assert lipid.get_sum_formula() == "C42H83NO4"
        assert abs(lipid.get_mass() - 665.632209) < 1e-3
        
        
        lipid = swiss_lipids_parser.parse("Cer(d18:1(4E)/24:0(2OH))")
        assert lipid != None
        assert lipid.get_lipid_string() == "Cer 18:1(4);(OH)2/24:0;OH"
        assert lipid.get_sum_formula() == "C42H83NO4"
        assert abs(lipid.get_mass() - 665.632209) < 1e-3
        
        
        lipid = lipid_maps_parser.parse("Cer(d18:1(4E)/24:0(2OH))")
        assert lipid != None
        assert lipid.get_lipid_string() == "Cer 18:1(4E);1OH,3OH/24:0;2OH"
        assert lipid.get_sum_formula() == "C42H83NO4"
        assert abs(lipid.get_mass() - 665.632209) < 1e-3
        
        
        lipid = goslin_parser.parse("Cer 18:1(4E);2/24:0;1")
        assert lipid != None
        assert lipid.get_lipid_string() == "Cer 18:1(4);(OH)2/24:0;OH"
        assert lipid.get_sum_formula() == "C42H83NO4"
        assert abs(lipid.get_mass() - 665.632209) < 1e-3
        
        lipid = hmdb_parser.parse("SM(d18:1(4)/16:1(9Z)(OH))")
        assert lipid != None
        assert lipid.get_lipid_string() == "SM 18:1(4);OH/16:1(9);OH"
        assert lipid.get_sum_formula() == "C39H77N2O7P"
        assert abs(lipid.get_mass() - 716.546841) < 1e-3


    def test_lyso(self):
        global lipid_parser
        
        lipid_name = "LPA O-16:1a"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "LPA O-16:1", "w: %s" % lipid.get_lipid_string()
        
        lipid_name = "LPC O-16:1a"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "LPC O-16:1"
        
        lipid_name = "LPE O-16:1p"
        lipid = lipid_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "LPE P-16:0", "w: %s" % lipid.get_lipid_string()
            
            
            
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
        assert lipid.get_lipid_string(LipidLevel.FULL_STRUCTURE) == "PG 22:1(5Z)/12:0"
        assert lipid.get_lipid_string(LipidLevel.STRUCTURE_DEFINED) == "PG 22:1(5)/12:0"
        assert lipid.get_lipid_string(LipidLevel.SN_POSITION) == "PG 22:1/12:0"
        assert lipid.get_lipid_string(LipidLevel.MOLECULAR_SPECIES) == "PG 22:1_12:0"
        assert lipid.get_lipid_string(LipidLevel.SPECIES) == "PG 34:1"
        assert lipid.get_lipid_string(LipidLevel.CLASS) == "PG"
        assert lipid.get_lipid_string(LipidLevel.CATEGORY) == "GP"
        
        lipid_name = "Cer 28:1;2"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "Cer 28:1;O2"
        
        lipid_name = "DAG 38:1"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "DG 38:1"
        
        
        
        
    def test_info_level(self):
        global swiss_lipids_parser
        global lipid_maps_parser
        
        lipid_name = "PG(22:1(5Z)/12:0)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string(LipidLevel.FULL_STRUCTURE) == "PG 22:1(5Z)/12:0"
        assert lipid.get_lipid_string(LipidLevel.STRUCTURE_DEFINED) == "PG 22:1(5)/12:0"
        assert lipid.get_lipid_string(LipidLevel.SN_POSITION) == "PG 22:1/12:0"
        assert lipid.get_lipid_string(LipidLevel.MOLECULAR_SPECIES) == "PG 22:1_12:0"
        assert lipid.get_lipid_string(LipidLevel.SPECIES) == "PG 34:1"
        assert lipid.get_lipid_string(LipidLevel.CLASS) == "PG"
        assert lipid.get_lipid_string(LipidLevel.CATEGORY) == "GP"
        
        lipid_name = "PG(22:1(5Z)/12:0)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.lipid.info.level == LipidLevel.FULL_STRUCTURE
        assert lipid.get_lipid_string() == "PG 22:1(5Z)/12:0"
        
        lipid_name = "PG(22:1(5)/12:0)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid
        assert lipid.lipid.info.level == LipidLevel.STRUCTURE_DEFINED
        assert lipid.get_lipid_string() == "PG 22:1(5)/12:0"
        
        lipid_name = "PG(22:1/12:0)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid
        assert lipid.lipid.info.level == LipidLevel.SN_POSITION
        assert lipid.get_lipid_string() == "PG 22:1/12:0"
        
        lipid_name = "PG(22:1_12:0)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid
        assert lipid.lipid.info.level == LipidLevel.MOLECULAR_SPECIES
        assert lipid.get_lipid_string() == "PG 22:1_12:0"
        
        lipid_name = "LPG(O-22:1)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid
        assert lipid.lipid.info.level == LipidLevel.SPECIES
        assert lipid.get_lipid_string() == "LPG O-22:1"
        
        
        lipid_name = "PG(22:1(5Z)/12:0)"
        lipid = lipid_maps_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string(LipidLevel.FULL_STRUCTURE) == "PG 22:1(5Z)/12:0"
        assert lipid.get_lipid_string(LipidLevel.STRUCTURE_DEFINED) == "PG 22:1(5)/12:0"
        assert lipid.get_lipid_string(LipidLevel.SN_POSITION) == "PG 22:1/12:0"
        assert lipid.get_lipid_string(LipidLevel.MOLECULAR_SPECIES) == "PG 22:1_12:0"
        assert lipid.get_lipid_string(LipidLevel.SPECIES) == "PG 34:1"
        assert lipid.get_lipid_string(LipidLevel.CLASS) == "PG"
        assert lipid.get_lipid_string(LipidLevel.CATEGORY) == "GP"
        
        
        
        
        
    def test_lipid_names(self):
        file_name = os.path.join("pygoslin", "data", "goslin", "Goslin.g4")
        goslin_parser = Parser(GoslinParserEventHandler(), file_name, ParserTest.PARSER_QUOTE)
        
        ## glycerophospholipid
        lipid_name = "PE 16:1(3E)/12:0"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURE_DEFINED) == "PE 16:1(3)/12:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.SN_POSITION) == "PE 16:1/12:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SPECIES) == "PE 16:1_12:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "PE 28:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "PE"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "GP"
        
        ## sphingolipid
        lipid_name = "Cer 16:1(5Z);2/12:0"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURE_DEFINED) == "Cer 16:1(5);(OH)2/12:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.SN_POSITION) == "Cer 16:1;O2/12:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SPECIES) == "Cer 16:1;O2/12:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "Cer 28:1;O2"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "Cer"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "SP"
        
        ## glycerolipid
        lipid_name = "TAG 16:1(5Z)/12:0/20:2(8Z,11Z)"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURE_DEFINED) == "TG 16:1(5)/12:0/20:2(8,11)"
        assert lipid.get_lipid_string(LipidLevel.SN_POSITION) == "TG 16:1/12:0/20:2"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SPECIES) == "TG 16:1_12:0_20:2"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "TG 48:3"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "TG"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "GL"
        
        
        
        
        ## sterol
        lipid_name = "ChE 16:1(4Z)"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURE_DEFINED) == "SE 27:1/16:1(4)"
        assert lipid.lipid.get_lipid_string(LipidLevel.SN_POSITION) == "SE 27:1/16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SPECIES) == "SE 27:1/16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "SE 27:1/16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "SE 27:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "ST"
        
        ## sterol
        lipid_name = "SE 27:1/16:1(4Z)"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURE_DEFINED) == "SE 27:1/16:1(4)"
        assert lipid.lipid.get_lipid_string(LipidLevel.SN_POSITION) == "SE 27:1/16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SPECIES) == "SE 27:1/16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "SE 27:1/16:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "SE 27:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "ST"
        
        
        
        
        ## PC
        lipid_name = "PC O-18:1(4E)a/16:0"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        
        assert lipid.lipid.get_lipid_string(LipidLevel.STRUCTURE_DEFINED) == "PC O-18:1(4)/16:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.SN_POSITION) == "PC O-18:1/16:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.MOLECULAR_SPECIES) == "PC O-18:1_16:0"
        assert lipid.lipid.get_lipid_string(LipidLevel.SPECIES) == "PC O-34:1"
        assert lipid.lipid.get_lipid_string(LipidLevel.CLASS) == "PC"
        assert lipid.lipid.get_lipid_string(LipidLevel.CATEGORY) == "GP"
        
        
    def test_adduct(self):
        goslin_parser_event_handler = GoslinParserEventHandler()
        file_name = os.path.join("pygoslin", "data", "goslin", "Goslin.g4")
        goslin_parser = Parser(goslin_parser_event_handler, file_name)

        lipid_name = "PE 16:1/12:0[M+H]1+"
        lipid = goslin_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE 16:1/12:0[M+H]1+"
        
        
        
        
      
    def test_swiss_lipids(self):
        global swiss_lipids_parser
        
        lipid_name = "TG(O-16:0/18:3(6Z,9Z,12Z)/18:1(11Z))"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "TG O-16:0/18:3(6Z,9Z,12Z)/18:1(11Z)"
        
        lipid_name = "PIP2[4,5](21:0/24:4(9Z,12Z,15Z,18Z))"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PIP2(4',5') 21:0/24:4(9Z,12Z,15Z,18Z)"
        
        lipid_name = "GalGb3Cer(d18:0/14:0)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "GalGb3Cer 18:0;OH/14:0", "wrong %s" % lipid.get_lipid_string()
        
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
        assert lipid.get_lipid_string() == "PE-N(FA 12:0) 12:0/30:4(15Z,18Z,21Z,24Z)"
        
        lipid_name = "NAPE (12:0/0:0/12:0)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "LPE-N(FA 12:0) 12:0/0:0"
        
        lipid_name = "NAPE (30:5(12Z,15Z,18Z,21Z,24Z)/18:0/8:0)"
        lipid = swiss_lipids_parser.parse(lipid_name)
        assert lipid != None
        assert lipid.get_lipid_string() == "PE-N(FA 8:0) 30:5(12Z,15Z,18Z,21Z,24Z)/18:0"
        assert lipid.get_sum_formula() == "C61H110NO9P"
      
      
      
      
    
    def test_parser_read(self):
        global lipid_parser
        lipidnames = []
        file_name = os.path.join("pygoslin", "data", "goslin", "testfiles", "lipidnames.csv")
        with open(file_name, mode = "rt", encoding= "utf-8") as infile:
            for line in infile:
                line = line.strip().strip(" ")
                if len(line) < 2: continue
                if line[0] == "#": continue
                lipidnames.append(Parser.split_string(line, ",", "\"")[0].strip("\""))
        
        for lipid_name in lipidnames:
            try:
                lipid = lipid_parser.parse(lipid_name)
            except Exception as e:
                print("Lipid could not be parsed: '%s'" % lipid_name)
                raise e
            assert lipid != None
