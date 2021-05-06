from pygoslin.domain.FattyAcid import FattyAcid
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.Element import Element
from pygoslin.domain.FunctionalGroup import *
from pygoslin.domain.Cycle import Cycle
from pygoslin.domain.LipidAdduct import LipidAdduct
from pygoslin.domain.LipidIsomericSubspecies import LipidIsomericSubspecies
from pygoslin.domain.LipidMolecularSubspecies import LipidMolecularSubspecies
from pygoslin.domain.LipidSpecies import LipidSpecies
from pygoslin.domain.LipidLevel import LipidLevel

from pygoslin.parser.Parser import ShorthandParser


import unittest


class ShorthandTest(unittest.TestCase):
    
    def test_nomenclature(self):
        parser = ShorthandParser()
        
        data = {"PC 18:1(11Z)/16:0": ["PC 18:1(11Z)/16:0", "PC 18:1(11)/16:0", "PC 18:1_16:0", "PC 34:1", "C42H82NO8P"],
                
                "CL O-16:0;3Me,7Me,11Me,15Me/O-16:0;3Me,7Me,11Me,15Me/O-16:0;3Me,7Me,11Me,15Me/O-16:0;3Me,7Me,11Me,15Me": ["CL O-16:0;3Me,7Me,11Me,15Me/O-16:0;3Me,7Me,11Me,15Me/O-16:0;3Me,7Me,11Me,15Me/O-16:0;3Me,7Me,11Me,15Me", "CL O-16:0;(Me)4/O-16:0;(Me)4/O-16:0;(Me)4/O-16:0;(Me)4", "CL O-20:0_O-20:0_O-20:0_O-20:0", "CL eO-80:0"], # CL eO-80:0
                
                "PC O-16:0/O-18:1(9Z)": ["PC O-16:0/O-18:1(9Z)", "PC O-16:0/O-18:1(9)", "PC O-16:0_O-18:1", "PC dO-34:1"], # PC dO-34:1
                
                "Cer 18:0;1OH,3OH/16:0": ["Cer 18:0;1OH,3OH/16:0", "Cer 18:0;(OH)2/16:0", "Cer 18:0;O2/16:0", "Cer 34:0;O2"], # Cer 34:0;O2
                
                "IPC(1) 18:1(8E);3OH,4OH/24:0;2OH": ["IPC(1) 18:1(8E);3OH,4OH/24:0;2OH", "IPC 18:1(8);(OH)2/24:0;OH", "IPC 18:1;O2/24:0;O", "IPC 42:1;O4"], # IPC 42:1;O4
                
                "CerP(1) 18:1(4E);3OH/16:0;2OH": ["CerP(1) 18:1(4E);3OH/16:0;2OH", "CerP 18:1(4);OH/16:0;OH", "CerP 18:1;O/16:0;O", "CerP 34:1;O3", "C34H68NO7P"], # CerP 34:1;O3 / C34H68NO7P
                
                "Gal-Gal-Glc-Cer(1) 18:1(4E);3OH/24:0": ["Gal-Gal-Glc-Cer(1) 18:1(4E);3OH/24:0", "Gal-Gal-Glc-Cer 18:1(4);OH/24:0", "GalGalGlcCer 18:1;O/24:0", "GalGalGlcCer 42:1;O2", "C60H113NO18"], # GalGalGlcCer 42:1;O2 / C60H113NO18
                
                "PC 16:0/20:2(5Z,13E);[8-12cy5;11OH;9oxo];15OH": ["PC 16:0/20:2(5Z,13E);15OH;[8-12cy5:0;11OH;9oxo]", "PC 16:0/20:2(5,13);OH;[cy5:0;OH;oxo]", "PC 16:0_20:4;O3", "PC 36:4;O3"], # PC 36:4;O3
                
                "PE-N(FA 18:1(9Z)) 16:0/20:4(4Z,8Z,11Z,14Z)": ["PE-N(FA 18:1(9Z)) 16:0/20:4(4Z,8Z,11Z,14Z)", "PE-N(FA 18:1(9)) 16:0/20:4(4,8,11,14)", "PE-N(FA 18:1) 16:0_20:4", "PE-N(FA) 54:5"], # PE-N(FA) 54:5
                
                "LPC O-16:1(11Z)/0:0": ["LPC O-16:1(11Z)/0:0", "LPC O-16:1(11)/0:0", "LPC O-16:1", "LPC O-16:1"], # LPC O-16:1 / C24H50NO6P
                
                }
        """
                "PE P-16:0/18:1(9Z)": [], # PE O-34:2
                "PE O-18:1(11E);5OMe/22:0;3OH": [], # PE O-41:1;O2
                "M(IP)2C(1) 20:0;3OH,4OH/26:0;2OH": [], # M(IP)2C 46:0;O4
                "SPB 18:0;3OH": [], # SPB 18:0;O
                "SPB 18:0;1OH;3oxo": [], # SPB 18:1;O2
                "PIP(3') 16:0/18:1(9Z)": [], # PIP 34:1
                "Cer 18:0;1OH,3OH,4OH/26:0;2OH,3OH": [], # Cer 44:0;O5
                "Cer 18:1(5Z);1OH,3OH/14:0": [], # Cer 32:1;O2 / C32H63NO3
                "PIP2(3',5') 17:0/20:4(5Z,8Z,11Z,14Z)": [], # PIP2 37:4 / C46H83O19P3
                "PE O-16:0/18:2(9Z,12Z)": [], # PE O-34:2 / C39H76NO7P
                "PS-N(6:0) 16:0/18:3(9Z,12Z,15Z)": [], # PS-N(Alk) 40:3
                "TG 16:0;5O(FA 16:0[R])/18:1(9Z)/18:1(9Z)": [], # TG 68:3;O2
                "TG O-18:1(9Z)/O-16:0/O-18:1(9Z)": [], # TG tO-52:2
                "FA 22:4(4Z,7Z,10Z,18E);[13-17cy5;14OH,16OH];20OH": [], # FA 22:5;O3 / C22H34O5
                "FA 20:2(5Z,13E);[8-13cy6;9OH,11OH;11oxy];15OH": [], # FA 20:3;O4 / C20H34O6
                "FA 20:4(6Z,8E,10E,14Z);5OH,12OH": [], # FA 20:4;O2 / C20H32O4
                "FA 22:0;4O(FA 10:0),5O(FA 10:0)": [], #  42:0;O4
                "FA 22:0;4O(FA 10:0);5O(10:0)": [] #  42:0;O4
                }
        """

        for lipid_name in data:
            print(lipid_name)
            
            results = data[lipid_name]
            lipid = parser.parse(lipid_name)
            formula = results[4] if len(results) > 4 else lipid.get_sum_formula()
            
            if len(results) > 4:
                formula = results[4]
                self.assertEqual(formula, lipid.get_sum_formula())
            else:
                formula = lipid.get_sum_formula()
            

            for l, lipid_level in enumerate([LipidLevel.ISOMERIC_SUBSPECIES, LipidLevel.STRUCTURAL_SUBSPECIES, LipidLevel.MOLECULAR_SUBSPECIES, LipidLevel.SPECIES]):
                n = lipid.get_lipid_string(lipid_level)
                self.assertEqual(results[l], n)
                self.assertEqual(formula, lipid.get_sum_formula())


                lipid2 = parser.parse(n)
                self.assertEqual(results[l], lipid2.get_lipid_string())
                self.assertEqual(formula, lipid2.get_sum_formula())
        
