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
import os
from pygoslin.parser.Parser import *
import time
import csv

import unittest


class ShorthandTest(unittest.TestCase):
    
    def test_sphingolipids(self):
        parser = ShorthandParser()
        
        l = parser.parse("Cer 18:1(8Z);1OH,3OH/24:0")
        self.assertEqual(l.get_lipid_string(), "Cer 18:1(8Z);1OH,3OH/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES), "Cer 18:1(8);(OH)2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES), "Cer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "Cer 42:1;O2")
        self.assertEqual(l.get_sum_formula(), "C42H83NO3")
        
        l = parser.parse("Cer 18:1(8);(OH)2/24:0")
        self.assertEqual(l.get_lipid_string(), "Cer 18:1(8);(OH)2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES), "Cer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "Cer 42:1;O2")
        self.assertEqual(l.get_sum_formula(), "C42H83NO3")
        
        l = parser.parse("Cer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(), "Cer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES), "Cer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "Cer 42:1;O2")
        self.assertEqual(l.get_sum_formula(), "C42H83NO3")
        
        l = parser.parse("Cer 42:1;O2")
        self.assertEqual(l.get_lipid_string(), "Cer 42:1;O2")
        self.assertEqual(l.get_sum_formula(), "C42H83NO3")
        
        
        
        l = parser.parse("Gal-Cer(1) 18:1(5Z);3OH/24:0")
        self.assertEqual(l.get_lipid_string(), "Gal-Cer(1) 18:1(5Z);3OH/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES), "Gal-Cer 18:1(5);OH/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES), "GalCer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "GalCer 42:1;O2")
        self.assertEqual(l.get_sum_formula(), "C48H93NO8")
        
        l = parser.parse("Gal-Cer 18:1(5);OH/24:0")
        self.assertEqual(l.get_lipid_string(), "Gal-Cer 18:1(5);OH/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES), "GalCer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "GalCer 42:1;O2")
        self.assertEqual(l.get_sum_formula(), "C48H93NO8")
        
        l = parser.parse("GalCer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(), "GalCer 18:1;O2/24:0")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "GalCer 42:1;O2")
        self.assertEqual(l.get_sum_formula(), "C48H93NO8")
        
        l = parser.parse("GalCer 42:1;O2")
        self.assertEqual(l.get_lipid_string(), "GalCer 42:1;O2")
        self.assertEqual(l.get_sum_formula(), "C48H93NO8")
        
        
        
        l = parser.parse("SPB 18:1(4Z);1OH,3OH")
        self.assertEqual(l.get_lipid_string(), "SPB 18:1(4Z);1OH,3OH")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES), "SPB 18:1(4);(OH)2")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES), "SPB 18:1;O2")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "SPB 18:1;O2")
        self.assertEqual(l.get_sum_formula(), "C18H37NO2")
        
        l = parser.parse("SPB 18:1(4);(OH)2")
        self.assertEqual(l.get_lipid_string(), "SPB 18:1(4);(OH)2")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES), "SPB 18:1;O2")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "SPB 18:1;O2")
        self.assertEqual(l.get_sum_formula(), "C18H37NO2")
        
        l = parser.parse("SPB 18:1;O2")
        self.assertEqual(l.get_lipid_string(), "SPB 18:1;O2")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "SPB 18:1;O2")
        self.assertEqual(l.get_sum_formula(), "C18H37NO2")
        

        
        l = parser.parse("LSM(1) 17:1(4E);3OH")
        self.assertEqual(l.get_lipid_string(), "LSM(1) 17:1(4E);3OH")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES), "LSM 17:1(4);OH")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES), "LSM 17:1;O2")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "LSM 17:1;O2")
        self.assertEqual(l.get_sum_formula(), "C22H47N2O5P")
        
        l = parser.parse("LSM 17:1(4);OH")
        self.assertEqual(l.get_lipid_string(), "LSM 17:1(4);OH")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES), "LSM 17:1;O2")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "LSM 17:1;O2")
        self.assertEqual(l.get_sum_formula(), "C22H47N2O5P")
        
        l = parser.parse("LSM 17:1;O2")
        self.assertEqual(l.get_lipid_string(), "LSM 17:1;O2")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "LSM 17:1;O2")
        self.assertEqual(l.get_sum_formula(), "C22H47N2O5P")
        

        
        l = parser.parse("EPC(1) 14:1(4E);3OH/20:1(11Z)")
        self.assertEqual(l.get_lipid_string(), "EPC(1) 14:1(4E);3OH/20:1(11Z)")
        self.assertEqual(l.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES), "EPC 14:1(4);OH/20:1(11)")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES), "EPC 14:1;O2/20:1")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "EPC 34:2;O2")
        self.assertEqual(l.get_sum_formula(), "C36H71N2O6P")
        
        l = parser.parse("EPC 14:1(4);OH/20:1(11)")
        self.assertEqual(l.get_lipid_string(), "EPC 14:1(4);OH/20:1(11)")
        self.assertEqual(l.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES), "EPC 14:1;O2/20:1")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "EPC 34:2;O2")
        self.assertEqual(l.get_sum_formula(), "C36H71N2O6P")
        
        l = parser.parse("EPC 14:1;O2/20:1")
        self.assertEqual(l.get_lipid_string(), "EPC 14:1;O2/20:1")
        self.assertEqual(l.get_lipid_string(LipidLevel.SPECIES), "EPC 34:2;O2")
        self.assertEqual(l.get_sum_formula(), "C36H71N2O6P")
        
        l = parser.parse("EPC 34:2;O2")
        self.assertEqual(l.get_lipid_string(), "EPC 34:2;O2")
        self.assertEqual(l.get_sum_formula(), "C36H71N2O6P")

        
        
    def test_nomenclature(self):
        data = {}
        file_name = os.path.join("pygoslin", "data", "goslin", "testfiles", "shorthand-test.csv")
        with open(file_name, mode = "rt") as infile:
            lipidreader = csv.reader(infile, delimiter=',', quotechar='"')
            for row in lipidreader:
                row[-1] = row[-1].strip(" ")
                data[row[0]] = row[1:]
                
        parser = ShorthandParser()

        for lipid_name in data:
            
            results = data[lipid_name]
            lipid = parser.parse(lipid_name)
            formula = results[4] if len(results) > 4 else lipid.get_sum_formula()
            
            if len(results) > 4 and len(results[4]) > 0:
                formula = results[4]
                self.assertEqual(formula, lipid.get_sum_formula(), "test on lipid '%s'" % lipid_name)
            else:
                formula = lipid.get_sum_formula()
            
            levels = [LipidLevel.ISOMERIC_SUBSPECIES, LipidLevel.STRUCTURAL_SUBSPECIES, LipidLevel.MOLECULAR_SUBSPECIES, LipidLevel.SPECIES]
            for l, lipid_level in enumerate(levels):
                n = lipid.get_lipid_string(lipid_level)
                self.assertEqual(results[l], n, "test on lipid '%s' and level '%s'" % (lipid_name, lipid_level))
                self.assertEqual(formula, lipid.get_sum_formula(), "test on lipid '%s' and level '%s'" % (lipid_name, lipid_level))


                lipid2 = parser.parse(n)
                for ll in range(l, 4):
                    self.assertEqual(results[ll], lipid2.get_lipid_string(levels[ll]), "test on lipid '%s' and level '%s'" % (lipid_name, lipid_level))
                    self.assertEqual(formula, lipid2.get_sum_formula(), "test on lipid '%s' and level '%s'" % (lipid_name, lipid_level))
        
        
        
    def teest_performance(self):
        cycles, parser = 50, ShorthandParser()
        length = 0
        
        start = time.time()
        for lipid_name in data:
            for i in range(cycles):
                lipid = parser.parse(lipid_name)
                length += len(lipid_name)
                
                lipid_name2 = lipid.get_lipid_string(LipidLevel.STRUCTURAL_SUBSPECIES)
                lipid = parser.parse(lipid_name2)
                length += len(lipid_name2)
                
                lipid_name2 = lipid.get_lipid_string(LipidLevel.MOLECULAR_SUBSPECIES)
                lipid = parser.parse(lipid_name2)
                length += len(lipid_name2)
                
                lipid_name2 = lipid.get_lipid_string(LipidLevel.SPECIES)
                lipid = parser.parse(lipid_name2)
                length += len(lipid_name2)
                
        elapsed = time.time() - start
        print("time elapsed: %f" % elapsed)
        print("lipid / second: %f" % (len(data) * cycles * 4 / elapsed))
        print("Avg. lipid name length: %f" % (length / (len(data) * cycles * 4)))
