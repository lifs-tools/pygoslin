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
from pygoslin.parser.Parser import *
from pygoslin.tests.pysmiles import read_smiles
import logging

logging.getLogger('pygoslin.tests.smiles_helper').setLevel(logging.CRITICAL) 
logging.getLogger('pygoslin.tests.pysmiles').setLevel(logging.CRITICAL) 



def clean_SMILES(smiles):
    pos = smiles.find("[")
    while pos > -1:
        end = smiles.find("]", pos)
        smiles = smiles[:pos] + smiles[pos + 1] + smiles[end + 1:]
        pos = smiles.find("[", pos + 1)
    
    return smiles



"""
parser = SwissLipidsParser()
lipid = parser.parse("PA(P-20:1(11Z)/28:6(10Z,13Z,16Z,19Z,22Z,25Z))")
print(lipid.get_lipid_string())
print(lipid.get_smiles())
exit()
"""




class SMILESTest(unittest.TestCase):
    
    def test_goodSMILES(self):
        parser = SwissLipidsParser()
        
        with open("pygoslin/tests/swiss-smiles.csv") as infile:
            for i, line in enumerate(infile):
                if i > 0 and i % 1000 == 0: print(i)
                
                tokens = line.strip().split("\t")
                #print(tokens[1])
                
                lipid = parser.parse(tokens[1])
                smiles = clean_SMILES(lipid.get_smiles())
                
                mol = read_smiles(smiles)
                
                
                
                edge_count = {n: 0 for n in mol.nodes()}
                for s, e in mol.edges():
                    edge_count[s] += 1
                    edge_count[e] += 1
                    
                fingerprint = {}
                for n in mol.nodes():
                    el = mol.nodes()[n]
                    ed = edge_count[n]
                    key = "%s-%i" % (el, ed)
                    if key not in fingerprint: fingerprint[key] = 0
                    fingerprint[key] += 1
                
                smiles_swiss = clean_SMILES(tokens[2])
                mol_swiss = read_smiles(smiles_swiss)
                
                edge_count = {n: 0 for n in mol_swiss.nodes()}
                for s, e in mol_swiss.edges():
                    edge_count[s] += 1
                    edge_count[e] += 1
                    
                fingerprint_swiss = {}
                for n in mol_swiss.nodes():
                    el = mol_swiss.nodes()[n]
                    ed = edge_count[n]
                    key = "%s-%i" % (el, ed)
                    if key not in fingerprint_swiss: fingerprint_swiss[key] = 0
                    fingerprint_swiss[key] += 1
                    
                for key, num in fingerprint.items():
                    if key not in fingerprint_swiss or fingerprint_swiss[key] != num:
                        print("Error on lipid %s (%s, %i):" % (tokens[1], tokens[0], i))
                        print("Computed SMILES: %s" % lipid.get_smiles())
                        print("SwissLipids SMILES: %s" % tokens[2])
                        break
                        #raise Exception()
        
            
    """
    @unittest.expectedFailure
    def test_wrong_bond_type(self):
        instanceZero = FattyAcid("FA1", 2, -1, 0, LipidFaBondType.UNDEFINED, False, 0)
    """
    
if __name__ == '__main__':
    unittest.main()
