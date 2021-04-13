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

from pygoslin.domain.FattyAcid import FattyAcid
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.Element import Element

class LipidSpeciesInfo(FattyAcid):
    
    def __init__(self, fa = None):
        if fa != None:
            super().__init__(fa.name, fa.num_carbon, fa.double_bonds, fa.functional_groups, fa.lipid_FA_bond_type, fa.lcb, fa.position)
        else:
            super().__init__("", 2, 0, {}, LipidFaBondType.ESTER, False, 0)
        
        self.level = None
        self.num_oxygen = 0
        self.num_esters = 0
        self.ester_prefix = ["", "O-", "dO-", "tO-", "eO-"]
        
    def add(self, fa):
        self.num_carbon += fa.num_carbon
        self.double_bonds += fa.double_bonds if type(fa.double_bonds) == int else len(fa.double_bonds)
        if fa.lipid_FA_bond_type in {LipidFaBondType.ETHER_PLASMENYL, LipidFaBondType.ETHER_PLASMANYL}:
            self.num_esters += 1
            
            if self.lipid_FA_bond_type in {LipidFaBondType.ESTER, LipidFaBondType.ETHER_PLASMENYL}:
                self.lipid_FA_bond_type = fa.lipid_FA_bond_type
        
        self.num_oxygen += fa.get_num_oxygens()
        
    
    def to_string(self):
        info_string = [self.ester_prefix[self.num_esters]]
        info_string.append("%i:%i" % (self.num_carbon, self.double_bonds))
        if self.num_oxygen > 0:
            info_string.append(";O%s" % (str(self.num_oxygen) if self.num_oxygen > 1 else ""))
            
        return "".join(info_string)
