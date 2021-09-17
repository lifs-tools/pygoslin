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
from pygoslin.domain.Element import *
from pygoslin.domain.LipidClass import all_lipids
from pygoslin.domain.FunctionalGroup import FunctionalGroup

ether_prefix = ["", "O-", "dO-", "tO-", "eO-"] 
        
class LipidSpeciesInfo(FattyAcid):
    
    def __init__(self, lipid_class):
        super().__init__("info")
        
        self.level = None
        self.num_ethers = 0
        self.num_specified_fa = 0
        self.total_fa = all_lipids[lipid_class]["max_fa"]
        self.extended_class = LipidFaBondType.ESTER
        
        
    def add(self, fa):
        
        if fa.lipid_FA_bond_type in {LipidFaBondType.ETHER_PLASMENYL, LipidFaBondType.ETHER_PLASMANYL} and fa.lipid_FA_bond_type not in {LipidFaBondType.LCB_REGULAR, LipidFaBondType.LCB_EXCEPTION}:
            self.num_ethers += 1
            self.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMANYL
            self.extended_class = fa.lipid_FA_bond_type
            
        elif fa.lipid_FA_bond_type in {LipidFaBondType.LCB_REGULAR, LipidFaBondType.LCB_EXCEPTION}:
            self.lipid_FA_bond_type = fa.lipid_FA_bond_type
        
                
        else:
            self.num_specified_fa += 1
        
        for fg, fg_list in fa.functional_groups.items():
            if fg not in self.functional_groups: self.functional_groups[fg] = []
            for func_group in fg_list:
                self.functional_groups[fg].append(func_group)
         
        self.num_carbon += fa.get_elements()[Element.C]
        self.double_bonds += fa.get_double_bonds()
        
        
        
    def get_elements(self):
        elements = super().get_elements()
        if self.lipid_FA_bond_type != LipidFaBondType.LCB_EXCEPTION:
            elements[Element.O] -= (self.num_ethers == 0)
        elements[Element.H] += 1 if self.num_ethers == 0 else -1
        
        return elements
        
    
    def to_string(self):
        global ether_prefix
        
        info_string = [ether_prefix[self.num_ethers]]
        info_string.append("%i:%i" % (self.num_carbon, self.double_bonds))
        

        elements = self.get_functional_group_elements()
        additional_elements = ";".join("%s%s" % (element_shortcut[e], str(elements[e]) if elements[e] > 1 else "") for e in element_order[2:] if e in elements and elements[e] > 0)
        
        if len(additional_elements) > 0:
            info_string.append(";")
            info_string.append(additional_elements)
            
        return "".join(info_string)
