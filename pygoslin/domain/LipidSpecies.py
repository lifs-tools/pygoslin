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


from pygoslin.domain.LipidExceptions import RuntimeException
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.LipidCategory import LipidCategory
from pygoslin.domain.FunctionalGroup import *
from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.LipidSpeciesInfo import LipidSpeciesInfo
from pygoslin.domain.LipidClass import *
from pygoslin.domain.LipidFaBondType import *
from pygoslin.domain.Element import Element

class LipidSpecies:
    def __init__(self, head_group, fa = []):
        self.head_group = head_group.strip(" ")
        self.lipid_category = get_category(self.head_group)
        self.lipid_class = get_class(self.head_group)
        self.use_head_group = False
        self.headgroup_decorators = []
        
        self.info = LipidSpeciesInfo()
        self.info.level = LipidLevel.SPECIES
    
        for fas in fa: self.info.add(fas)
        
        
    def clone(self, fa):
        self.head_group = fa.head_group
        self.lipid_category = fa.lipid_category
        self.lipid_class = fa.lipid_class
        self.info = LipidSpeciesInfo(fa)
        self.use_head_group = fa.use_head_group
        
        
        
    def get_extended_class(self):
        special_case = self.lipid_category == LipidCategory.GP if self.info != None and self.info.num_carbon > 0 else False
        if special_case and self.info.lipid_FA_bond_type in [LipidFaBondType.ETHER_PLASMANYL, LipidFaBondType.ETHER_UNSPECIFIED]:
            return all_lipids[self.lipid_class]["name"] + "-O"
        
        if special_case and self.info.lipid_FA_bond_type == LipidFaBondType.ETHER_PLASMENYL:
            return all_lipids[self.lipid_class]["name"] + "-p"
        
        return all_lipids[self.lipid_class]["name"]
            


    def get_lipid_string(self, level = None):
        if level == LipidLevel.CATEGORY:
            return self.lipid_category.name
        
        elif level == LipidLevel.CLASS:
            return all_lipids[self.lipid_class]["name"] if not self.use_head_group else self.head_group
        
        elif level == None or level == LipidLevel.SPECIES:
            lipid_string = [all_lipids[self.lipid_class]["name"] if not self.use_head_group else self.head_group]
            if self.info != None and self.info.num_carbon > 0:
                #special_case = self.lipid_class in self.special_cases
                special_case = self.lipid_category == LipidCategory.GP
                
                lipid_string += " " if all_lipids[self.lipid_class]["category"] != LipidCategory.ST else "/"
                lipid_string += self.info.to_string()
            return "".join(lipid_string)
        
        else:
            raise RuntimeException("LipidSpecies does not know how to create a lipid string for level %s" + level)
        
        
        
    def get_elements(self):
        if self.use_head_group or self.info.level not in {LipidLevel.STRUCTURAL_SUBSPECIES, LipidLevel.ISOMERIC_SUBSPECIES}:
            raise LipidException("Element table cannot be computed for lipid level '%s'" % self.info.level)
        
        dummy = FunctionalGroup("dummy", elements = all_lipids[self.lipid_class]["elements"])
        for hgd in self.headgroup_decorators: dummy += hgd
        
        return dummy.elements
        
