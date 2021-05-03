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
        
        self.info = LipidSpeciesInfo(self.lipid_class)
        self.info.level = LipidLevel.SPECIES
        if self.lipid_category == LipidCategory.SP and all_lipids[self.lipid_class]["name"] not in {"Cer", "SPB"}:
            if "OH" not in self.info.functional_groups: self.info.functional_groups["OH"] = []
            self.info.functional_groups["OH"].append(get_functional_group("OH").copy())
        
        for fas in fa: self.info.add(fas)
        
        
        
    def add_decorator(self, decorator):
        self.headgroup_decorators.append(decorator)
        if decorator.name in {"decorator_alkyl", "decorator_acyl"}:
            if decorator.name not in self.info.functional_groups: self.info.functional_groups[decorator.name] = []
            self.info.functional_groups[decorator.name].append(decorator)
            
            self.info.num_carbon += decorator.get_elements()[Element.C]
            self.info.double_bonds += decorator.get_double_bonds()
        
        
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
                
            prefix = sorted(hgd.to_string(level) for hgd in self.headgroup_decorators if not hgd.suffix)
            suffix = [hgd.to_string(level) for hgd in self.headgroup_decorators if hgd.suffix]
            lipid_string = prefix + lipid_string + suffix
            
            
            
            if self.info != None and (self.info.elements[Element.C] > 0 or self.info.num_carbon > 0):
                lipid_string += " " if all_lipids[self.lipid_class]["category"] != LipidCategory.ST else "/"
                lipid_string += self.info.to_string()
                
            return "".join(lipid_string)
        
        else:
            raise RuntimeException("LipidSpecies does not know how to create a lipid string for level %s" + level)
        
        
        
    def get_elements_headgroup(self):
        if self.use_head_group or self.info.level not in {LipidLevel.STRUCTURAL_SUBSPECIES, LipidLevel.ISOMERIC_SUBSPECIES, LipidLevel.SPECIES, LipidLevel.MOLECULAR_SUBSPECIES}:
            raise LipidException("Element table cannot be computed for lipid level '%s'" % self.info.level)
        
        dummy = FunctionalGroup("dummy", elements = all_lipids[self.lipid_class]["elements"])
        for hgd in self.headgroup_decorators: dummy += hgd
        
        return dummy.elements
    
        
        
    def get_elements(self):
        
        dummy = FunctionalGroup("dummy", elements = self.get_elements_headgroup())
        dummy += self.info
        
        # since only one FA info is provided, we have to treat this single information as
        # if we would have the complete information about all possible FAs in that lipid
        additional_fa = all_lipids[self.lipid_class]["poss_fa"]
        is_sp = self.lipid_category == LipidCategory.SP
        dummy.elements[Element.O] += additional_fa - max(0, self.info.num_ethers) - (self.lipid_category == LipidCategory.SP)
        dummy.elements[Element.H] -= additional_fa - 2 * max(0, self.info.num_ethers) + (self.lipid_category == LipidCategory.SP)
        
        return dummy.elements
        
