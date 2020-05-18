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
from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.LipidClass import *
from pygoslin.domain.Element import Element

class LipidSpecies:
    def __init__(self, head_group, lipid_category = None, lipid_class = None, lipid_species_info = None):
        self.head_group = head_group.strip(" ")
        self.lipid_category = lipid_category if lipid_category != None else get_category(self.head_group)
        
        self.lipid_class = lipid_class if lipid_class != None else get_class(self.head_group)
        self.info = lipid_species_info
        self.use_head_group = False
        self.special_cases = {class_string_to_class["PC"], class_string_to_class["PE"], class_string_to_class["LPC"], class_string_to_class["LPE"]}
        
        
    def clone(self, fa):
        self.head_group = fa.head_group
        self.lipid_category = fa.lipid_category
        self.lipid_class = fa.lipid_class
        self.info = LipidSpeciesInfo(fa)
        self.use_head_group = fa.use_head_group
            
            
    def validate(self):
        return True
    
        """
        if self.use_head_group: return True
        if len(all_lipids) <= self.lipid_class: return False
        return all_lipids[self.lipid_class][3] == 0 or (all_lipids[self.lipid_class][3] > 0 and self.info.num_carbon >= 2)
        """


    def get_lipid_string(self, level = None):
        
        if level == None:
            if self.info != None:
                level = self.info.level
            else:
                raise RuntimeException("LipidSpecies does not know how to create a lipid string for level %s" % (level if level != None else " unknown"))
        
        if level == LipidLevel.CATEGORY:
            return self.lipid_category.name
        
        elif level == LipidLevel.CLASS:
            return all_lipids[self.lipid_class]["name"] if not self.use_head_group else self.head_group
        
        elif level == LipidLevel.SPECIES:
            if not self.validate(): raise ConstraintViolationException("No fatty acly chain information present for lipid '%s'" % all_lipids[self.lipid_class]["name"])
            lipid_string = [all_lipids[self.lipid_class]["name"] if not self.use_head_group else self.head_group]
            if self.info != None and self.info.num_carbon > 0:
                special_case = self.lipid_class in self.special_cases
                
                lipid_string += " " if all_lipids[self.lipid_class]["category"] != LipidCategory.ST else "/"
                lipid_string += self.info.to_string(special_case)
            return "".join(lipid_string)
        
        else:
            raise RuntimeException("LipidSpecies does not know how to create a lipid string for level %s" + level)
        
        
    def get_elements(self):
        if self.use_head_group:
            return {e: 0 for e in Element}
        
        
        hg_elements = all_lipids[self.lipid_class]["elements"]
        try:
            elements = {e: hg_elements[e] for e in Element}
        except:
            raise LipidException("Inconsistant element tables")
            
        if self.info.level in {LipidLevel.MOLECULAR_SUBSPECIES, LipidLevel.STRUCTURAL_SUBSPECIES, LipidLevel.ISOMERIC_SUBSPECIES}:
            num_true_fa = 0
            for fa in self.fa_list:
                fa_elements = fa.get_elements()
                #print("%s: %s" % (self.head_group, fa_elements))
                if fa.num_carbon != 0 or fa.num_double_bonds != 0: num_true_fa += 1
                for e in Element:
                    elements[e] += fa_elements[e]
            if all_lipids[self.lipid_class]["max_fa"] < num_true_fa:
                raise LipidException("Inconsistancy in number of fatty acyl chains for lipid '%s'" % self.head_group)
            elements[Element.H] += all_lipids[self.lipid_class]["max_fa"] - num_true_fa # adding hydrogens for absent fatty acyl chains
            
        
        elif self.info.level == LipidLevel.SPECIES:
            fa_elements = self.info.get_elements(max(all_lipids[self.lipid_class]["poss_fa"]))
            for e in Element:
                elements[e] += fa_elements[e]
            elements[Element.H] += all_lipids[self.lipid_class]["max_fa"] - max(all_lipids[self.lipid_class]["poss_fa"]) # adding hydrogens for absent fatty acyl chains
            
        
        return elements
        
