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


from pygoslin.domain.LipidSpecies import LipidSpecies
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.LipidSpeciesInfo import LipidSpeciesInfo
from pygoslin.domain.FunctionalGroup import FunctionalGroup
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.LipidClass import *

class LipidMolecularSubspecies(LipidSpecies):


    def __init__(self, head_group, fa = []):
        super().__init__(head_group)
        self.fa = {}
        self.fa_list = []
        self.info = LipidSpeciesInfo()
        self.info.level = LipidLevel.MOLECULAR_SUBSPECIES
        
        for fas in fa:
            if fas.name in self.fa:
                raise ConstraintViolationException("FA names must be unique! FA with name %s was already added!" % fas.name)
            
            else:
                self.fa[fas.name] = fas
                self.fa_list.append(fas)
                self.info.add(fas)
    


    def get_extended_class(self):
        return super().get_extended_class()
    
    

    def build_lipid_subspecies_name(self, fa_separator, level):

        fa_headgroup_separator = " " if all_lipids[self.lipid_class]["category"] != LipidCategory.ST else "/"
        
        if level in {LipidLevel.ISOMERIC_SUBSPECIES, LipidLevel.STRUCTURAL_SUBSPECIES}:
            fa_string = fa_separator.join(fatty_acid.to_string(level) for fatty_acid in self.fa_list)
            if len(fa_string) > 0: fa_string = fa_headgroup_separator + fa_string
        else:
            fa_string = fa_separator.join(fatty_acid.to_string(level) for fatty_acid in self.fa_list if fatty_acid.num_carbon > 0)
            if len(fa_string) > 0: fa_string = fa_headgroup_separator + fa_string
            
        head_group = [(all_lipids[self.lipid_class]["name"] if not self.use_head_group else self.head_group)]
        
        
        if level == LipidLevel.ISOMERIC_SUBSPECIES and self.lipid_category == LipidCategory.SP and head_group[0] not in {"Cer", "SPB"}:
            head_group.append("(1)")
            
        if level == LipidLevel.ISOMERIC_SUBSPECIES:
            head_group = ["%s-" % hgd.to_string(level) for hgd in self.headgroup_decorators] + head_group
                
        else:
            head_group = sorted(hgd.to_string(level) for hgd in self.headgroup_decorators) + head_group
        
        return "".join(head_group) + fa_string
    
    
    
    def get_elements(self):
        dummy = FunctionalGroup("dummy", elements = super().get_elements()) # get elements from head group + all decorators
        # add elements from all fatty acyl chains
        
        head_group = (all_lipids[self.lipid_class]["name"] if not self.use_head_group else self.head_group)
        if self.lipid_category == LipidCategory.SP and head_group == "Cer" and len(self.headgroup_decorators) == 0:
            dummy.elements[Element.O] -= 1
        
        for fa in self.fa_list:
            fa.compute_elements()
            dummy += fa
        
        return dummy.elements
    
    
    def get_lipid_string(self, level = None):
        if level == None or level == LipidLevel.MOLECULAR_SUBSPECIES:
            return self.build_lipid_subspecies_name("_", LipidLevel.MOLECULAR_SUBSPECIES)
        
        elif level in (LipidLevel.CATEGORY, LipidLevel.CLASS, LipidLevel.SPECIES):
            return super().get_lipid_string(level)
        else:
            raise Exception("LipidMolecularSubspecies does not know how to create a lipid string for level %s" % level)
    
    
