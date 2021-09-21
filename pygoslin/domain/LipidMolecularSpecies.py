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
from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.FunctionalGroup import FunctionalGroup
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.LipidClass import *
from pygoslin.domain.FattyAcid import FattyAcid

class LipidMolecularSpecies(LipidSpecies):


    def __init__(self, head_group, fa = []):
        super().__init__(head_group)
        self.fa = {}
        self.fa_list = []
        
        self.info.level = LipidLevel.MOLECULAR_SPECIES
        
        for fas in fa:
            if fas.name in self.fa:
                raise ConstraintViolationException("FA names must be unique! FA with name %s was already added!" % fas.name)
            
            else:
                self.fa[fas.name] = fas
                self.fa_list.append(fas)
                self.info.add(fas)
                
                
        # add 0:0 dummys
        for i in range(len(fa), self.info.total_fa):
            fas = FattyAcid("FA%i" % (i + len(fa) + 1))
            self.info.add(fas)
            self.fa[fas.name] = fas
            self.fa_list.append(fas)
    

    def get_extended_class(self):
        return super().get_extended_class()
    
    

    def build_lipid_subspecies_name(self, level):
        if level == None: level = LipidLevel.MOLECULAR_SPECIES
        fa_separator = "/" if level != LipidLevel.MOLECULAR_SPECIES or all_lipids[self.headgroup.lipid_class]["category"] == LipidCategory.SP else "_"

        lipid_name = [self.headgroup.get_lipid_string(level)]

        fa_headgroup_separator = " " if all_lipids[self.headgroup.lipid_class]["category"] != LipidCategory.ST else "/"
        
        if level in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE, LipidLevel.STRUCTURE_DEFINED, LipidLevel.SN_POSITION}:
            fa_string = fa_separator.join(fatty_acid.to_string(level) for fatty_acid in self.fa_list)
            if len(fa_string) > 0: lipid_name += [fa_headgroup_separator, fa_string]
        else:
            fa_string = fa_separator.join(fatty_acid.to_string(level) for fatty_acid in self.fa_list if fatty_acid.num_carbon > 0)
            if len(fa_string) > 0: lipid_name += [fa_headgroup_separator, fa_string]
            
        return "".join(lipid_name)
    
    
    
    def get_elements(self):
        dummy = FunctionalGroup("headgroup", elements = self.headgroup.get_elements())
        
        # add elements from all fatty acyl chains
        for fa in self.fa_list: dummy += fa
        
        return dummy.elements
    
    
    def get_lipid_string(self, level = None):
        if level == None or level == LipidLevel.MOLECULAR_SPECIES:
            return self.build_lipid_subspecies_name(LipidLevel.MOLECULAR_SPECIES)
        
        elif level in (LipidLevel.CATEGORY, LipidLevel.CLASS, LipidLevel.SPECIES):
            return super().get_lipid_string(level)
        else:
            raise Exception("LipidMolecularSpecies does not know how to create a lipid string for level %s" % level)
    
    