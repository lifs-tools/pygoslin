"""
MIT License

Copyright (c) 2020 Dominik Kopczynski   -   dominik.kopczynski {at} univie.ac.at
                   Nils Hoffmann  -  nils.hoffmann {at} cebitec.uni-bielefeld.de

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
from pygoslin.domain.LipidSpeciesInfo import LipidSpeciesInfo
from pygoslin.domain.LipidClass import *
from pygoslin.domain.FunctionalGroup import *
from pygoslin.domain.LipidFaBondType import *
from pygoslin.domain.Element import Element

glyco_table = {"ga1": ["Gal", "GalNAc", "Gal", "Glc"],
               "ga2": ["GalNAc", "Gal", "Glc"],
               "gb3": ["Gal", "Gal", "Glc"],
               "gb4": ["GalNAc", "Gal", "Gal", "Glc"],
               "gd1": ["Gal", "GalNAc", "NeuAc", "NeuAc", "Gal", "Glc"],
               "gd1a": ["Hex", "Hex", "Hex", "HexNAc", "NeuAc", "NeuAc"],
               "gd2": ["GalNAc", "NeuAc", "NeuAc", "Gal", "Glc"],
               "gd3": ["NeuAc", "NeuAc", "Gal", "Glc"],
               "gm1": ["Gal", "GalNAc", "NeuAc", "Gal", "Glc"],
               "gm2": ["GalNAc", "NeuAc", "Gal", "Glc"],
               "gm3": ["NeuAc", "Gal", "Glc"],
               "gm4": ["NeuAc", "Gal"],
               "gp1": ["NeuAc", "NeuAc", "Gal", "GalNAc", "NeuAc", "NeuAc", "NeuAc", "Gal", "Glc"],
               "gq1": ["NeuAc", "Gal", "GalNAc", "NeuAc", "NeuAc", "NeuAc", "Gal", "Glc"],
               "gt1": ["Gal", "GalNAc", "NeuAc", "NeuAc", "NeuAc", "Gal", "Glc"],
               "gt2": ["GalNAc", "NeuAc", "NeuAc", "NeuAc", "Gal", "Glc"],
               "gt3": ["NeuAc", "NeuAc", "NeuAc", "Gal", "Glc"]
               }

class HeadGroup:
    def __init__(self, headgroup, decorators = None, use_headgroup = False):
        self.decorators = [d for d in decorators] if decorators != None else []
        self.parsed_headgroup = get_class(headgroup)
        if self.parsed_headgroup == UNDEFINED_LIPID_CLASS:
            self.parsed_headgroup = headgroup
        else:
            self.parsed_headgroup = all_lipids[self.parsed_headgroup]["name"]
        
        # checking if head group is a glyco-sphingolipid
        hg = headgroup.strip(" ").lower()
        if hg in glyco_table and not use_headgroup:
            for carbohydrate in glyco_table[hg]:
                try:
                    
                    functional_group = get_functional_group(carbohydrate)
                    functional_group.elements[Element.O] -= 1
                    self.decorators.append(functional_group)
                except Exception:
                    raise LipidParsingException("Carbohydrate '%s' unknown" % carbohydrate)
            headgroup = "Cer"
        
        self.headgroup = headgroup.strip(" ")
        self.lipid_category = get_category(self.headgroup)
        self.lipid_class = get_class(self.headgroup)
        self.use_headgroup = use_headgroup
        self.sp_exception = self.lipid_category == LipidCategory.SP and all_lipids[self.lipid_class]["name"] in {"Cer", "SPB"} and len(self.decorators) == 0
        
        
    def get_lipid_string(self, level = None):
        if level == LipidLevel.CATEGORY:
            return self.lipid_category.name
        
        headgoup_string = [all_lipids[self.lipid_class]["name"] if not self.use_headgroup else self.headgroup]
        prefix = []
                
        if len(self.decorators) > 0:
            if level not in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE, LipidLevel.STRUCTURE_DEFINED}:
                decorators_sorted = []
                for hgd in self.decorators:
                    if hgd.suffix: continue
                    
                    hgd_copy = hgd.copy()
                    if hgd_copy.name.find("Gal") > -1 or hgd_copy.name.find("Glc") > -1 or hgd_copy.name.find("S(3')") > -1:
                        hgd_copy.name = hgd_copy.name.replace("Gal", "Hex").replace("Glc", "Hex").replace("S(3')", "S")
                    decorators_sorted.append(hgd_copy)
                
                if len(decorators_sorted) > 0:
                    decorators_sorted.sort(key = lambda x: x.name)
                    name, rep = "", 0
                    
                    for hgd in decorators_sorted:
                        if name != hgd.name:
                            if rep > 1: prefix.append("%i" % rep)
                            name, rep = hgd.name, hgd.count
                            prefix.append(name)
                        else:
                            rep += hgd.count
                        
                    if rep > 1: prefix.append("%i" % rep)
                
            else:
                prefix = ["%s-" % hgd.to_string(level) for hgd in self.decorators if not hgd.suffix]
            
        suffix = [hgd.to_string(level) for hgd in self.decorators if hgd.suffix]
        if level in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE} and self.lipid_category == LipidCategory.SP and not self.sp_exception: suffix.append("(1)")
            
        return "".join(prefix + headgoup_string + suffix)
        
        
        
    def get_elements(self):
        if self.use_headgroup:
            raise LipidException("Element table cannot be computed for lipid '%s'" % self.headgroup)
        
        dummy = FunctionalGroup("dummy", elements = all_lipids[self.lipid_class]["elements"])
        for hgd in self.decorators: dummy += hgd
        
        return dummy.elements
