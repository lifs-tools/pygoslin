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

class HeadGroup:
    def __init__(self, headgroup, decorators = None, use_headgroup = False):
        self.headgroup = headgroup.strip(" ")
        self.lipid_category = get_category(self.headgroup)
        self.lipid_class = get_class(self.headgroup)
        self.use_headgroup = use_headgroup
        self.decorators = [d for d in decorators] if decorators != None else []
        self.sp_exception = self.lipid_category == LipidCategory.SP and all_lipids[self.lipid_class]["name"] in {"Cer", "SPB"} and len(self.decorators) == 0
        
        
    def get_lipid_string(self, level = None):
        if level == LipidLevel.CATEGORY:
            return self.lipid_category.name
        
        elif level == LipidLevel.CLASS:
            return all_lipids[self.lipid_class]["name"] if not self.use_headgroup else self.headgroup
        
        headgoup_string = [all_lipids[self.lipid_class]["name"] if not self.use_headgroup else self.headgroup]
                
        if level not in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE, LipidLevel.STRUCTURE_DEFINED}:
            decorators_tmp = [hgd.copy() for hgd in self.decorators if not hgd.suffix]
            decorators_tmp.sort(key = lambda x: x.name)
            for i in range(len(decorators_tmp) - 1, 0, - 1):
                hge = decorators_tmp[i]
                hge_before = decorators_tmp[i - 1]
                if hge_before.name == hge.name:
                    hge_before.count += hge.count
                    decorators_tmp.pop(i)
            prefix = [hgd.to_string(level) for hgd in decorators_tmp]
            
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
