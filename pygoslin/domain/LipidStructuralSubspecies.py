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


from pygoslin.domain.LipidMolecularSubspecies import LipidMolecularSubspecies
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidSpeciesInfo import LipidSpeciesInfo
from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.LipidClass import *
from pygoslin.domain.FattyAcid import FattyAcid

class LipidStructuralSubspecies(LipidMolecularSubspecies):


    def __init__(self, head_group, fa = []):
        super().__init__(head_group, fa)
        
        self.info.level = LipidLevel.STRUCTURAL_SUBSPECIES   
    

    def get_extended_class(self):
        return super().get_extended_class()
    
    
    def get_lipid_string(self, level = None):
        if level == None or level == LipidLevel.STRUCTURAL_SUBSPECIES:
            return self.build_lipid_subspecies_name(LipidLevel.STRUCTURAL_SUBSPECIES)
        
        elif level in (LipidLevel.MOLECULAR_SUBSPECIES, LipidLevel.CATEGORY, LipidLevel.CLASS, LipidLevel.SPECIES):
            return super().get_lipid_string(level)
            
        else:
            raise RuntimeException("LipidStructuralSubspecies does not know how to create a lipid string for level %s" % level)
