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


from pygoslin.domain.LipidFullStructure import LipidFullStructure
from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.LipidLevel import LipidLevel

class LipidCompleteStructure(LipidFullStructure):


    def __init__(self, head_group, fa = []):
        super().__init__(head_group, fa)
                
        self.info.level = LipidLevel.COMPLETE_STRUCTURE
        
    
    def get_extended_class(self):
        return super().get_extended_class()
    

    def get_lipid_string(self, level = None):
        
        if level == None or level == LipidLevel.COMPLETE_STRUCTURE:
            return self.build_lipid_subspecies_name(LipidLevel.COMPLETE_STRUCTURE)
        
        elif level in (LipidLevel.FULL_STRUCTURE, LipidLevel.STRUCTURE_DEFINED, LipidLevel.SN_POSITION, LipidLevel.COMPLETE_STRUCTURE, LipidLevel.MOLECULAR_SPECIES, LipidLevel.CATEGORY, LipidLevel.CLASS, LipidLevel.SPECIES):
            return super().get_lipid_string(level)
        
        else:
            raise Exception("LipidIsomericSubspecies does not know how to create a lipid string for level %s" % level)
