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


from pygoslin.parser.BaseParserEventHandler import BaseParserEventHandler
from pygoslin.domain.LipidAdduct import LipidAdduct
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.LipidClass import *
from pygoslin.domain.Adduct import Adduct

from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.FattyAcid import FattyAcid
from pygoslin.domain.HeadGroup import HeadGroup
from pygoslin.domain.FunctionalGroup import *

from pygoslin.domain.LipidCompleteStructure import LipidCompleteStructure
from pygoslin.domain.LipidFullStructure import LipidFullStructure
from pygoslin.domain.LipidStructureDefined import LipidStructureDefined
from pygoslin.domain.LipidSnPosition import LipidSnPosition
from pygoslin.domain.LipidMolecularSpecies import LipidMolecularSpecies
from pygoslin.domain.LipidSpecies import LipidSpecies

from pygoslin.domain.LipidExceptions import *

class LipidBaseParserEventHandler(BaseParserEventHandler):
    SP_EXCEPTION_CLASSES = set(all_lipids[get_class("Cer")]["synonyms"]) | set(all_lipids[get_class("SPB")]["synonyms"]) | set(["Cer"]) | set(["SPB"])
    
    def __init__(self):
        super().__init__()
        self.level = LipidLevel.FULL_STRUCTURE
        self.head_group = ""
        self.lcb = None
        self.fa_list = []
        self.current_fa = None
        self.adduct = None
        self.headgroup_decorators = []
        self.use_head_group = False
        
    
    
    def set_lipid_level(self, level):
        self.level = self.level if self.level.value < level.value else level
        
        
        
    def sp_regular_lcb(self):
        return get_category(self.head_group) == LipidCategory.SP and self.current_fa.lipid_FA_bond_type in FattyAcid.LCB_STATES and not (self.head_group in self.SP_EXCEPTION_CLASSES and len(self.headgroup_decorators) == 0)
    
    
    
    def prepare_headgroup_and_checks(self):
        headgroup = HeadGroup(self.head_group, self.headgroup_decorators, self.use_head_group)
        if use_head_group: return headgroup
    
        max_num_fa = all_lipids[headgroup.lipid_class]["max_fa"]
        if max_num_fa != len(self.fa_list): self.set_lipid_level(LipidLevel.STRUCTURE_DEFINED)
    
        true_fa = sum(1 for fa in self.fa_list if fa.num_carbon > 0 or (fa.double_bonds > 0 if type(fa.double_bonds) == int else len(fa.double_bonds)) > 0)
        
        poss_fa = all_lipids[headgroup.lipid_class]["poss_fa"]
        
        # make lyso
        can_be_lyso = "Lyso" in all_lipids[get_class("L" + self.head_group)]["specials"] if get_class("L" + self.head_group) < len(all_lipids) else False
        
        if true_fa + 1 == poss_fa and self.level != LipidLevel.SPECIES and headgroup.lipid_category == LipidCategory.GP and can_be_lyso:
            self.head_group = "L" + self.head_group
            headgroup = HeadGroup(self.head_group, self.headgroup_decorators, self.use_head_group)
            poss_fa = all_lipids[headgroup.lipid_class]["poss_fa"] if headgroup.lipid_class < len(all_lipids) else 0
        
        elif true_fa + 2 == poss_fa and self.level != LipidLevel.SPECIES and headgroup.lipid_category == LipidCategory.GP and self.head_group == "CL":
            self.head_group = "DL" + self.head_group
            headgroup = HeadGroup(self.head_group, self.headgroup_decorators, self.use_head_group)
            poss_fa = all_lipids[headgroup.lipid_class]["poss_fa"] if headgroup.lipid_class < len(all_lipids) else 0

        if self.level == LipidLevel.SPECIES:
            if true_fa == 0 and poss_fa != 0:
                raise ConstraintViolationException("No fatty acyl information lipid class '%s' provided." % headgroup.headgroup)
            
        elif true_fa != poss_fa and self.level in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE, LipidLevel.STRUCTURE_DEFINED, LipidLevel.SN_POSITION}:
            raise ConstraintViolationException("Number of described fatty acyl chains (%i) not allowed for lipid class '%s' (having %i fatty aycl chains)." % (true_fa, headgroup.headgroup, poss_fa))
        
        if "HC" in all_lipids[headgroup.lipid_class]["specials"] and len(self.fa_list) > 0:
            self.fa_list[0].lipid_FA_bond_type = LipidFaBondType.AMINE
        
        # make LBC exception
        if len(self.fa_list) > 0 and headgroup.sp_exception:
            self.fa_list[0].set_type(LipidFaBondType.LCB_EXCEPTION)
            
        return headgroup
    
    
    
    def assemble_lipid(self, headgroup):
        lipid_level_class = None
        if self.level == LipidLevel.COMPLETE_STRUCTURE: lipid_level_class = LipidCompleteStructure
        if self.level == LipidLevel.FULL_STRUCTURE: lipid_level_class = LipidFullStructure
        elif self.level == LipidLevel.STRUCTURE_DEFINED: lipid_level_class = LipidStructureDefined
        elif self.level == LipidLevel.SN_POSITION: lipid_level_class = LipidSnPosition
        elif self.level == LipidLevel.MOLECULAR_SPECIES: lipid_level_class = LipidMolecularSpecies
        elif self.level == LipidLevel.SPECIES: lipid_level_class = LipidSpecies
        
        return lipid_level_class(headgroup, self.fa_list)
