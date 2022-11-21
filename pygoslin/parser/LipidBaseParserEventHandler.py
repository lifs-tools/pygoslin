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
        # checking if head group is a glyco-sphingolipid
        hg = self.head_group.lower()
        if hg in glyco_table:
            for carbohydrate in glyco_table[hg]:
                try:
                    
                    functional_group = get_functional_group(carbohydrate)
                    functional_group.elements[Element.O] -= 1
                    self.headgroup_decorators.append(functional_group)
                except Exception:
                    raise LipidParsingException("Carbohydrate '%s' unknown" % carbohydrate)
            self.head_group = "Cer"
        
        
        headgroup = HeadGroup(self.head_group, self.headgroup_decorators, self.use_head_group)
        if self.use_head_group: return headgroup
    
        self.head_group = all_lipids[headgroup.lipid_class]["name"]
        max_num_fa = all_lipids[headgroup.lipid_class]["max_fa"]
        if max_num_fa != len(self.fa_list): self.set_lipid_level(LipidLevel.STRUCTURE_DEFINED)
    
        true_fa = sum(1 for fa in self.fa_list if fa.num_carbon > 0 or (fa.double_bonds > 0 if type(fa.double_bonds) == int else len(fa.double_bonds)) > 0)
        
        poss_fa = all_lipids[headgroup.lipid_class]["poss_fa"]
        max_fa = all_lipids[headgroup.lipid_class]["max_fa"]
        
        # make lyso
        can_be_lyso = "Lyso" in all_lipids[get_class("L" + self.head_group)]["specials"] if get_class("L" + self.head_group) < len(all_lipids) else False
        
        if (true_fa + 1 == poss_fa or true_fa + 2 == poss_fa) and self.level != LipidLevel.SPECIES and headgroup.lipid_category == LipidCategory.GP and can_be_lyso:
            if true_fa + 1 == poss_fa: self.head_group = "L" + self.head_group
            else: self.head_group = "DL" + self.head_group
            headgroup = HeadGroup(self.head_group, self.headgroup_decorators, self.use_head_group)
            poss_fa = all_lipids[headgroup.lipid_class]["poss_fa"] if headgroup.lipid_class < len(all_lipids) else 0
            
        elif (true_fa + 1 == poss_fa or true_fa + 2 == poss_fa) and self.level != LipidLevel.SPECIES and headgroup.lipid_category == LipidCategory.GL and self.head_group == "TG":
            if true_fa + 1 == poss_fa: self.head_group = "DG"
            else: self.head_group = "MG"
            headgroup = HeadGroup(self.head_group, self.headgroup_decorators, self.use_head_group)
            poss_fa = all_lipids[headgroup.lipid_class]["poss_fa"] if headgroup.lipid_class < len(all_lipids) else 0



        if self.level == LipidLevel.SPECIES:
            if true_fa == 0 and poss_fa != 0:
                raise ConstraintViolationException("No fatty acyl information lipid class '%s' provided." % headgroup.headgroup)
            
        elif true_fa != poss_fa and self.level in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE, LipidLevel.STRUCTURE_DEFINED, LipidLevel.SN_POSITION, LipidLevel.MOLECULAR_SPECIES}:
            raise ConstraintViolationException("Number of specified fatty acyl chains (%i) not allowed for lipid class '%s' (having %i fatty aycl chains)." % (true_fa, headgroup.headgroup, poss_fa))
        
        elif "Lyso" in all_lipids[get_class(self.head_group)]["specials"] and true_fa > poss_fa:
            raise ConstraintViolationException("Number of specified fatty acyl chains (%i) not allowed for lipid class '%s' (having %i fatty aycl chains)." % (true_fa, headgroup.headgroup, poss_fa))
            
        
        if "HC" in all_lipids[headgroup.lipid_class]["specials"] and len(self.fa_list) > 0:
            self.fa_list[0].lipid_FA_bond_type = LipidFaBondType.ETHER
            
        if "Amide" in all_lipids[headgroup.lipid_class]["specials"] and len(self.fa_list) > 0:
            for fatty in self.fa_list: fatty.lipid_FA_bond_type = LipidFaBondType.AMIDE
        
        # make LBC exception
        if len(self.fa_list) > 0 and headgroup.sp_exception:
            self.fa_list[0].set_type(LipidFaBondType.LCB_EXCEPTION)
            
        return headgroup
    
    
    
    def assemble_lipid(self, headgroup):
        
        for fa in self.fa_list:
            if fa.stereo_information_missing():
                self.set_lipid_level(LipidLevel.FULL_STRUCTURE)
                break
        
        lipid_level_class = None
        if self.level == LipidLevel.COMPLETE_STRUCTURE: lipid_level_class = LipidCompleteStructure
        elif self.level == LipidLevel.FULL_STRUCTURE: lipid_level_class = LipidFullStructure
        elif self.level == LipidLevel.STRUCTURE_DEFINED: lipid_level_class = LipidStructureDefined
        elif self.level == LipidLevel.SN_POSITION: lipid_level_class = LipidSnPosition
        elif self.level == LipidLevel.MOLECULAR_SPECIES: lipid_level_class = LipidMolecularSpecies
        elif self.level == LipidLevel.SPECIES: lipid_level_class = LipidSpecies
        
        return lipid_level_class(headgroup, self.fa_list)
