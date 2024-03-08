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
from pygoslin.domain.HeadGroup import HeadGroup, glyco_table
from pygoslin.domain.FunctionalGroup import *

from pygoslin.domain.LipidCompleteStructure import LipidCompleteStructure
from pygoslin.domain.LipidFullStructure import LipidFullStructure
from pygoslin.domain.LipidStructureDefined import LipidStructureDefined
from pygoslin.domain.LipidSnPosition import LipidSnPosition
from pygoslin.domain.LipidMolecularSpecies import LipidMolecularSpecies
from pygoslin.domain.LipidSpecies import LipidSpecies

from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.Cycle import *


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
        
        
        
    def check_full_structure(obj):
        full = True
        if type(obj) == FattyAcid and obj.num_carbon == 0: return True
        if type(obj) == FattyAcid and type(obj.double_bonds) == int and obj.double_bonds > 0: return False
        if type(obj) == FattyAcid and type(obj.double_bonds) == dict:
            full &= len(obj.double_bonds) == sum(1 if (v in {"E", "Z"} or (v == '' and k == obj.num_carbon - 1)) else 0 for k, v in obj.double_bonds.items())
        for fg_name in obj.functional_groups:
            for fg in obj.functional_groups[fg_name]:
                if fg.name == "X": continue
                if fg.position < 0: return False
                if obj.functional_groups != None: full &= LipidBaseParserEventHandler.check_full_structure(fg)
        return full
    
    
    def set_lipid_level(self, level):
        self.level = self.level if self.level.value < level.value else level
        
        
        
    def sp_regular_lcb(self):
        return get_category(self.head_group) == LipidCategory.SP and self.current_fa.lipid_FA_bond_type in FattyAcid.LCB_STATES and not (self.head_group in self.SP_EXCEPTION_CLASSES and len(self.headgroup_decorators) == 0)
    
    
    
    def prepare_headgroup_and_checks(self, allow_class_shift = True):
        
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
        
        if allow_class_shift:
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


        # check if all functional groups have a position to be full structure
        if self.level in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE}:
            for fa in self.fa_list:
                if not LipidBaseParserEventHandler.check_full_structure(fa):
                    self.set_lipid_level(LipidLevel.STRUCTURE_DEFINED)
                    break
                    

        if self.level == LipidLevel.SPECIES:
            if true_fa == 0 and poss_fa != 0:
                raise ConstraintViolationException("No fatty acyl information lipid class '%s' provided." % headgroup.headgroup)
            
        elif true_fa > poss_fa and self.level in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE, LipidLevel.STRUCTURE_DEFINED, LipidLevel.SN_POSITION, LipidLevel.MOLECULAR_SPECIES}:
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
    
    
    
    
    
    def resolve_fa_synonym(self, mediator_name):
        
        if mediator_name == "Palmitic acid": # FA 16:0
            return FattyAcid("FA", 16)
            
        elif mediator_name in {"Linoleic acid", "LA"}:
            return FattyAcid("FA", 18, {9: "Z", 12: "Z"})
            
        elif mediator_name in {"Arachidonic acid", "AA", "ARA"}:
            return FattyAcid("FA", 20, {5: "Z", 8: "Z", 11: "Z", 14: "Z"})
            
        elif mediator_name == "ALA":
            return FattyAcid("FA", 18, {9: "Z", 12: "Z", 15: "Z"})
            
        elif  mediator_name == "EPA":
            return FattyAcid("FA", 20, {5: "Z", 8: "Z", 11: "Z", 14: "Z", 17: "Z"})
            
        elif  mediator_name == "DHA":
            return FattyAcid("FA", 22, {4: "Z", 7: "Z", 10: "Z", 13: "Z", 16: "Z", 19: "Z"})
            
        elif mediator_name == "LTB4":
            f1, f2 = get_functional_group("OH"), get_functional_group("OH")
            f1.position = 5
            f2.position = 12
            fg = {"OH": [f1, f2]}
            return FattyAcid("FA", 20, {6: "Z", 8: "E", 10: "E", 14: "Z"}, functional_groups = fg)
            
        elif mediator_name in {"Resolvin D3", "RvD3"}:
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("OH"), get_functional_group("OH")
            f1.position = 4
            f2.position = 11
            f3.position = 17
            fg = {"OH": [f1, f2, f3]}
            return FattyAcid("FA", 22, {5: "E", 7: "E", 9: "E", 13: "Z", 15: "E", 19: "Z"}, functional_groups = fg)
            
        elif mediator_name in {"Maresin 1", "Mar1"}:
            f1, f2 = get_functional_group("OH"), get_functional_group("OH")
            f1.position = 7
            f2.position = 14
            fg = {"OH": [f1, f2]}
            return FattyAcid("FA", 22, {4: "Z", 8: "E", 10: "E", 12: "Z", 16: "Z", 19: "Z"}, functional_groups = fg)
            
        elif mediator_name == "Resolvin D1":
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("OH"), get_functional_group("OH")
            f1.position = 7
            f2.position = 8
            f3.position = 17
            fg = {"OH": [f1, f2, f3]}
            return FattyAcid("FA", 22, {4: "Z", 9: "E", 11: "E", 13: "Z", 15: "E", 19: "Z"}, functional_groups = fg)
            
        elif mediator_name == "Resolvin D2":
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("OH"), get_functional_group("OH")
            f1.position = 4
            f2.position = 16
            f3.position = 17
            fg = {"OH": [f1, f2, f3]}
            return FattyAcid("FA", 22, {4: "Z", 8: "E", 10: "Z", 12: "E", 14: "E", 19: "Z"}, functional_groups = fg)
            
        elif mediator_name == "Resolvin D5":
            f1, f2 = get_functional_group("OH"), get_functional_group("OH")
            f1.position = 7
            f2.position = 17
            fg = {"OH": [f1, f2]}
            return FattyAcid("FA", 22, {4: "Z", 8: "E", 10: "Z", 13: "Z", 15: "E", 19: "Z"}, functional_groups = fg)
            
        elif mediator_name == "Resolvin E1":
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("OH"), get_functional_group("OH")
            f1.position = 5
            f2.position = 13
            f3.position = 18
            fg = {"OH": [f1, f2, f3]}
            return FattyAcid("FA", 20, {6: "Z", 8: "E", 10: "E", 14: "Z", 16: "E"}, functional_groups = fg)
            
        elif mediator_name == "Resolvin E2":
            f1, f2 = get_functional_group("OH"), get_functional_group("OH")
            f1.position = 5
            f2.position = 18
            fg = {"OH": [f1, f2]}
            return FattyAcid("FA", 20, {6: "E", 8: "Z", 11: "Z", 14: "Z", 16: "E"}, functional_groups = fg)
            
        elif mediator_name == "TXB1":
            f1, f2, f3, f4 = get_functional_group("OH"), get_functional_group("OH"), get_functional_group("OH"), get_functional_group("oxy")
            f1.position = 15
            f2.position = 9
            f3.position = 11
            f4.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"OH": [f2, f3], "oxy": [f4]})]}
            return FattyAcid("FA", 20, {13: "E"}, functional_groups = fg)
            
        elif mediator_name == "TXB2":
            f1, f2, f3, f4 = get_functional_group("OH"), get_functional_group("OH"), get_functional_group("OH"), get_functional_group("oxy")
            f1.position = 15
            f2.position = 9
            f3.position = 11
            f4.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"OH": [f2, f3], "oxy": [f4]})]}
            return FattyAcid("FA", 20, {5: "Z", 13: "E"}, functional_groups = fg)
            
        elif mediator_name == "TXB3":
            f1, f2, f3, f4 = get_functional_group("OH"), get_functional_group("OH"), get_functional_group("OH"), get_functional_group("oxy")
            f1.position = 15
            f2.position = 9
            f3.position = 11
            f4.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"OH": [f2, f3], "oxy": [f4]})]}
            return FattyAcid("FA", 20, {5: "Z", 13: "E", 17: "Z"}, functional_groups = fg)
            
        elif mediator_name in {"PGF2alpha", "PGF2-alpha", "PGF2 alpha"}:
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("OH"), get_functional_group("OH")
            f1.position = 15
            f2.position = 9
            f3.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"OH": [f2, f3]})]}
            return FattyAcid("FA", 20, {5: "Z", 13: "E"}, functional_groups = fg)
            
        elif mediator_name == "PGD2":
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("OH"), get_functional_group("oxo")
            f1.position = 15
            f2.position = 9
            f3.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"OH": [f2], "oxo": [f3]})]}
            return FattyAcid("FA", 20, {5: "Z", 13: "E"}, functional_groups = fg)
            
        elif mediator_name == "PGE2":
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("oxo"), get_functional_group("OH")
            f1.position = 15
            f2.position = 9
            f3.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"OH": [f3], "oxo": [f2]})]}
            return FattyAcid("FA", 20, {5: "Z", 13: "E"}, functional_groups = fg)
            
        elif mediator_name == "PGB2":
            f1, f2 = get_functional_group("OH"), get_functional_group("OH")
            f1.position = 15
            f2.position = 9
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 1, functional_groups = {"OH": [f2]})]}
            return FattyAcid("FA", 20, {5: "Z", 13: "E"}, functional_groups = fg)
            
        elif mediator_name == "15d-PGJ2":
            f1, f2 = get_functional_group("OH"), get_functional_group("oxo")
            f1.position = 15
            f2.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 1, functional_groups = {"oxo": [f2]})]}
            return FattyAcid("FA", 20, {5: "Z", 12: "E", 14: "E"}, functional_groups = fg)
        
        elif mediator_name == "PGJ2":
            f1, f2 = get_functional_group("OH"), get_functional_group("oxo")
            f1.position = 15
            f2.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 1, functional_groups = {"OH": [f2]})]}
            return FattyAcid("FA", 20, {5: "Z", 13: "E"}, functional_groups = fg)
        
        elif mediator_name in {"PGEM", "PGE-M"}:
            f1, f2, f3 = get_functional_group("oxo"), get_functional_group("oxo"), get_functional_group("OH")
            f1.position = 15
            f2.position = 9
            f3.position = 11
            fg = {"oxo": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"OH": [f3], "oxo": [f2]})]}
            return FattyAcid("FA", 20, {5: "Z"}, functional_groups = fg)
            
        elif mediator_name in {"PGF1alpha", "PGF1-alpha", "PGF1 alpha"}:
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("OH"), get_functional_group("OH")
            f1.position = 15
            f2.position = 9
            f3.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"OH": [f2, f3]})]}
            return FattyAcid("FA", 20, {13: "E"}, functional_groups = fg)
            
        elif mediator_name == "PDX":
            f1, f2 = get_functional_group("OH"), get_functional_group("OH")
            f1.position = 0
            f2.position = 17
            fg = {"OH": [f1, f2]}
            return FattyAcid("FA", 22, {4: "Z", 7: "Z", 11: "E", 13: "Z", 15: "E", 19: "Z"}, functional_groups = fg)
            
        elif mediator_name in {"Oleic acid", "OA"}:
            return FattyAcid("FA", 18, {9: "Z"})
        
        elif mediator_name == "DGLA":
            return FattyAcid("FA", 20, {8: "Z", 11: "Z", 14: "Z"})
        
        elif mediator_name in {"iPF2alpha-VI", "iPF2-alpha-VI", "iPF2 alpha-VI"}:
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("oxo"), get_functional_group("OH")
            f1.position = 5
            f2.position = 9
            f3.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"oxo": [f2], "OH": [f3]})]}
            return FattyAcid("FA", 20, {6: "E", 14: "Z"}, functional_groups = fg)
            
        elif mediator_name == "PGE1":
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("oxo"), get_functional_group("OH")
            f1.position = 15
            f2.position = 9
            f3.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"OH": [f3], "oxo": [f2]})]}
            return FattyAcid("FA", 20, {13: "E"}, functional_groups = fg)
            
        elif mediator_name == "PGE3":
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("oxo"), get_functional_group("OH")
            f1.position = 15
            f2.position = 9
            f3.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"OH": [f3], "oxo": [f2]})]}
            return FattyAcid("FA", 20, {5: "Z", 13: "E", 17: "Z"}, functional_groups = fg)
        

            
        
        return None
