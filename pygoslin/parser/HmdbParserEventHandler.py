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
from pygoslin.domain.Adduct import Adduct
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidSpeciesInfo import LipidSpeciesInfo
from pygoslin.domain.LipidSpecies import LipidSpecies
from pygoslin.domain.FattyAcid import FattyAcid
from pygoslin.domain.LipidMolecularSubspecies import LipidMolecularSubspecies
from pygoslin.domain.LipidStructuralSubspecies import LipidStructuralSubspecies
from pygoslin.domain.LipidIsomericSubspecies import LipidIsomericSubspecies
from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.HeadGroup import HeadGroup
from pygoslin.domain.FunctionalGroup import *
from pygoslin.domain.LipidClass import *

class HmdbParserEventHandler(BaseParserEventHandler):
    
    def __init__(self):
        super().__init__()
        
        self.registered_events["lipid_pre_event"] = self.reset_lipid
        self.registered_events["lipid_post_event"] = self.build_lipid
        self.registered_events["fa_hg_pre_event"] = self.set_head_group_name
        self.registered_events["gl_hg_pre_event"] = self.set_head_group_name
        self.registered_events["gl_mono_hg_pre_event"] = self.set_head_group_name
        self.registered_events["gl_molecular_hg_pre_event"] = self.set_head_group_name
        self.registered_events["mediator_pre_event"] = self.mediator_event
        
        self.registered_events["pl_hg_pre_event"] = self.set_head_group_name
        self.registered_events["pl_three_hg_pre_event"] = self.set_head_group_name
        self.registered_events["pl_four_hg_pre_event"] = self.set_head_group_name
        self.registered_events["sl_hg_pre_event"] = self.set_head_group_name
        self.registered_events["st_species_hg_pre_event"] = self.set_head_group_name
        self.registered_events["st_sub1_hg_pre_event"] = self.set_head_group_name
        self.registered_events["st_sub2_hg_pre_event"] = self.set_head_group_name
        self.registered_events["fa_species_pre_event"] = self.set_species_level
        self.registered_events["gl_molecular_pre_event"] = self.set_molecular_level
        self.registered_events["unsorted_fa_separator_pre_event"] = self.set_molecular_level
        self.registered_events["fa2_unsorted_pre_event"] = self.set_molecular_level
        self.registered_events["fa3_unsorted_pre_event"] = self.set_molecular_level
        self.registered_events["fa4_unsorted_pre_event"] = self.set_molecular_level
        
        self.registered_events["db_single_position_pre_event"] = self.set_isomeric_level
        self.registered_events["db_single_position_post_event"] = self.add_db_position
        self.registered_events["db_position_number_pre_event"] = self.add_db_position_number
        self.registered_events["cistrans_pre_event"] = self.add_cistrans
        
        self.registered_events["lcb_pre_event"] = self.new_lcb
        self.registered_events["lcb_post_event"] = self.clean_lcb
        self.registered_events["fa_pre_event"] = self.new_fa
        self.registered_events["fa_post_event"] = self.append_fa
        self.registered_events["fa_lcb_suffix_type_pre_event"] = self.add_one_hydroxyl
        
        self.registered_events["ether_pre_event"] = self.add_ether
        self.registered_events["hydroxyl_pre_event"] = self.add_hydroxyl
        self.registered_events["db_count_pre_event"] = self.add_double_bonds
        self.registered_events["carbon_pre_event"] = self.add_carbon
        
        self.registered_events["furan_fa_pre_event"] = self.furan_fa
        self.registered_events["interlink_fa_pre_event"] = self.interlink_fa
        self.registered_events["lipid_suffix_pre_event"] = self.lipid_suffix
        self.registered_events["methyl_pre_event"] = self.add_methyl
        
        
    def reset_lipid(self, node):
        self.level = LipidLevel.ISOMERIC_SUBSPECIES
        self.lipid = None
        self.head_group = ""
        self.lcb = None
        self.fa_list = []
        self.current_fa = None
        self.db_position = 0
        self.db_cistrans = ""
        self.use_head_group = False
        self.db_number = -1
        

    def add_db_position(self, node):
        if self.current_fa != None:
            if type(self.current_fa.double_bonds) == int:
                self.db_numbers = self.current_fa.double_bonds
                self.current_fa.double_bonds = {}
            self.current_fa.double_bonds[self.db_position] = self.db_cistrans
        

    def add_db_position_number(self, node):
        self.db_position = int(node.get_text())
        

    def add_cistrans(self, node):
        self.db_cistrans = node.get_text()
        

    def set_head_group_name(self, node):
        self.head_group = node.get_text()
        
        
    def set_isomeric_level(self, node):
        self.db_position = 0
        self.db_cistrans = ""
        
    
    def set_species_level(self, node):
        self.level = self.level if self.level.value < LipidLevel.SPECIES.value else LipidLevel.SPECIES
        
        
    def set_molecular_level(self, node):
        self.level = self.level if self.level.value < LipidLevel.MOLECULAR_SUBSPECIES.value else LipidLevel.MOLECULAR_SUBSPECIES
        
        
    def mediator_event(self, node):
        self.use_head_group = True
        self.head_group = node.get_text()
        
        
    def new_fa(self, node):
        self.db_numbers = -1
        self.current_fa = FattyAcid("FA%i" % (len(self.fa_list) + 1))
        
        
    def new_lcb(self, node):
        self.lcb = FattyAcid("LCB")
        self.lcb.lcb = True
        self.current_fa = self.lcb
        self.set_structural_subspecies_level(node)
            
            
    def clean_lcb(self, node):
        if type(self.current_fa.double_bonds) != int:
            if self.db_numbers > -1 and self.db_numbers != len(self.current_fa.double_bonds):
                raise LipidException("Double bond count does not match with number of double bond positions")
        elif self.current_fa.double_bonds > 0:
                self.set_structural_subspecies_level(node)
            
        self.current_fa = None
        
        
            
    def append_fa(self, node):
        if type(self.current_fa.double_bonds) != int:
            if self.db_numbers > -1 and self.db_numbers != len(self.current_fa.double_bonds):
                raise LipidException("Double bond count does not match with number of double bond positions")
        elif self.current_fa.double_bonds > 0:
                self.set_structural_subspecies_level(node)
            
        if self.level in {LipidLevel.STRUCTURAL_SUBSPECIES, LipidLevel.ISOMERIC_SUBSPECIES}:
            self.current_fa.position = len(self.fa_list) + 1
            
        self.fa_list.append(self.current_fa)
        self.current_fa = None
        
    
        
    def build_lipid(self, node):
        if self.lcb != None:
            for fa in self.fa_list: fa.position += 1
            self.fa_list = [self.lcb] + self.fa_list
        
        headgroup = HeadGroup(self.head_group, use_headgroup = self.use_head_group)
    
        max_num_fa = all_lipids[headgroup.lipid_class]["max_fa"]
        if max_num_fa != len(self.fa_list): self.set_structural_subspecies_level(node)

        lipid_level_class = None
        if self.level == LipidLevel.ISOMERIC_SUBSPECIES: lipid_level_class = LipidIsomericSubspecies
        if self.level == LipidLevel.STRUCTURAL_SUBSPECIES: lipid_level_class = LipidStructuralSubspecies
        if self.level == LipidLevel.MOLECULAR_SUBSPECIES: lipid_level_class = LipidMolecularSubspecies
        if self.level == LipidLevel.SPECIES: lipid_level_class = LipidSpecies
        
        self.lipid = LipidAdduct()
        self.lipid.lipid = lipid_level_class(headgroup, self.fa_list)
        
        self.content = self.lipid
        
        
    def add_ether(self, node):
        ether = node.get_text()
        if ether in {"o-", "O-"}: self.current_fa.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMANYL
        elif ether == "P-": self.current_fa.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMENYL
        else: raise UnsupportedLipidException("Fatty acyl chain of type '%s' is currently not supported" % ether)
    
    def add_methyl(self, node):
        functional_group = get_functional_group("Me").copy()
        functional_group.position = self.current_fa.num_carbon - (1 if node.get_text() == "i-" else 2)
        self.current_fa.num_carbon -= 1
        if "Me" not in self.current_fa.functional_groups: self.current_fa.functional_groups["Me"] = []
        self.current_fa.functional_groups["Me"].append(functional_group)
        
        
    def set_structural_subspecies_level(self, node):
        self.level = self.level if self.level.value < LipidLevel.STRUCTURAL_SUBSPECIES.value else LipidLevel.STRUCTURAL_SUBSPECIES
            
        
    def add_hydroxyl(self, node):
        hydroxyl, num_h = node.get_text(), 0
        if hydroxyl == "m": num_h = 1
        if hydroxyl == "d": num_h = 2
        if hydroxyl == "t": num_h = 3
        
        if get_category(self.head_group) == LipidCategory.SP and self.current_fa.lcb and self.head_group not in {"Cer", "LCB"}: num_h -= 1
        
        functional_group = get_functional_group("OH").copy()
        functional_group.count = num_h
        if "OH" not in self.current_fa.functional_groups: self.current_fa.functional_groups["OH"] = []
        self.current_fa.functional_groups["OH"].append(functional_group)
        
        
    def add_one_hydroxyl(self, node):
        functional_group = get_functional_group("OH").copy()
        if "OH" not in self.current_fa.functional_groups: self.current_fa.functional_groups["OH"] = []
        self.current_fa.functional_groups["OH"].append(functional_group)
        
        
    def add_double_bonds(self, node):
        self.current_fa.double_bonds = int(node.get_text())
        
        
    def add_carbon(self, node):
        self.current_fa.num_carbon += int(node.get_text())
        
        
    def furan_fa(self, node):
        raise UnsupportedLipidException("Furan fatty acyl chains are currently not supported")
        
        
    def interlink_fa(self, node):
        raise UnsupportedLipidException("Interconnected fatty acyl chains are currently not supported")


    def lipid_suffix(self, node):
        pass
        #raise UnsupportedLipidException("Lipids with suffix '%s' are currently not supported" % node.get_text())
