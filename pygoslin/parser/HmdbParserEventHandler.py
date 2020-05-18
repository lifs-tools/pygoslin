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
from pygoslin.domain.LipidExceptions import UnsupportedLipidException

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
        
        
    def reset_lipid(self, node):
        self.level = LipidLevel.STRUCTURAL_SUBSPECIES
        self.lipid = None
        self.head_group = ""
        self.lcb = None
        self.fa_list = []
        self.current_fa = None
        self.db_position = 0
        self.db_cistrans = ""
        self.use_head_group = False
        

    def add_db_position(self, node):
        if self.current_fa != None: self.current_fa.double_bond_positions[self.db_position] = self.db_cistrans
        

    def add_db_position_number(self, node):
        self.db_position = int(node.get_text())
        

    def add_cistrans(self, node):
        self.db_cistrans = node.get_text()
        

    def set_head_group_name(self, node):
        self.head_group = node.get_text()
        
        
    def set_isomeric_level(self, node):
        self.level = LipidLevel.ISOMERIC_SUBSPECIES
        self.db_position = 0
        self.db_cistrans = ""
        
    
    def set_species_level(self, node):
        self.level = LipidLevel.SPECIES
        
        
    def set_molecular_level(self, node):
        self.level = LipidLevel.MOLECULAR_SUBSPECIES
        
        
    def mediator_event(self, node):
        self.use_head_group = True
        self.head_group = node.get_text()
        
        
    def new_fa(self, node):
        self.current_fa = FattyAcid("FA%i" % (len(self.fa_list) + 1), 2, 0, 0, LipidFaBondType.ESTER, False, 0, {})
        
        
    def new_lcb(self, node):
        self.lcb = FattyAcid("LCB", 2, 0, 0, LipidFaBondType.ESTER, True, 1, {})
        self.current_fa = self.lcb
            
            
    def clean_lcb(self, node):
        lcb = self.lcb
        if self.level == LipidLevel.SPECIES:
            self.lcb = LipidSpeciesInfo(lcb)
            self.lcb.lipid_FA_bond_type = LipidFaBondType.ESTER
        
        self.current_fa = None
        
        
            
    def append_fa(self, node):
        current_fa = self.current_fa
        if self.level == LipidLevel.SPECIES:
            self.current_fa = LipidSpeciesInfo(current_fa)
            
        elif self.level in {LipidLevel.STRUCTURAL_SUBSPECIES, LipidLevel.ISOMERIC_SUBSPECIES}:
            self.current_fa.position = len(self.fa_list) + 1
            
            
            
            
        self.fa_list.append(self.current_fa)
        self.current_fa = None
        
        
        
        
        
    def build_lipid(self, node):
        if self.lcb != None:
            for fa in self.fa_list: fa.position += 1
            self.fa_list = [self.lcb] + self.fa_list
        
        lipid = None
        
        
        
        if self.level == LipidLevel.SPECIES:
            if len(self.fa_list) > 0:
                lipid_species_info = LipidSpeciesInfo(self.fa_list[0])
                lipid_species_info.level = LipidLevel.SPECIES
                lipid = LipidSpecies(self.head_group, lipid_species_info = lipid_species_info)
            else:
                lipid = LipidSpecies(self.head_group)
            
        elif self.level == LipidLevel.MOLECULAR_SUBSPECIES:
            lipid = LipidMolecularSubspecies(self.head_group, self.fa_list)
            
        elif self.level == LipidLevel.STRUCTURAL_SUBSPECIES:
            lipid = LipidStructuralSubspecies(self.head_group, self.fa_list)
            
        elif self.level == LipidLevel.ISOMERIC_SUBSPECIES:
            lipid = LipidIsomericSubspecies(self.head_group, self.fa_list)
            
        else:
            raise UnsupportedLipidException("Lipid level not supported")
    
        lipid.use_head_group = self.use_head_group
        self.lipid = LipidAdduct()
        self.lipid.lipid = lipid
        self.content = self.lipid
        
        
    def add_ether(self, node):
        ether = node.get_text()
        if ether in {"o-", "O-"}: self.current_fa.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMANYL
        elif ether == "P-": self.current_fa.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMENYL
        else: raise UnsupportedLipidException("Fatty acyl chain of type '%s' is currently not supported" % ether)
            
        
    def add_hydroxyl(self, node):
        hydroxyl = node.get_text()
        if hydroxyl == "m": self.current_fa.num_hydroxyl = 1
        if hydroxyl == "d": self.current_fa.num_hydroxyl = 2
        if hydroxyl == "t": self.current_fa.num_hydroxyl = 3
        
        
    def add_one_hydroxyl(self, node):
        self.current_fa.num_hydroxyl += 1
        
        
    def add_double_bonds(self, node):
        self.current_fa.num_double_bonds = int(node.get_text())
        
        
    def add_carbon(self, node):
        self.current_fa.num_carbon = int(node.get_text())
        
        
    def furan_fa(self, node):
        raise UnsupportedLipidException("Furan fatty acyl chains are currently not supported")
        
        
    def interlink_fa(self, node):
        raise UnsupportedLipidException("Interconnected fatty acyl chains are currently not supported")


    def lipid_suffix(self, node):
        raise UnsupportedLipidException("Lipids with suffix '%s' are currently not supported" % node.get_text())
