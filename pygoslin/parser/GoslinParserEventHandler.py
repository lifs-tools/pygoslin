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


from pygoslin.parser.LipidBaseParserEventHandler import LipidBaseParserEventHandler
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
from pygoslin.domain.LipidSpeciesInfo import LipidSpeciesInfo

from pygoslin.domain.LipidExceptions import *

class GoslinParserEventHandler(LipidBaseParserEventHandler):
    
    def __init__(self):
        super().__init__()
        self.reset_lipid(None)
        
        self.registered_events["lipid_pre_event"] = self.reset_lipid
        self.registered_events["lipid_post_event"] = self.build_lipid
        
        self.registered_events["hg_cl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_mlcl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_pl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_lpl_pre_event"] = self.set_lpl_head_group_name
        self.registered_events["hg_lsl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_dsl_pre_event"] = self.set_head_group_name
        self.registered_events["mediator_pre_event"] = self.set_head_group_name
        self.registered_events["hg_mgl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_dgl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_sgl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_tgl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_dlcl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_sac_di_pre_event"] = self.set_head_group_name
        self.registered_events["hg_sac_f_pre_event"] = self.set_head_group_name
        self.registered_events["hg_tpl_pre_event"] = self.set_head_group_name
        self.registered_events["st_pre_event"] = self.set_head_group_name
        self.registered_events["hg_ste_pre_event"] = self.set_head_group_name
        self.registered_events["hg_stes_pre_event"] = self.set_head_group_name
        
        self.registered_events["gl_species_pre_event"] = self.set_species_level
        self.registered_events["pl_species_pre_event"] = self.set_species_level
        self.registered_events["sl_species_pre_event"] = self.set_species_level
        self.registered_events["fa2_unsorted_pre_event"] = self.set_molecular_subspecies_level
        self.registered_events["fa3_unsorted_pre_event"] = self.set_molecular_subspecies_level
        self.registered_events["fa4_unsorted_pre_event"] = self.set_molecular_subspecies_level
        self.registered_events["slbpa_pre_event"] = self.set_molecular_subspecies_level
        self.registered_events["dlcl_pre_event"] = self.set_molecular_subspecies_level
        self.registered_events["mlcl_pre_event"] = self.set_molecular_subspecies_level
        
        self.registered_events["lcb_pre_event"] = self.new_lcb
        self.registered_events["lcb_post_event"] = self.clean_lcb
        self.registered_events["fa_pre_event"] = self.new_fa
        self.registered_events["fa_post_event"] = self.append_fa
        
        self.registered_events["db_single_position_pre_event"] = self.set_isomeric_level
        self.registered_events["db_single_position_post_event"] = self.add_db_position
        self.registered_events["db_position_number_pre_event"] = self.add_db_position_number
        self.registered_events["cistrans_pre_event"] = self.add_cistrans
        
        self.registered_events["ether_pre_event"] = self.add_ether
        self.registered_events["old_hydroxyl_pre_event"] = self.add_old_hydroxyl
        self.registered_events["db_count_pre_event"] = self.add_double_bonds
        self.registered_events["carbon_pre_event"] = self.add_carbon
        self.registered_events["hydroxyl_pre_event"] = self.add_hydroxyl
        self.registered_events["tpl_post_event"] = self.set_nape
        
        
        self.registered_events["adduct_info_pre_event"] = self.new_adduct
        self.registered_events["adduct_pre_event"] = self.add_adduct
        self.registered_events["charge_pre_event"] = self.add_charge
        self.registered_events["charge_sign_pre_event"] = self.add_charge_sign
        self.registered_events["plasmalogen_pre_event"] = self.add_plasmalogen
    
        self.debug = ""



    def reset_lipid(self, node):
        self.level = LipidLevel.FULL_STRUCTURE
        self.head_group = ""
        self.lcb = None
        self.fa_list = []
        self.current_fa = None
        self.adduct = None
        self.db_position = 0
        self.db_cistrans = ""
        self.unspecified_ether = False
        self.db_numbers = -1
        self.headgroup_decorators = []
        self.plasmalogen = ""
        
        
        
    def add_plasmalogen(self, node):
        plasmalogen = node.get_text().upper()
        if plasmalogen in {"O", "P"}:
            self.plasmalogen = plasmalogen
        
        
    def set_lpl_head_group_name(self, node):
        self.set_lipid_level(LipidLevel.MOLECULAR_SPECIES)
        self.head_group = node.get_text()
        
        

    def set_head_group_name(self, node):
        self.head_group = node.get_text()
        
        
        
    def set_unspecified_ether(self, node):
        self.unspecified_ether = True
        

        
    def set_nape(self, node):
        if self.head_group == "NAPE":
            self.head_group = "PE-N"
            hgd = HeadgroupDecorator("decorator_acyl", suffix = True)
            self.headgroup_decorators.append(hgd)
            hgd.functional_groups["decorator_acyl"] = [self.fa_list[-1]]
            self.fa_list.pop()
        
    
    
    def set_species_level(self, node):
        self.set_lipid_level(LipidLevel.SPECIES)
        
        
        
    def set_isomeric_level(self, node):
        self.db_position = 0
        self.db_cistrans = ""
        
        
        
    def set_structural_subspecies_level(self, node):
        self.set_lipid_level(LipidLevel.STRUCTURE_DEFINED)
        
        

    def add_db_position(self, node):
        if self.current_fa != None:
            if type(self.current_fa.double_bonds) == int:
                self.db_numbers = self.current_fa.double_bonds
                self.current_fa.double_bonds = {}
            self.current_fa.double_bonds[self.db_position] = self.db_cistrans
            if self.db_cistrans not in {"E", "Z"}: self.set_lipid_level(LipidLevel.STRUCTURE_DEFINED)
        


    def add_db_position_number(self, node):
        self.db_position = int(node.get_text())
        
        

    def add_cistrans(self, node):
        self.db_cistrans = node.get_text()
        
        
        
    def set_molecular_subspecies_level(self, node):
        self.set_lipid_level(LipidLevel.MOLECULAR_SPECIES)
        
        
        
    def new_fa(self, node):
        self.db_numbers = -1
        
        if self.unspecified_ether:
            lipid_FA_bond_type = LipidFaBondType.ETHER_UNSPECIFIED
            self.unspecified_ether = False
        else:
            lipid_FA_bond_type = LipidFaBondType.ESTER
        self.current_fa = FattyAcid("FA%i" % (len(self.fa_list) + 1), lipid_FA_bond_type = lipid_FA_bond_type)

        
            
    def append_fa(self, node):
        if type(self.current_fa.double_bonds) != int:
            if self.db_numbers > -1 and self.db_numbers != len(self.current_fa.double_bonds):
                raise LipidException("Double bond count does not match with number of double bond positions")
        elif self.current_fa.double_bonds > 0:
                self.set_lipid_level(LipidLevel.SN_POSITION)
            
        if self.current_fa.lipid_FA_bond_type == LipidFaBondType.ETHER_UNSPECIFIED:
            raise LipidException("Lipid with unspecified ether bond cannot be treated properly.")
        
        if self.level in {LipidLevel.SN_POSITION, LipidLevel.STRUCTURE_DEFINED, LipidLevel.FULL_STRUCTURE, LipidLevel.COMPLETE_STRUCTURE}:
            self.current_fa.position = len(self.fa_list) + 1
            
            
        self.fa_list.append(self.current_fa)
        self.current_fa = None
        
        
        
        
    def new_lcb(self, node):
        self.lcb = FattyAcid("LCB")
        self.current_fa = self.lcb
        self.set_structural_subspecies_level(node)
        self.lcb.set_type(LipidFaBondType.LCB_REGULAR)
        self.set_lipid_level(LipidLevel.STRUCTURE_DEFINED)
            
            
            
    def clean_lcb(self, node):
        if type(self.current_fa.double_bonds) != int:
            if self.db_numbers > -1 and self.db_numbers != len(self.current_fa.double_bonds):
                raise LipidException("Double bond count does not match with number of double bond positions")
                
        elif self.current_fa.double_bonds > 0:
                self.set_lipid_level(LipidLevel.SN_POSITION)
        self.current_fa = None
        
        
        
    def build_lipid(self, node):
        if self.lcb != None:
            for fa in self.fa_list: fa.position += 1
            self.fa_list = [self.lcb] + self.fa_list
            
        if self.plasmalogen != "" and self.lcb == None and len(self.fa_list) > 0:
            self.fa_list[0].lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMANYL if self.plasmalogen == "O" else LipidFaBondType.ETHER_PLASMENYL
        
        headgroup = self.prepare_headgroup_and_checks()

        lipid = LipidAdduct()
        lipid.adduct = self.adduct
        lipid.lipid = self.assemble_lipid(headgroup)
        
        self.content = lipid
        
        
        
    def add_ether(self, node):
        ether = node.get_text()
        if ether == "a": self.current_fa.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMANYL
        elif ether == "p":
            self.current_fa.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMENYL
            if type(self.current_fa.double_bonds) == int:
                self.current_fa.double_bonds = max(0, self.current_fa.double_bonds - 1)
        
        
    def add_old_hydroxyl(self, node):
        old_hydroxyl = node.get_text()
        
        num_h = 0
        if old_hydroxyl == "d": num_h = 2
        if old_hydroxyl == "t": num_h = 3
        
        if self.sp_regular_lcb(): num_h -= 1
        
        functional_group = get_functional_group("OH").copy()
        functional_group.count = num_h
        if "OH" not in self.current_fa.functional_groups: self.current_fa.functional_groups["OH"] = []
        self.current_fa.functional_groups["OH"].append(functional_group)
        
        
        
    def add_double_bonds(self, node):
        self.current_fa.double_bonds = int(node.get_text())
        
        
    def add_carbon(self, node):
        self.current_fa.num_carbon = int(node.get_text())
        
        
    def add_hydroxyl(self, node):
        num_h = int(node.get_text())
        functional_group = get_functional_group("OH").copy()
        self.set_lipid_level(LipidLevel.STRUCTURE_DEFINED)
        
        if self.sp_regular_lcb(): num_h -= 1
            
        functional_group.count = num_h
        if "OH" not in self.current_fa.functional_groups: self.current_fa.functional_groups["OH"] = []
        self.current_fa.functional_groups["OH"].append(functional_group)
        
        
    def new_adduct(self, node):
        self.adduct = Adduct("", "")
        
        
    def add_adduct(self, node):
        self.adduct.adduct_string = node.get_text()
        
        
    def add_charge(self, node):
        
        self.adduct.charge = int (node.get_text())
        
        
    def add_charge_sign(self, node):
        sign = node.get_text()
        if sign == "+": self.adduct.set_charge_sign(1)
        if sign == "-": self.adduct.set_charge_sign(-1)
        
        
        
