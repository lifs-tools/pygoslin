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


from pygoslin.parser.LipidBaseParserEventHandler import LipidBaseParserEventHandler
from pygoslin.domain.LipidAdduct import LipidAdduct
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.LipidClass import *
from pygoslin.domain.Adduct import Adduct

from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.FattyAcid import FattyAcid
from pygoslin.domain.HeadGroup import HeadGroup
from pygoslin.domain.FunctionalGroup import *
from pygoslin.domain.Cycle import Cycle

from pygoslin.domain.LipidCompleteStructure import LipidCompleteStructure
from pygoslin.domain.LipidFullStructure import LipidFullStructure
from pygoslin.domain.LipidStructureDefined import LipidStructureDefined
from pygoslin.domain.LipidSnPosition import LipidSnPosition
from pygoslin.domain.LipidMolecularSpecies import LipidMolecularSpecies
from pygoslin.domain.LipidSpecies import LipidSpecies
from pygoslin.domain.LipidSpeciesInfo import LipidSpeciesInfo

from pygoslin.domain.LipidExceptions import *

class GoslinParserEventHandler(LipidBaseParserEventHandler):
    
    mediator_FA = {'H': 17, 'O': 18, 'E': 20, 'Do': 22, 'D': 22}
    mediator_DB = {'M': 1, 'D': 2, 'Tr': 3, 'tr': 3, 'T': 4, 'P': 5, 'H': 6}
    
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
        self.registered_events["hg_so_lsl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_dsl_pre_event"] = self.set_head_group_name
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
        
        self.registered_events["mediator_pre_event"] = self.set_mediator
        self.registered_events["mediator_post_event"] = self.add_mediator
        self.registered_events["unstructured_mediator_pre_event"] = self.set_unstructured_mediator
        self.registered_events["trivial_mediator_pre_event"] = self.set_trivial_mediator
        self.registered_events["mediator_carbon_pre_event"] = self.set_mediator_carbon
        self.registered_events["mediator_db_pre_event"] = self.set_mediator_db
        self.registered_events["mediator_mono_functions_pre_event"] = self.set_mediator_function
        self.registered_events["mediator_di_functions_pre_event"] = self.set_mediator_function
        self.registered_events["mediator_tri_functions_pre_event"] = self.set_mediator_function
        self.registered_events["mediator_position_pre_event"] = self.set_mediator_function_position
        self.registered_events["mediator_functional_group_post_event"] = self.add_mediator_function
        self.registered_events["mediator_suffix_pre_event"] = self.add_mediator_suffix
        self.registered_events["mediator_tetranor_pre_event"] = self.set_mediator_tetranor
        
        
        self.registered_events["isotope_pair_pre_event"] = self.new_adduct
        self.registered_events["isotope_element_pre_event"] = self.set_heavy_d_element
        self.registered_events["isotope_number_pre_event"] = self.set_heavy_d_number
        self.registered_events["heavy_pre_event"] = self.new_adduct
        self.registered_events["adduct_heavy_element_pre_event"] = self.set_heavy_element
        self.registered_events["adduct_heavy_number_pre_event"] = self.set_heavy_number
        self.registered_events["adduct_heavy_component_post_event"] = self.add_heavy_component

        self.registered_events["prostaglandin_post_event"] = self.add_prostaglandin
        self.registered_events["prostaglandin_type_pre_event"] = self.set_prostaglandin_type
        self.registered_events["prostaglandin_number_pre_event"] = self.set_prostaglandin_number

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
        self.mediator_function = ""
        self.mediator_function_positions = []
        self.mediator_suffix = False
        self.use_head_group = False
        self.heavy_element = None
        self.heavy_element_number = 0
        self.trivial_mediator = False
        self.prostaglandin = {"type": "", "number": ""}
        
        
        
    def add_prostaglandin(self, node):
        if self.prostaglandin["type"] not in {"B", "D", "E", "F", "J", "K"} or self.prostaglandin["number"] not in {"1", "2", "3"}: return
    
        db = {}
        tmp_fa = self.current_fa
        if self.prostaglandin["number"] == "1": db = {13: "E"}
        elif self.prostaglandin["number"] == "2": db = {5: "Z", 13: "E"}
        elif self.prostaglandin["number"] == "3": db = {5: "Z", 13: "E", 17: "Z"}
            
        if self.prostaglandin["type"] == "B":
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("oxo")
            f1.position = 15
            f2.position = 9
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"oxo": [f2]})]}
            self.current_fa = FattyAcid("FA", 20, db, functional_groups = fg)
            
        elif self.prostaglandin["type"] == "D":
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("OH"), get_functional_group("oxo")
            f1.position = 15
            f2.position = 9
            f3.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"OH": [f2], "oxo": [f3]})]}
            self.current_fa = FattyAcid("FA", 20, db, functional_groups = fg)
        
        elif self.prostaglandin["type"] == "E":
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("oxo"), get_functional_group("OH")
            f1.position = 15
            f2.position = 9
            f3.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"OH": [f3], "oxo": [f2]})]}
            self.current_fa = FattyAcid("FA", 20, db, functional_groups = fg)
            
        elif self.prostaglandin["type"] == "F":
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("OH"), get_functional_group("OH")
            f1.position = 15
            f2.position = 9
            f3.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"OH": [f2, f3]})]}
            self.current_fa = FattyAcid("FA", 20, db, functional_groups = fg)
        
        
        elif self.prostaglandin["type"] == "J":
            f1, f2 = get_functional_group("OH"), get_functional_group("oxo")
            f1.position = 15
            f2.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 1, functional_groups = {"OH": [f2]})]}
            self.current_fa = FattyAcid("FA", 20, db, functional_groups = fg)
            
        elif self.prostaglandin["type"] == "K":
            f1, f2, f3 = get_functional_group("OH"), get_functional_group("oxo"), get_functional_group("oxo")
            f1.position = 15
            f2.position = 9
            f3.position = 11
            fg = {"OH": [f1], "cy": [Cycle(5, 8, 12, 0, functional_groups = {"oxo": [f2, f3]})]}
            self.current_fa = FattyAcid("FA", 20, db, functional_groups = fg)
            
        else:
            return
        
        self.clean_mediator(tmp_fa)
        self.fa_list = [self.current_fa]
        self.mediator_suffix = True
        self.prostaglandin = {"type": "", "number": ""}
        
        
        
    def set_prostaglandin_type(self, node):
        self.prostaglandin["type"] = node.get_text()
        
        
        
    def set_prostaglandin_number(self, node):
        self.prostaglandin["number"] = node.get_text()
        
        
        
    def set_mediator(self, node):
        self.head_group = "FA"
        self.current_fa = FattyAcid("FA")
        self.fa_list.append(self.current_fa)
        
        
        
    def set_unstructured_mediator(self, node):
        self.head_group = node.get_text()
        self.use_head_group = True
        self.fa_list = []
        
        
    def remove_deoxy(self, functional_groups):
        if "d" in functional_groups: del functional_groups["d"]
        for fg_name, fg_list in functional_groups.items():
            for fg in fg_list:
                self.remove_deoxy(fg.functional_groups)
        
        
    def clean_mediator(self, tmp_fa):
        if len(tmp_fa.functional_groups) > 0:
            for fg_name, fg_list in tmp_fa.functional_groups.items():
                for fg in fg_list:
                    if fg.position >= -1:
                        # erase prevously added functional group, if findable
                        def recursive_deletion(fg_name, fg, functional_groups):
                            for fg_name_curr, fg_list_curr in functional_groups.items():
                                del_fg = [i for i, fg_curr in enumerate(fg_list_curr) if fg_curr.position == fg.position]
                                if len(del_fg) > 0:
                                    for i in del_fg[::-1]: fg_list_curr.pop(i)
                                    if fg_name not in functional_groups: functional_groups[fg_name] = []
                                    functional_groups[fg_name].append(fg)
                                    return True
                                
                                for fg_curr in fg_list_curr:
                                    if recursive_deletion(fg_name, fg, fg_curr.functional_groups): return True
                            return False
                                
                        if not recursive_deletion(fg_name, fg, self.current_fa.functional_groups):
                            if fg_name not in self.current_fa.functional_groups: self.current_fa.functional_groups[fg_name] = []
                            self.current_fa.functional_groups[fg_name].append(fg)
                                    
                    else:
                        if fg_name not in self.current_fa.functional_groups: self.current_fa.functional_groups[fg_name] = []
                        self.current_fa.functional_groups[fg_name].append(fg)
        
        self.remove_deoxy(self.current_fa.functional_groups)
        
        
    def set_trivial_mediator(self, node):
        self.head_group = "FA"
        mediator_name = node.get_text()
            
        tmp_fa = self.current_fa
        self.current_fa = self.resolve_fa_synonym(mediator_name)
        self.clean_mediator(tmp_fa)
        self.fa_list = [self.current_fa]
        self.mediator_suffix = True
        
        
            
        
        
        
    def set_mediator_tetranor(self, node):
        self.current_fa.num_carbon -= 4
    
    
    
    def set_mediator_carbon(self, node):
        self.trivial_mediator = True
        self.current_fa.num_carbon += GoslinParserEventHandler.mediator_FA[node.get_text()]
        
        

    def set_mediator_db(self, node):
        self.current_fa.double_bonds = GoslinParserEventHandler.mediator_DB[node.get_text()]
        
        
        
    def set_mediator_function(self, node):
        self.mediator_function = node.get_text()
        
        
        
    def set_mediator_function_position(self, node):
        self.mediator_function_positions.append(int(node.get_text()))
        
        
        
    def add_mediator_function(self, node):
        functional_group, fg = None, ""
        if self.mediator_function in {"H", "hydroxy"}:
            functional_group, fg = get_functional_group("OH"), "OH"
            if len(self.mediator_function_positions) > 0: functional_group.position = self.mediator_function_positions[0]
            
        elif self.mediator_function.lower() in {"d", "deoxy"}:
            functional_group, fg = get_functional_group("d"), "d"
            if len(self.mediator_function_positions) == 0: functional = None
            else: functional_group.position = self.mediator_function_positions[0]
            
        elif self.mediator_function.lower() in {"oxo", "k", "keto"}:
            functional_group, fg = get_functional_group("oxo"), "oxo"
            if len(self.mediator_function_positions) > 0: functional_group.position = self.mediator_function_positions[0]
            
        elif self.mediator_function.lower() == "hp":
            functional_group, fg = get_functional_group("OOH"), "OOH"
            if len(self.mediator_function_positions) > 0: functional_group.position = self.mediator_function_positions[0]
            
        elif self.mediator_function == "NO2":
            functional_group, fg = get_functional_group("NO2"), "NO2"
            if len(self.mediator_function_positions) > 0: functional_group.position = self.mediator_function_positions[0]
            
        elif self.mediator_function in {"E", "Ep", "EP"}:
            functional_group, fg = get_functional_group("Ep"), "Ep"
            if len(self.mediator_function_positions) > 0: functional_group.position = self.mediator_function_positions[0]
            
        elif self.mediator_function.lower() in {"dh", "dih", "dihydro", 'dh'}:
            functional_group, fg = get_functional_group("OH"), "OH"
            if len(self.mediator_function_positions) >= 2:
                functional_group.position = self.mediator_function_positions[0]
                functional_group2 = get_functional_group("OH")
                functional_group2.position = self.mediator_function_positions[1]
                if "OH" not in self.current_fa.functional_groups: self.current_fa.functional_groups["OH"] = []
                self.current_fa.functional_groups["OH"].append(functional_group2)
            else:
                functional_group.count = 2
            
        elif self.mediator_function.lower() in {'trih', 'trihydroxy'}:
            functional_group, fg = get_functional_group("OH"), "OH"
            if len(self.mediator_function_positions) >= 3:
                functional_group.position = self.mediator_function_positions[0]
                functional_group2 = get_functional_group("OH")
                functional_group2.position = self.mediator_function_positions[1]
                functional_group3 = get_functional_group("OH")
                functional_group3.position = self.mediator_function_positions[2]
                if "OH" not in self.current_fa.functional_groups: self.current_fa.functional_groups["OH"] = []
                self.current_fa.functional_groups["OH"].append(functional_group2)
                self.current_fa.functional_groups["OH"].append(functional_group3)
            else:
                functional_group.count = 3
                
        if functional_group != None:
            if fg not in self.current_fa.functional_groups: self.current_fa.functional_groups[fg] = []
            self.current_fa.functional_groups[fg].append(functional_group)
        
        self.mediator_function_positions = []
        self.mediator_function = ""
        
    def add_mediator_suffix(self, node):
        self.mediator_suffix = True
        
        
        
    def set_heavy_d_element(self, node):
        self.adduct.heavy_elements[Element.H2] = 1
        
        
        
    def set_heavy_d_number(self, node):
        self.adduct.heavy_elements[Element.H2] = int(node.get_text())
        
        
        
    def add_mediator(self, node):
        if not self.mediator_suffix:
            if type(self.current_fa.double_bonds) == int: self.current_fa.double_bonds -= 1
        
        
        
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
        self.current_fa = FattyAcid("FA", lipid_FA_bond_type = lipid_FA_bond_type)

        
            
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
        
        self.remove_deoxy(self.current_fa.functional_groups)    
        self.fa_list.append(self.current_fa)
        
        oh_cnt = 1
        if self.head_group in {'Sa', 'So', 'S1P', 'Sa1P'}:
            #self.set_lipid_level(LipidLevel.STRUCTURE_DEFINED)
            functional_group = get_functional_group("OH").copy()
            if self.head_group in {'Sa', 'So'}:
                self.fa_list[0].lipid_FA_bond_type = LipidFaBondType.LCB_EXCEPTION
                functional_group.count = 2
            else:
                self.fa_list[0].lipid_FA_bond_type = LipidFaBondType.LCB_REGULAR
                functional_group.count = 1
            if "OH" not in self.current_fa.functional_groups: self.current_fa.functional_groups["OH"] = []
            self.current_fa.functional_groups["OH"].append(functional_group)
        
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
            self.fa_list = [self.lcb] + self.fa_list
            
        if self.plasmalogen != "" and self.lcb == None and len(self.fa_list) > 0:
            self.fa_list[0].lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMANYL if self.plasmalogen == "O" else LipidFaBondType.ETHER_PLASMENYL
        
        
        headgroup = self.prepare_headgroup_and_checks()
        
        if self.trivial_mediator and node.get_text() in trivial_mediators:
            if type(self.fa_list[0].double_bonds) != dict: self.fa_list[0].double_bonds = {}
            for p in trivial_mediators[node.get_text()]:
                if len(p) == 0: continue
                if p[-1] in {"E", "Z"}:
                    self.fa_list[0].double_bonds[int(p[:-1])] = p[-1]
                else:
                    self.fa_list[0].double_bonds[int(p)] = ""
            self.level = LipidLevel.FULL_STRUCTURE
            
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
        self.plasmalogen = ""
        
        
        
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
        if self.sp_regular_lcb(): num_h -= 1
        if num_h <= 0: return
        
        functional_group = get_functional_group("OH").copy()
        self.set_lipid_level(LipidLevel.STRUCTURE_DEFINED)
        functional_group.count = num_h
        if "OH" not in self.current_fa.functional_groups: self.current_fa.functional_groups["OH"] = []
        self.current_fa.functional_groups["OH"].append(functional_group)
        
        
        
    def new_adduct(self, node):
        if self.adduct == None: self.adduct = Adduct("", "")
        
        
        
    def add_adduct(self, node):
        self.adduct.adduct_string = node.get_text()
        
        
        
    def add_charge(self, node):
        self.adduct.charge = int (node.get_text())
        
        
        
    def add_charge_sign(self, node):
        sign = node.get_text()
        if sign == "+": self.adduct.set_charge_sign(1)
        if sign == "-": self.adduct.set_charge_sign(-1)
        if self.adduct.charge == 0: self.adduct.charge = 1
        
        
        
    def set_heavy_element(self, node):
        self.heavy_element = heavy_element_table[node.get_text()]
        self.heavy_element_number = 1
        
        
        
    def set_heavy_number(self, node):
        self.heavy_element_number = int(node.get_text())
        
        
        
    def add_heavy_component(self, node):
        self.adduct.heavy_elements[self.heavy_element] += self.heavy_element_number
        
        
trivial_mediators = {}
tm_dir_name = path.dirname(pygoslin.__file__)
tm_file_name = path.join(tm_dir_name, "data", "goslin", "trivial_mediators.csv")
with open(tm_file_name, mode = "rt", encoding= "utf-8") as tm_infile:
    for line in tm_infile:
        tokens = line.strip().split("\t")
        if len(tokens) != 2: continue
        trivial_mediators[tokens[0]] = [p for p in tokens[1].split(",")]
