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
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.LipidClass import *
from pygoslin.domain.FattyAcid import FattyAcid
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidMolecularSubspecies import LipidMolecularSubspecies
from pygoslin.domain.LipidStructuralSubspecies import LipidStructuralSubspecies
from pygoslin.domain.LipidIsomericSubspecies import LipidIsomericSubspecies
from pygoslin.domain.LipidAdduct import LipidAdduct
from pygoslin.domain.LipidSpeciesInfo import LipidSpeciesInfo
from pygoslin.domain.LipidSpecies import LipidSpecies
from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.HeadGroup import HeadGroup
from pygoslin.domain.FunctionalGroup import *


class LipidMapsParserEventHandler(BaseParserEventHandler):
    def __init__(self):
        super().__init__()
        self.reset_lipid(None)
        
        self.registered_events["lipid_pre_event"] = self.reset_lipid
        self.registered_events["lipid_post_event"] = self.build_lipid
        
        self.registered_events["mediator_pre_event"] = self.mediator_event
        self.registered_events["sphingoxine_pre_event"] = self.mediator_event
        self.registered_events["fa_no_hg_pre_event"] = self.pure_fa
        
        self.registered_events["sgl_species_pre_event"] = self.set_species_level
        self.registered_events["tgl_species_pre_event"] = self.set_species_level
        self.registered_events["dpl_species_pre_event"] = self.set_species_level
        self.registered_events["cl_species_pre_event"] = self.set_species_level
        self.registered_events["dsl_species_pre_event"] = self.set_species_level
        self.registered_events["species_fa_pre_event"] = self.set_species_level
        self.registered_events["fa2_unsorted_pre_event"] = self.set_molecular_subspecies_level
        self.registered_events["fa3_unsorted_pre_event"] = self.set_molecular_subspecies_level
        self.registered_events["fa4_unsorted_pre_event"] = self.set_molecular_subspecies_level
        
        self.registered_events["hg_sgl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_gl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_cl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_dpl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_lpl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_fourpl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_threepl_pre_event"] = self.set_head_group_name
        self.registered_events["sphingosine_name_pre_event"] = self.set_head_group_name
        self.registered_events["sphinganine_name_pre_event"] = self.set_head_group_name
        self.registered_events["hg_dsl_pre_event"] = self.set_head_group_name
        self.registered_events["hg_lsl_pre_event"] = self.set_head_group_name
        self.registered_events["ch_pre_event"] = self.set_head_group_name
        self.registered_events["pk_hg_pre_event"] = self.set_head_group_name
        self.registered_events["hg_fa_pre_event"] = self.set_head_group_name
        self.registered_events["hg_che_pre_event"] = self.set_head_group_name
        self.registered_events["special_cer_pre_event"] = self.set_head_group_name
        self.registered_events["special_cer_hg_pre_event"] = self.set_head_group_name
        self.registered_events["omega_linoleoyloxy_Cer_pre_event"] = self.set_omega_head_group_name
        
        self.registered_events["glyco_struct_pre_event"] = self.add_glyco
        
        self.registered_events["lcb_pre_event"] = self.new_lcb
        self.registered_events["lcb_post_event"] = self.clean_lcb
        self.registered_events["fa_pre_event"] = self.new_fa
        self.registered_events["fa_post_event"] = self.append_fa
        
        self.registered_events["db_single_position_pre_event"] = self.set_isomeric_level
        self.registered_events["db_single_position_post_event"] = self.add_db_position
        self.registered_events["db_position_number_pre_event"] = self.add_db_position_number
        self.registered_events["cistrans_pre_event"] = self.add_cistrans
        
        self.registered_events["ether_pre_event"] = self.add_ether
        self.registered_events["hydroxyl_pre_event"] = self.add_hydroxyl
        self.registered_events["hydroxyl_lcb_pre_event"] = self.add_hydroxyl_lcb
        self.registered_events["db_count_pre_event"] = self.add_double_bonds
        self.registered_events["carbon_pre_event"] = self.add_carbon
    
        self.registered_events["structural_mod_pre_event"] = self.set_structural_subspecies_level
        self.registered_events["single_mod_pre_event"] = self.set_mod
        self.registered_events["mod_text_pre_event"] = self.set_mod_text
        self.registered_events["mod_pos_pre_event"] = self.set_mod_pos
        self.registered_events["mod_num_pre_event"] = self.set_mod_num
        self.registered_events["single_mod_post_event"] = self.add_functional_group
        
        self.debug = ""
        
        
    def reset_lipid(self, node):
        self.level = LipidLevel.ISOMERIC_SUBSPECIES
        self.lipid = None
        self.head_group = ""
        self.lcb = None
        self.fa_list = []
        self.current_fa = None
        self.use_head_group = False
        self.omit_fa = False
        self.db_position = 0
        self.db_cistrans = ""
        self.db_numbers = -1
        self.mod_text = ""
        self.mod_pos = -1
        self.mod_num = 1
        self.headgroup_decorators = []
        self.add_omega_linoleoyloxy_Cer = False
        
        
    def set_molecular_subspecies_level(self, node):
        self.level = self.level if self.level.value < LipidLevel.MOLECULAR_SUBSPECIES.value else LipidLevel.MOLECULAR_SUBSPECIES
        
        
    def mediator_event(self, node):
        self.use_head_group = True
        self.head_group = node.get_text()
        
        
        
    def set_isomeric_level(self, node):
        self.db_position = 0
        self.db_cistrans = ""
        

    def add_db_position(self, node):
        if self.current_fa != None:
            if type(self.current_fa.double_bonds) == int:
                self.db_numbers = self.current_fa.double_bonds
                self.current_fa.double_bonds = {self.db_position: self.db_cistrans}
            else:
                self.current_fa.double_bonds[self.db_position] = self.db_cistrans
            
            if self.db_cistrans not in {'E', 'Z'}:
                self.set_structural_subspecies_level(node)
        

    def add_db_position_number(self, node):
        self.db_position = int(node.get_text())
        

    def add_cistrans(self, node):
        self.db_cistrans = node.get_text()
        
        
        
    def pure_fa(self, node):
        self.head_group = "FA"
        
        
        
    def set_head_group_name(self, node):
        self.head_group = node.get_text()
        
        
        
    def set_omega_head_group_name(self, node):
        self.set_head_group_name(node)
        self.add_omega_linoleoyloxy_Cer = True
        
        
        
    def set_species_level(self, node):
        self.level = self.level if self.level.value < LipidLevel.SPECIES.value else LipidLevel.SPECIES
        
        
    def set_structural_subspecies_level(self, node):
        self.level = self.level if self.level.value < LipidLevel.STRUCTURAL_SUBSPECIES.value else LipidLevel.STRUCTURAL_SUBSPECIES


    def set_mod(self, node):
        self.mod_text = ""
        self.mod_pos = -1
        self.mod_num = 1


    def set_mod_text(self, node):
        self.mod_text = node.get_text()
        
        
    def add_glyco(self, node):
        glyco = node.get_text()
        try:
            functional_group = get_functional_group(glyco).copy()
        except Exception:
            raise LipidParsingException("Carbohydrate '%s' unknown" % glyco)
        
        functional_group.elements[Element.O] -= 1
        print("adding", functional_group.elements[Element.O])
        self.headgroup_decorators.append(functional_group)
        


    def set_mod_pos(self, node):
        self.mod_pos = int(node.get_text())


    def set_mod_num(self, node):
        self.mod_num = int(node.get_text())
        
        
        
    def add_functional_group(self, node):
        functional_group = get_functional_group(self.mod_text).copy()
        functional_group.position = self.mod_pos
        functional_group.count = self.mod_num
        self.current_fa.add_functional_group(functional_group)
          
          
    def new_fa(self, node):
        self.db_numbers = -1
        self.current_fa = FattyAcid("FA%i" % (len(self.fa_list) + 1))

        
            
    def append_fa(self, node):
        if type(self.current_fa.double_bonds) != int:
            if self.db_numbers > -1 and self.db_numbers != len(self.current_fa.double_bonds):
                raise LipidException("Double bond count does not match with number of double bond positions")
        elif self.current_fa.double_bonds > 0:
                self.set_structural_subspecies_level(node)
            
        
        if self.level in {LipidLevel.STRUCTURAL_SUBSPECIES, LipidLevel.ISOMERIC_SUBSPECIES}:
            self.current_fa.position = len(self.fa_list) + 1
            
        if self.current_fa.num_carbon == 0: self.omit_fa = True
        self.fa_list.append(self.current_fa)
        
        self.current_fa = None
        
        
        
        
        
    def new_lcb(self, node):
        self.lcb = FattyAcid("LCB")
        self.current_fa = self.lcb
        self.set_structural_subspecies_level(node)
        self.current_fa.lcb = True
            
            
    def clean_lcb(self, node):
        if type(self.current_fa.double_bonds) != int:
            if self.db_numbers > -1 and self.db_numbers != len(self.current_fa.double_bonds):
                raise LipidException("Double bond count does not match with number of double bond positions")
        elif self.current_fa.double_bonds > 0:
                self.set_structural_subspecies_level(node)
        self.current_fa = None
        
        
    def add_ether(self, node):
        ether = node.get_text()
        
        if ether == "O-": self.current_fa.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMANYL
        elif ether == "P-": self.current_fa.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMENYL
        
        
        
    def add_hydroxyl(self, node):
        num_h = int(node.get_text())
        functional_group = get_functional_group("OH").copy()
        
        if get_category(self.head_group) == LipidCategory.SP and self.current_fa.lcb and (self.head_group not in {"Cer", "LCB"} or len(self.headgroup_decorators) > 0): num_h -= 1
        functional_group.count = num_h
        if "OH" not in self.current_fa.functional_groups: self.current_fa.functional_groups["OH"] = []
        self.current_fa.functional_groups["OH"].append(functional_group)
        
        
        
    def add_hydroxyl_lcb(self, node):
        hydroxyl = node.get_text()
        num_h = 0
        if hydroxyl == "m": num_h = 1
        elif hydroxyl == "d": num_h = 2
        elif hydroxyl == "t": num_h = 3
        
        if get_category(self.head_group) == LipidCategory.SP and self.current_fa.lcb and (self.head_group not in {"Cer", "LCB"} or len(self.headgroup_decorators) > 0): num_h -= 1
        functional_group = get_functional_group("OH").copy()
        functional_group.count = num_h
        if "OH" not in self.current_fa.functional_groups: self.current_fa.functional_groups["OH"] = []
        self.current_fa.functional_groups["OH"].append(functional_group)
        
        
    def add_double_bonds(self, node):
        self.current_fa.double_bonds += int(node.get_text())
        
        
    def add_carbon(self, node):
        self.current_fa.num_carbon = int(node.get_text())
        
        

    def build_lipid(self, node):
        if self.omit_fa and self.head_group in set(["PA", "PC", "PE", "PG", "PI", "PS"]):
            self.head_group = "L" + self.head_group
        
        
        if self.lcb != None:
            for fa in self.fa_list: fa.position += 1
            self.fa_list = [self.lcb] + self.fa_list
            
            
        if self.add_omega_linoleoyloxy_Cer:
            if len(self.fa_list) != 2:
                raise LipidException("omega-linoleoyloxy-Cer with a different combination to one long chain base and one fatty acyl chain unknown")
            if "acyl" not in self.fa_list[-1].functional_groups: self.fa_list[-1].functional_groups["acyl"] = []
            self.fa_list[-1].functional_groups["acyl"].append(AcylAlkylGroup(FattyAcid("FA", 18, {9: "Z", 12: "Z"})))
        
        headgroup = HeadGroup(self.head_group, self.headgroup_decorators, self.use_head_group)
    
        max_num_fa = all_lipids[headgroup.lipid_class]["max_fa"]
        if max_num_fa != len(self.fa_list): self.level = self.level if self.level.value < LipidLevel.MOLECULAR_SUBSPECIES.value else LipidLevel.MOLECULAR_SUBSPECIES

        lipid_level_class = None
        if self.level == LipidLevel.ISOMERIC_SUBSPECIES: lipid_level_class = LipidIsomericSubspecies
        if self.level == LipidLevel.STRUCTURAL_SUBSPECIES: lipid_level_class = LipidStructuralSubspecies
        if self.level == LipidLevel.MOLECULAR_SUBSPECIES: lipid_level_class = LipidMolecularSubspecies
        if self.level == LipidLevel.SPECIES: lipid_level_class = LipidSpecies
        
        self.lipid = LipidAdduct()
        self.lipid.lipid = lipid_level_class(headgroup, self.fa_list)
        
        self.content = self.lipid
        
        
