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
from pygoslin.domain.LipidCategory import LipidCategory
from pygoslin.domain.LipidClass import all_lipids
from pygoslin.domain.Adduct import Adduct
from pygoslin.domain.HeadGroup import HeadGroup
from pygoslin.domain.Element import element_positions

from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.FattyAcid import FattyAcid
from pygoslin.domain.FunctionalGroup import *
from pygoslin.domain.Cycle import *

from pygoslin.domain.LipidSpeciesInfo import LipidSpeciesInfo
from pygoslin.domain.LipidSpecies import LipidSpecies
from pygoslin.domain.LipidMolecularSubspecies import LipidMolecularSubspecies
from pygoslin.domain.LipidStructuralSubspecies import LipidStructuralSubspecies
from pygoslin.domain.LipidIsomericSubspecies import LipidIsomericSubspecies

from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.LipidClass import *

class ShorthandParserEventHandler(BaseParserEventHandler):
    
    def __init__(self):
        super().__init__()
        self.reset_lipid(None)
        
        self.registered_events["lipid_pre_event"] = self.reset_lipid
        self.registered_events["lipid_post_event"] = self.build_lipid
        
        ## set adduct events
        self.registered_events["adduct_info_pre_event"] = self.new_adduct
        self.registered_events["adduct_pre_event"] = self.add_adduct
        self.registered_events["charge_pre_event"] = self.add_charge
        self.registered_events["charge_sign_pre_event"] = self.add_charge_sign
        
        ## set species events
        self.registered_events["med_species_pre_event"] = self.set_species_level
        self.registered_events["gl_species_pre_event"] = self.set_species_level
        self.registered_events["gl_molecular_species_pre_event"] = self.set_molecular_level
        self.registered_events["pl_species_pre_event"] = self.set_species_level
        self.registered_events["pl_molecular_species_pre_event"] = self.set_molecular_level
        self.registered_events["sl_species_pre_event"] = self.set_species_level
        
        self.registered_events["pl_single_pre_event"] = self.set_molecular_level
        self.registered_events["unsorted_fa_separator_pre_event"] = self.set_molecular_level
        self.registered_events["ether_num_pre_event"] = self.set_ether_num
        
        ## set head groups events
        self.registered_events["med_hg_single_pre_event"] = self.set_headgroup_name
        self.registered_events["med_hg_double_pre_event"] = self.set_headgroup_name
        self.registered_events["med_hg_triple_pre_event"] = self.set_headgroup_name
        self.registered_events["gl_hg_single_pre_event"] = self.set_headgroup_name
        self.registered_events["gl_hg_double_pre_event"] = self.set_headgroup_name
        self.registered_events["gl_hg_true_double_pre_event"] = self.set_headgroup_name
        self.registered_events["gl_hg_triple_pre_event"] = self.set_headgroup_name
        self.registered_events["pl_hg_single_pre_event"] = self.set_headgroup_name
        self.registered_events["pl_hg_double_pre_event"] = self.set_headgroup_name
        self.registered_events["pl_hg_quadro_pre_event"] = self.set_headgroup_name
        self.registered_events["sl_hg_single_pre_event"] = self.set_headgroup_name
        self.registered_events["pl_hg_double_fa_hg_pre_event"] = self.set_headgroup_name
        self.registered_events["sl_hg_double_name_pre_event"] = self.set_headgroup_name
        self.registered_events["st_hg_pre_event"] = self.set_headgroup_name
        self.registered_events["st_hg_ester_pre_event"] = self.set_headgroup_name
        self.registered_events["hg_pip_pure_m_pre_event"] = self.set_headgroup_name
        self.registered_events["hg_pip_pure_d_pre_event"] = self.set_headgroup_name
        self.registered_events["hg_pip_pure_t_pre_event"] = self.set_headgroup_name
        self.registered_events["hg_PE_PS_pre_event"] = self.set_headgroup_name

        ## set head group headgroup_decorators
        self.registered_events["carbohydrate_pre_event"] = self.set_carbohydrate
        self.registered_events["carbohydrate_structural_pre_event"] = self.set_carbohydrate_structural
        self.registered_events["carbohydrate_isomeric_pre_event"] = self.set_carbohydrate_isomeric
        
        # fatty acyl events
        self.registered_events["lcb_post_event"] = self.set_lcb
        self.registered_events["fatty_acyl_chain_pre_event"] = self.new_fatty_acyl_chain
        self.registered_events["fatty_acyl_chain_post_event"] = self.add_fatty_acyl_chain
        self.registered_events["carbon_pre_event"] = self.set_carbon
        self.registered_events["db_count_pre_event"] = self.set_double_bond_count
        self.registered_events["db_position_number_pre_event"] = self.set_double_bond_position
        self.registered_events["db_single_position_pre_event"] = self.set_double_bond_information
        self.registered_events["db_single_position_post_event"] = self.add_double_bond_information
        self.registered_events["cistrans_pre_event"] = self.set_cistrans
        self.registered_events["ether_type_pre_event"] = self.set_ether_type
        
        ## set functional group events
        self.registered_events["func_group_data_pre_event"] = self.set_functional_group
        self.registered_events["func_group_data_post_event"] = self.add_functional_group
        self.registered_events["func_group_pos_number_pre_event"] = self.set_functional_group_position
        self.registered_events["func_group_name_pre_event"] = self.set_functional_group_name
        self.registered_events["func_group_count_pre_event"] = self.set_functional_group_count
        self.registered_events["stereo_type_pre_event"] = self.set_functional_group_stereo
        self.registered_events["molecular_func_group_name_pre_event"] = self.set_molecular_func_group
        
        ## set cycle events
        self.registered_events["func_group_cycle_pre_event"] = self.set_cycle
        self.registered_events["func_group_cycle_post_event"] = self.add_cycle
        self.registered_events["cycle_start_pre_event"] = self.set_cycle_start
        self.registered_events["cycle_end_pre_event"] = self.set_cycle_end
        self.registered_events["cycle_number_pre_event"] = self.set_cycle_number
        self.registered_events["cycle_db_cnt_pre_event"] = self.set_cycle_db_count
        self.registered_events["cycle_db_positions_pre_event"] = self.set_cycle_db_positions
        self.registered_events["cycle_db_positions_post_event"] = self.check_cycle_db_positions
        self.registered_events["cycle_db_position_number_pre_event"] = self.set_cycle_db_position
        self.registered_events["cycle_db_position_cis_trans_pre_event"] = self.set_cycle_db_position_cistrans
        self.registered_events["cylce_element_pre_event"] = self.add_cycle_element
        
        ## set linkage events
        self.registered_events["fatty_acyl_linkage_pre_event"] = self.set_acyl_linkage
        self.registered_events["fatty_acyl_linkage_post_event"] = self.add_acyl_linkage
        self.registered_events["fatty_alkyl_linkage_pre_event"] = self.set_alkyl_linkage
        self.registered_events["fatty_alkyl_linkage_post_event"] = self.add_alkyl_linkage
        self.registered_events["fatty_linkage_number_pre_event"] = self.set_fatty_linkage_number
        self.registered_events["fatty_acyl_linkage_sign_pre_event"] = self.set_linkage_type
        self.registered_events["hydrocarbon_chain_pre_event"] = self.set_hydrocarbon_chain
        self.registered_events["hydrocarbon_chain_post_event"] = self.add_hydrocarbon_chain
        self.registered_events["hydrocarbon_number_pre_event"] = self.set_fatty_linkage_number
        
        ## set remaining events
        self.registered_events["ring_stereo_pre_event"] = self.set_ring_stereo
        self.registered_events["pl_hg_fa_pre_event"] = self.set_hg_acyl
        self.registered_events["pl_hg_fa_post_event"] = self.add_hg_acyl
        self.registered_events["pl_hg_alk_pre_event"] = self.set_hg_alkyl
        self.registered_events["pl_hg_alk_post_event"] = self.add_hg_alkyl
        self.registered_events["pl_hg_species_pre_event"] = self.add_pl_species_data
        self.registered_events["hg_pip_m_pre_event"] = self.suffix_decorator_molecular
        self.registered_events["hg_pip_d_pre_event"] = self.suffix_decorator_molecular
        self.registered_events["hg_pip_t_pre_event"] = self.suffix_decorator_molecular
        self.registered_events["hg_PE_PS_type_pre_event"] = self.suffix_decorator_species
        



    def reset_lipid(self, node):
        self.level = LipidLevel.ISOMERIC_SUBSPECIES
        self.lipid = None
        self.headgroup = ""
        self.adduct = None
        self.fa_list = []
        self.current_fa = []
        self.headgroup_decorators = []
        self.tmp = {}
        #self.debug = "full"
        
        
    def set_lipid_level(self, level):
        self.level = self.level if self.level.value < level.value else level
        
        
        
    def add_cycle_element(self, node):
        element = node.get_text()
        if element not in element_positions:
            raise LipidParsingException("Element '%s' unknown" % element);
            
        self.tmp["fa%i" % len(self.current_fa)]["cycle_elements"].append(element_positions[element])
        
        
        
    def set_headgroup_name(self, node):
        if len(self.headgroup) == 0: self.headgroup = node.get_text()
        
        
    def set_carbohydrate(self, node):
        carbohydrate = node.get_text()
        try:
            functional_group = get_functional_group(carbohydrate).copy()
        except Exception:
            raise LipidParsingException("Carbohydrate '%s' unknown" % carbohydrate)
        
        functional_group.elements[Element.O] -= 1
        if "func_group_head" in self.tmp and self.tmp["func_group_head"]:
            self.headgroup_decorators.append(functional_group)
        else:
            if carbohydrate not in self.current_fa[-1].functional_groups:
                self.current_fa[-1].functional_groups[carbohydrate] = []
            self.current_fa[-1].functional_groups[carbohydrate].append(functional_group)
        
        
    def set_carbohydrate_structural(self, node):
        self.set_lipid_level(LipidLevel.STRUCTURAL_SUBSPECIES)
        self.tmp["func_group_head"] = True
        
        
    def set_carbohydrate_isomeric(self, node):
        self.tmp["func_group_head"] = True
        
        
    def suffix_decorator_molecular(self, node):
        self.headgroup_decorators.append(HeadgroupDecorator(node.get_text(), suffix = True, level = LipidLevel.MOLECULAR_SUBSPECIES))
        
        
    def suffix_decorator_species(self, node):
        self.headgroup_decorators.append(HeadgroupDecorator(node.get_text(), suffix = True, level = LipidLevel.SPECIES))
        
        
    
    def set_pl_hg_triple(self, node):
        self.set_molecular_level(node)
        self.set_headgroup_name(node)
        
        
    def set_ring_stereo(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_ring_stereo"] = node.get_text()
        
        
    def set_lcb(self, node):
        self.fa_list[-1].set_type(LipidFaBondType.LCB_REGULAR)
        
        
    def add_pl_species_data(self, node):
        self.set_lipid_level(LipidLevel.SPECIES)
        hgd = HeadgroupDecorator("")
        hgd.elements[Element.O] += 1
        hgd.elements[Element.H] -= 1
        self.headgroup_decorators.append(hgd)
    
        
    
    def new_fatty_acyl_chain(self, node):
        self.current_fa.append(FattyAcid("FA"))
        self.tmp["fa%i" % len(self.current_fa)] = {}
        
        
    def add_fatty_acyl_chain(self, node):
        fg_i = "fa%i" % (len(self.current_fa) - 2)
        special_type = ""
        if len(self.current_fa) >= 2 and fg_i in self.tmp and "fg_name" in self.tmp[fg_i]:
            if self.tmp[fg_i]["fg_name"] in {"acyl", "alkyl", "decorator_acyl", "decorator_alkyl", "cc"}:
                special_type = self.tmp[fg_i]["fg_name"]
                
        
        fa_i = "fa%i" % len(self.current_fa)
        if type(self.current_fa[-1].double_bonds) != int:
            if len(self.current_fa[-1].double_bonds) != self.tmp[fa_i]["db_count"]:
                raise LipidException("Double bond count does not match with number of double bond positions")
            
        else:
            if self.current_fa[-1].double_bonds > 0:
                self.set_lipid_level(LipidLevel.STRUCTURAL_SUBSPECIES)
                
        del self.tmp[fa_i]
        
        if len(special_type) > 0:
            fg_acyl_alkyl = self.current_fa.pop()
            fg_acyl_alkyl.name = special_type
            if special_type not in self.current_fa[-1].functional_groups: self.current_fa[-1].functional_groups[special_type] = []
            self.current_fa[-1].functional_groups[special_type].append(fg_acyl_alkyl)
            
        else:
            self.fa_list.append(self.current_fa.pop())
        
        
        
    def set_double_bond_position(self, node):
        fa_i = "fa%i" % len(self.current_fa)
        self.tmp[fa_i]["db_position"] = int(node.get_text())
        
        
        
    def set_double_bond_information(self, node):
        fa_i = "fa%i" % len(self.current_fa)
        self.tmp[fa_i]["db_position"] = 0
        self.tmp[fa_i]["db_cistrans"] = ""
        
        
        
    def add_double_bond_information(self, node):
        fa_i = "fa%i" % len(self.current_fa)
        pos = self.tmp[fa_i]["db_position"]
        cistrans = self.tmp[fa_i]["db_cistrans"]
        
        if cistrans == "": self.set_lipid_level(LipidLevel.STRUCTURAL_SUBSPECIES)
        
        del self.tmp[fa_i]["db_position"]
        del self.tmp[fa_i]["db_cistrans"]
        if type(self.current_fa[-1].double_bonds) == int: self.current_fa[-1].double_bonds = {}
        self.current_fa[-1].double_bonds[pos] = cistrans
        
        
        
    def set_cistrans(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["db_cistrans"] = node.get_text()
        
        
        
    def set_functional_group(self, node):
        fa_i = "fa%i" % len(self.current_fa)
        self.tmp[fa_i]["fg_pos"] = -1
        self.tmp[fa_i]["fg_name"] = "O"
        self.tmp[fa_i]["fg_cnt"] = 1
        self.tmp[fa_i]["fg_stereo"] = ""
        self.tmp[fa_i]["fg_ring_stereo"] = ""
        
        
        
    def set_cycle(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_name"] = "cy"
        self.current_fa.append(Cycle(0))
        fa_i = "fa%i" % len(self.current_fa)
        self.tmp[fa_i] = {"cycle_elements": []}
        
        
        
    def add_cycle(self, node):
        cycle_elements = self.tmp["fa%i" % len(self.current_fa)]["cycle_elements"]
        del self.tmp["fa%i" % len(self.current_fa)]
        cycle = self.current_fa.pop()
        cycle.bridge_chain = cycle_elements
        
        
        if type(cycle.start) == int and type(cycle.end) == int:
            if cycle.end - cycle.start + 1 + len(cycle.bridge_chain) < cycle.cycle:
                raise ConstraintViolationException("Cycle length '%i' does not match with cycle description" % cycle.cycle)
        
        if "cy" not in self.current_fa[-1].functional_groups: self.current_fa[-1].functional_groups["cy"] = []
        self.current_fa[-1].functional_groups["cy"].append(cycle)
        
        
        
    def set_fatty_linkage_number(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["linkage_pos"] = int(node.get_text())
        
        
        
    def set_hg_acyl(self, node):
        self.tmp["fa%i" % len(self.current_fa)] = {"fg_name": "decorator_acyl"}
        self.current_fa.append(HeadgroupDecorator("decorator_acyl", suffix = True))
        self.tmp["fa%i" % len(self.current_fa)] = {}
        
        
    
    def add_hg_acyl(self, node):
        del self.tmp["fa%i" % len(self.current_fa)]
        self.headgroup_decorators.append(self.current_fa.pop())
        del self.tmp["fa%i" % len(self.current_fa)]
        
        
        
        
    def set_hg_alkyl(self, node):
        self.tmp["fa%i" % len(self.current_fa)] = {"fg_name": "decorator_alkyl"}
        self.current_fa.append(HeadgroupDecorator("decorator_alkyl", suffix = True))
        self.tmp["fa%i" % len(self.current_fa)] = {}
        
    
    
    def add_hg_alkyl(self, node):
        del self.tmp["fa%i" % len(self.current_fa)]
        self.headgroup_decorators.append(self.current_fa.pop())
        del self.tmp["fa%i" % len(self.current_fa)]
        
        
        
    def set_linkage_type(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["linkage_type"] = node.get_text() == "N"
        
        
        
    def set_hydrocarbon_chain(self, node):
        fa_i = "fa%i" % len(self.current_fa)
        self.tmp[fa_i]["fg_name"] = "cc"
        self.current_fa.append(CarbonChain(None))
        self.tmp["fa%i" % len(self.current_fa)] = {"linkage_pos": -1}
        
        
        
    def add_hydrocarbon_chain(self, node):
        fa_i = "fa%i" % len(self.current_fa)
        linkage_pos = self.tmp[fa_i]["linkage_pos"]
        
        del self.tmp[fa_i]
        cc = self.current_fa.pop()
        
        cc.position = linkage_pos
        if linkage_pos == -1: self.set_lipid_level(LipidLevel.STRUCTURAL_SUBSPECIES)
        
        if "cc" not in self.current_fa[-1].functional_groups: self.current_fa[-1].functional_groups["cc"] = []
        self.current_fa[-1].functional_groups["cc"].append(cc)
        
        
        
    def set_acyl_linkage(self, node):
        fa_i = "fa%i" % len(self.current_fa)
        self.tmp[fa_i]["fg_name"] = "acyl"
        self.current_fa.append(AcylAlkylGroup(None))
        self.tmp["fa%i" % len(self.current_fa)] = {"linkage_pos": -1}
        

        
    def add_acyl_linkage(self, node):
        fa_i = "fa%i" % len(self.current_fa)
        linkage_type = self.tmp[fa_i]["linkage_type"]
        linkage_pos = self.tmp[fa_i]["linkage_pos"]
        
        del self.tmp[fa_i]
        acyl = self.current_fa.pop()
        
        acyl.position = linkage_pos
        acyl.set_N_bond_type(linkage_type)
        if linkage_pos == -1: self.set_lipid_level(LipidLevel.STRUCTURAL_SUBSPECIES)
        
        if "acyl" not in self.current_fa[-1].functional_groups: self.current_fa[-1].functional_groups["acyl"] = []
        self.current_fa[-1].functional_groups["acyl"].append(acyl)
        
        
        
    def set_alkyl_linkage(self, node):
        fa_i = "fa%i" % len(self.current_fa)
        self.tmp[fa_i]["fg_name"] = "alkyl"
        self.current_fa.append(AcylAlkylGroup(None, alkyl = True))
        self.tmp["fa%i" % len(self.current_fa)] = {"linkage_pos": -1}
        
        
        
    def add_alkyl_linkage(self, node):
        linkage_pos = self.tmp["fa%i" % len(self.current_fa)]["linkage_pos"]
        
        del self.tmp["fa%i" % len(self.current_fa)]
        alkyl = self.current_fa.pop()
        
        alkyl.position = linkage_pos
        if linkage_pos == -1: self.set_lipid_level(LipidLevel.STRUCTURAL_SUBSPECIES)
        
        if "alkyl" not in self.current_fa[-1].functional_groups: self.current_fa[-1].functional_groups["alkyl"] = []
        self.current_fa[-1].functional_groups["alkyl"].append(alkyl)
        
        
        
        
    def set_cycle_start(self, node):
        self.current_fa[-1].start = int(node.get_text())
        
        
    def set_cycle_end(self, node):
        self.current_fa[-1].end = int(node.get_text())
        
        
    def set_cycle_number(self, node):
        self.current_fa[-1].cycle = int(node.get_text())
        
        
    def set_cycle_db_count(self, node):
        self.current_fa[-1].double_bonds = int(node.get_text())
    
    
    def set_cycle_db_positions(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["cycle_db"] = self.current_fa[-1].double_bonds
        self.current_fa[-1].double_bonds = {}
    
    
    def check_cycle_db_positions(self, node):
        if len(self.current_fa[-1].double_bonds) != self.tmp["fa%i" % len(self.current_fa)]["cycle_db"]:
            raise LipidException("Double bond number in cycle does not correspond to number of double bond positions.")
    
    def set_cycle_db_position(self, node):
        pos = int(node.get_text())
        self.current_fa[-1].double_bonds[pos] = ""
        self.tmp["fa%i" % len(self.current_fa)]["last_db_pos"] = pos
        
        
    def set_cycle_db_position_cistrans(self, node):
        pos = self.tmp["fa%i" % len(self.current_fa)]["last_db_pos"]
        self.current_fa[-1].double_bonds[pos] = node.get_text()
        
    
    def set_functional_group_position(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_pos"] = int(node.get_text())
        
    
    def set_functional_group_name(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_name"] = node.get_text()
        
    
    def set_functional_group_count(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_cnt"] = int(node.get_text())
        
    
    def set_functional_group_stereo(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_stereo"] = node.get_text()
        
        
    def set_molecular_func_group(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_name"] = node.get_text()
        self.set_lipid_level(LipidLevel.MOLECULAR_SUBSPECIES)
        
        
    def add_functional_group(self, node):
        
        fa_i = "fa%i" % len(self.current_fa)
        fg_pos = self.tmp[fa_i]["fg_pos"]
        fg_name = self.tmp[fa_i]["fg_name"]
        fg_cnt = self.tmp[fa_i]["fg_cnt"]
        fg_stereo = self.tmp[fa_i]["fg_stereo"]
        fg_ring_stereo = self.tmp[fa_i]["fg_ring_stereo"]
        
        if fg_name in {"cy", "cc", "acyl", "alkyl", "decorator_acyl", "decorator_alkyl"}: return
        
        if fg_pos == -1:
            self.set_lipid_level(LipidLevel.STRUCTURAL_SUBSPECIES)
        
        try:
            functional_group = get_functional_group(fg_name).copy()
        except Exception:
            raise LipidParsingException(" '%s' unknown" % fg_name)
        
        functional_group.position = fg_pos
        functional_group.count = fg_cnt
        functional_group.stereochemistry = fg_stereo
        functional_group.ring_stereo = fg_ring_stereo
        
        del self.tmp[fa_i]["fg_pos"]
        del self.tmp[fa_i]["fg_name"]
        del self.tmp[fa_i]["fg_cnt"]
        del self.tmp[fa_i]["fg_stereo"]
        
        if fg_name not in self.current_fa[-1].functional_groups: self.current_fa[-1].functional_groups[fg_name] = []
        self.current_fa[-1].functional_groups[fg_name].append(functional_group)
        
        
    def set_ether_type(self, node):
        ether_type = node.get_text()
        if ether_type == "O-": self.current_fa[-1].lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMANYL
        elif ether_type == "P-": self.current_fa[-1].lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMENYL
        
        
    def set_ether_num(self, node):
        num_ethers, ether = 0, node.get_text()
        if ether == "d": num_ethers = 2
        elif ether == "t": num_ethers = 3
        elif ether == "e": num_ethers = 4
        self.tmp["num_ethers"] = num_ethers
        
        
    def set_species_level(self, node):
        self.set_lipid_level(LipidLevel.SPECIES)
        
        
        
    def set_molecular_level(self, node):
        self.set_lipid_level(LipidLevel.MOLECULAR_SUBSPECIES)
        
        
    
    def set_carbon(self, node):
        self.current_fa[-1].num_carbon = int(node.get_text())
      
      
    def set_double_bond_count(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["db_count"] = self.current_fa[-1].double_bonds = int(node.get_text())
        
        
        
    def build_lipid(self, node):
        
        headgroup = HeadGroup(self.headgroup, self.headgroup_decorators)

        true_fa = sum(1 for fa in self.fa_list if fa.num_carbon > 0 or (fa.double_bonds > 0 if type(fa.double_bonds) == int else len(fa.double_bonds)) > 0)
        
        poss_fa = all_lipids[headgroup.lipid_class]["poss_fa"]
        
        # make lyso
        can_be_lyso = "Lyso" in all_lipids[get_class("L" + self.headgroup)]["specials"] if get_class("L" + self.headgroup) < len(all_lipids) else False
        
        if true_fa + 1 == poss_fa and self.level != LipidLevel.SPECIES and headgroup.lipid_category == LipidCategory.GP and can_be_lyso:
            self.headgroup = "L" + self.headgroup
            headgroup = HeadGroup(self.headgroup, self.headgroup_decorators)
            poss_fa = all_lipids[headgroup.lipid_class]["poss_fa"]
        
        elif true_fa + 2 == poss_fa and self.level != LipidLevel.SPECIES and headgroup.lipid_category == LipidCategory.GP and self.head_group == "CL":
            self.head_group = "DL" + self.head_group
            headgroup = HeadGroup(self.head_group, self.headgroup_decorators)
            poss_fa = all_lipids[headgroup.lipid_class]["poss_fa"] if headgroup.lipid_class < len(all_lipids) else 0

        
        if self.level == LipidLevel.SPECIES:
            if true_fa == 0 and poss_fa != 0:
                raise ConstraintViolationException("No fatty acyl information lipid class '%s' provided." % headgroup.headgroup)
            
        elif true_fa != poss_fa and self.level in {LipidLevel.ISOMERIC_SUBSPECIES, LipidLevel.STRUCTURAL_SUBSPECIES}:
            raise ConstraintViolationException("Number of described fatty acyl chains (%i) not allowed for lipid class '%s' (having %i fatty aycl chains)." % (true_fa, headgroup.headgroup, poss_fa))
        
        
        if "HC" in all_lipids[headgroup.lipid_class]["specials"]:
            self.fa_list[0].lipid_FA_bond_type = LipidFaBondType.AMINE
        
        
        # make LBC exception
        if len(self.fa_list) > 0 and headgroup.sp_exception:
            self.fa_list[0].set_type(LipidFaBondType.LCB_EXCEPTION)
        
        
        #if len(self.fa_list) > 0 and not headgroup.sp_exception and self.level.value <= LipidLevel.MOLECULAR_SUBSPECIES.value and "O" in self.fa_list[0].functional_groups:
        #    self.fa_list[0].functional_groups["O"][0].count -= 1
        
        
        # add count numbers for fatty acyl chains
        fa_it = len(self.fa_list) > 0 and self.fa_list[0].lipid_FA_bond_type in {LipidFaBondType.LCB_EXCEPTION, LipidFaBondType.LCB_REGULAR}
        for it in range(fa_it, len(self.fa_list)):
            self.fa_list[it].name += "%i" % (it + 1)
            
        
        
        lipid_level_class = None
        if self.level == LipidLevel.ISOMERIC_SUBSPECIES: lipid_level_class = LipidIsomericSubspecies
        if self.level == LipidLevel.STRUCTURAL_SUBSPECIES: lipid_level_class = LipidStructuralSubspecies
        if self.level == LipidLevel.MOLECULAR_SUBSPECIES: lipid_level_class = LipidMolecularSubspecies
        if self.level == LipidLevel.SPECIES: lipid_level_class = LipidSpecies
        
        
        self.lipid = LipidAdduct()
        self.lipid.adduct = self.adduct
        self.lipid.lipid = lipid_level_class(headgroup, self.fa_list)
            
        if "num_ethers" in self.tmp: self.lipid.lipid.info.num_ethers = self.tmp["num_ethers"]
        l = self.lipid.lipid
        
        
        #if self.level == LipidLevel.SPECIES and l.headgroup.sp_exception and "O" in l.info.functional_groups:
        #    l.info.functional_groups["O"][0].count -= 1
        
        
        self.content = self.lipid
        
        
    def new_adduct(self, node):
        self.adduct = Adduct("", "", 0, 0)
        
        
    def add_adduct(self, node):
        self.adduct.adduct_string = node.get_text()
        
        
    def add_charge(self, node):
        self.adduct.charge = int (node.get_text())
        
        
    def add_charge_sign(self, node):
        sign = node.get_text()
        if sign == "+": self.adduct.set_charge_sign(1)
        if sign == "-": self.adduct.set_charge_sign(-1)
        
        
