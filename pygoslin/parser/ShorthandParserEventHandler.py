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
from pygoslin.domain.FattyAcid import FattyAcid
from pygoslin.domain.FunctionalGroup import *
from pygoslin.domain.Cycle import *

from pygoslin.domain.LipidSpeciesInfo import LipidSpeciesInfo
from pygoslin.domain.LipidSpecies import LipidSpecies
from pygoslin.domain.LipidMolecularSubspecies import LipidMolecularSubspecies
from pygoslin.domain.LipidStructuralSubspecies import LipidStructuralSubspecies
from pygoslin.domain.LipidIsomericSubspecies import LipidIsomericSubspecies

from pygoslin.domain.LipidExceptions import *

class ShorthandParserEventHandler(BaseParserEventHandler):
    
    def __init__(self):
        super().__init__()
        self.reset_lipid(None)
        
        self.registered_events["lipid_pre_event"] = self.reset_lipid
        self.registered_events["lipid_post_event"] = self.build_lipid
        
        ## set categories
        self.registered_events["sl_pre_event"] = self.pre_sphingolipid
        self.registered_events["sl_post_event"] = self.post_sphingolipid
        self.registered_events["sl_hydroxyl_pre_event"] = self.set_hydroxyl
        
        ## set adduct events
        self.registered_events["adduct_info_pre_event"] = self.new_adduct
        self.registered_events["adduct_pre_event"] = self.add_adduct
        self.registered_events["charge_pre_event"] = self.add_charge
        self.registered_events["charge_sign_pre_event"] = self.add_charge_sign
        
        ## set species events
        self.registered_events["med_species_pre_event"] = self.set_species_level
        self.registered_events["gl_species_pre_event"] = self.set_species_level
        self.registered_events["gl_species_double_pre_event"] = self.set_species_level
        self.registered_events["pl_species_pre_event"] = self.set_species_level
        self.registered_events["pl_species_double_pre_event"] = self.set_species_level
        self.registered_events["pl_species_triple_pre_event"] = self.set_species_level
        self.registered_events["sl_species_pre_event"] = self.set_species_level
        
        ## set head groups events
        self.registered_events["med_hg_single_pre_event"] = self.set_head_group_name
        self.registered_events["med_hg_double_pre_event"] = self.set_head_group_name
        self.registered_events["med_hg_triple_pre_event"] = self.set_head_group_name
        self.registered_events["gl_hg_single_pre_event"] = self.set_head_group_name
        self.registered_events["gl_hg_double_pre_event"] = self.set_head_group_name
        self.registered_events["gl_hg_triple_pre_event"] = self.set_head_group_name
        self.registered_events["pl_hg_single_pre_event"] = self.set_head_group_name
        self.registered_events["pl_hg_double_pre_event"] = self.set_head_group_name
        self.registered_events["pl_hg_triple_pre_event"] = self.set_head_group_name
        self.registered_events["pl_hg_quadro_pre_event"] = self.set_head_group_name
        self.registered_events["sl_hg_single_pre_event"] = self.set_head_group_name
        self.registered_events["sl_hg_double_name_pre_event"] = self.set_head_group_name
        self.registered_events["st_hg_pre_event"] = self.set_head_group_name
        self.registered_events["st_hg_ester_pre_event"] = self.set_head_group_name

        ## set head group head_group_decorators
        self.registered_events["carbohydrate_pre_event"] = self.set_carbohydrate
        
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
        self.registered_events["func_group_placeholder_pre_event"] = self.set_placeholder
        self.registered_events["func_group_placeholder_number_pre_event"] = self.set_functional_group_count
        self.registered_events["func_group_cycle_pre_event"] = self.set_cycle
        self.registered_events["func_group_cycle_post_event"] = self.add_cycle
        self.registered_events["cycle_start_pre_event"] = self.set_cycle_start
        self.registered_events["cycle_end_pre_event"] = self.set_cycle_end
        self.registered_events["cycle_number_pre_event"] = self.set_cycle_number
        
        ## set linkage events
        
        self.registered_events["fatty_acyl_linkage_pre_event"] = self.set_acyl_linkage
        self.registered_events["fatty_acyl_linkage_post_event"] = self.add_acyl_linkage
        self.registered_events["fatty_alkyl_linkage_pre_event"] = self.set_alkyl_linkage
        self.registered_events["fatty_alkyl_linkage_post_event"] = self.add_alkyl_linkage
        self.registered_events["fatty_linkage_number_pre_event"] = self.set_fatty_linkage_number
        
        ## set remaining events
        self.registered_events["ring_stereo_pre_event"] = self.set_ring_stereo
        self.registered_events["pl_hg_fa_pre_event"] = self.set_hg_acyl
        self.registered_events["pl_hg_fa_post_event"] = self.add_hg_acyl
        self.registered_events["pl_hg_alk_pre_event"] = self.set_hg_alkyl
        self.registered_events["pl_hg_alk_post_event"] = self.add_hg_alkyl
        



    def reset_lipid(self, node):
        self.level = LipidLevel.ISOMERIC_SUBSPECIES
        self.lipid = None
        self.head_group = ""
        self.adduct = None
        self.fa_list = []
        self.current_fa = []
        self.head_group_decorators = []
        self.tmp = {}
        
        #self.debug = ""
        
        
    def set_lipid_level(self, level):
        self.level = self.level if self.level.value < level.value else level
        
        
    def set_head_group_name(self, node):
        if len(self.head_group) == 0: self.head_group = node.get_text()
        
        
    def set_carbohydrate(self, node):
        carbohydrate = node.get_text()
        try:
            functional_group = get_functional_group(carbohydrate).copy()
        except Exception:
            raise LipidParsingException("Carbohydrate '%s' unknown" % carbohydrate)
        
        self.head_group_decorators.append(functional_group)
        
        
    def pre_sphingolipid(self, node):
        self.tmp["sl_hydroxyl"] = False
        
        
    def set_ring_stereo(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_ring_stereo"] = node.get_text()
        
        
    def post_sphingolipid(self, node):
        if not self.tmp["sl_hydroxyl"] and self.head_group not in {"Cer", "SPB"}: self.set_lipid_level(LipidLevel.STRUCTURAL_SUBSPECIES)
        
        
    def set_hydroxyl(self, node):
        self.tmp["sl_hydroxyl"] = True
        
        
    def set_lcb(self, node):
        self.fa_list[-1].lcb = True
        self.fa_list[-1].name = "LCB"
        
    
    def new_fatty_acyl_chain(self, node):
        self.current_fa.append(FattyAcid("FA"))
        self.tmp["fa%i" % len(self.current_fa)] = {}
        
        
    def add_fatty_acyl_chain(self, node):
        
        fg_i = "fa%i" % (len(self.current_fa) - 2)
        special_type = ""
        if len(self.current_fa) > 2 and "fg_name" in self.tmp[fg_i]:
            if self.tmp[fg_i]["fg_name"] in {"acyl", "alkyl", "decorator_acyl", "decorator_alkyl"}:
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
        self.tmp[fa_i] = {}
        
        
    def add_cycle(self, node):
        del self.tmp["fa%i" % len(self.current_fa)]
        cycle = self.current_fa.pop()
        if "cy" not in self.current_fa[-1].functional_groups: self.current_fa[-1].functional_groups["cy"] = []
        self.current_fa[-1].functional_groups["cy"].append(cycle)
        
        
    def set_fatty_linkage_number(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["linkage_pos"] = int(node.get_text())
        
        
    def set_hg_acyl(self, node):
        self.tmp["fa%i" % len(self.current_fa)] = {"fg_name": "decorator_acyl"}
        self.current_fa.append(HeadGroupDecorator(None, suffix = True))
        self.tmp["fa%i" % len(self.current_fa)] = {}
        
    
    def add_hg_acyl(self, node):
        del self.tmp["fa%i" % len(self.current_fa)]
        
        self.head_group_decorators.append(self.current_fa.pop())
        del self.tmp["fa%i" % len(self.current_fa)]
        
        
        
        
    def set_hg_alkyl(self, node):
        self.tmp["fa%i" % len(self.current_fa)] = {"fg_name": "decorator_alkyl"}
        self.current_fa.append(HeadGroupDecorator(None, suffix = True))
        self.tmp["fa%i" % len(self.current_fa)] = {}
        
    
    def add_hg_alkyl(self, node):
        del self.tmp["fa%i" % len(self.current_fa)]
        
        self.head_group_decorators.append(self.current_fa.pop())
        del self.tmp["fa%i" % len(self.current_fa)]
        
        
        
    def set_acyl_linkage(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_name"] = "acyl"
        self.current_fa.append(AcylAlkylGroup(None))
        self.tmp["fa%i" % len(self.current_fa)] = {"linkage_pos": -1}

        
    def add_acyl_linkage(self, node):
        linkage_pos = self.tmp["fa%i" % len(self.current_fa)]["linkage_pos"]
        
        del self.tmp["fa%i" % len(self.current_fa)]
        acyl = self.current_fa.pop()
        
        acyl.position = linkage_pos
        if linkage_pos == -1: self.set_lipid_level(LipidLevel.STRUCTURAL_SUBSPECIES)
        
        if "acyl" not in self.current_fa[-1].functional_groups: self.current_fa[-1].functional_groups["acyl"] = []
        self.current_fa[-1].functional_groups["acyl"].append(acyl)
        
        
    def set_alkyl_linkage(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_name"] = "alkyl"
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
        
    
    def set_functional_group_position(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_pos"] = int(node.get_text())
        
    
    def set_functional_group_name(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_name"] = node.get_text()
        
    
    def set_functional_group_count(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_cnt"] = int(node.get_text())
        
    
    def set_functional_group_stereo(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_stereo"] = node.get_text()
        
        
    def set_placeholder(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["fg_name"] = "O"
        self.set_lipid_level(LipidLevel.MOLECULAR_SUBSPECIES)
        
        
    def add_functional_group(self, node):
        
        fa_i = "fa%i" % len(self.current_fa)
        fg_pos = self.tmp[fa_i]["fg_pos"]
        fg_name = self.tmp[fa_i]["fg_name"]
        fg_cnt = self.tmp[fa_i]["fg_cnt"]
        fg_stereo = self.tmp[fa_i]["fg_stereo"]
        fg_ring_stereo = self.tmp[fa_i]["fg_ring_stereo"]
        
        if fg_name in {"cy", "acyl", "alkyl", "decorator_acyl", "decorator_alkyl"}: return
        
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
        
        
    def set_species_level(self, node):
        self.level = self.level if self.level.value < LipidLevel.SPECIES.value else LipidLevel.SPECIES
        
    
    def set_carbon(self, node):
        self.current_fa[-1].num_carbon = int(node.get_text())
      
      
    def set_double_bond_count(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["db_count"] = self.current_fa[-1].double_bonds = int(node.get_text())
        
        
        
    def build_lipid(self, node):
        # add count numbers for fatty acyl chains
        fa_it = len(self.fa_list) > 0 and self.fa_list[0].lcb
        for it in range(fa_it, len(self.fa_list)):
            self.fa_list[it].name += "%i" % (it + 1)
        
        lipid_level_class = None
        if self.level == LipidLevel.ISOMERIC_SUBSPECIES: lipid_level_class = LipidIsomericSubspecies
        if self.level == LipidLevel.STRUCTURAL_SUBSPECIES: lipid_level_class = LipidStructuralSubspecies
        if self.level == LipidLevel.MOLECULAR_SUBSPECIES: lipid_level_class = LipidMolecularSubspecies
        if self.level == LipidLevel.SPECIES: lipid_level_class = LipidSpecies
        
        self.lipid = LipidAdduct()
        self.lipid.adduct = self.adduct
        self.lipid.lipid = lipid_level_class(self.head_group, self.fa_list)
        for decorator in self.head_group_decorators[::-1]:
            self.lipid.lipid.headgroup_decorators.append(decorator)
        
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
        
        
