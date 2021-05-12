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
from pygoslin.domain.Cycle import *

last_numbers = {'un': 1, 'do': 2, 'di': 2, 'tri': 3, 'buta': 4, 'but': 4, 'tetra': 4, 'penta': 5, 'pent': 5, 'hexa': 6, 'hex': 6, 'hepta': 7, 'hept': 7, 'octa': 8, 'oct': 8, 'nona': 9, 'non': 9}
second_numbers = {'deca': 10, 'dec': 10, 'cosa': 20, 'cos': 20, 'triaconta': 30, 'triacont': 30, 'tetraconta': 40, 'tetracont': 40, 'pentaconta': 50, 'pentacont': 50}
special_numbers = {'buta': 4, 'but': 4, 'eicosa': 20, 'eicos': 20, 'heneicosa': 21, 'heneicos': 21, 'prosta': 20, 'prost': 20, 'prostan': 20}
func_groups = {'keto': 'oxo', 'ethyl': 'Et', 'hydroxy': "OH", 'oxo': 'oxo', 'bromo': 'Br', 'methyl': 'Me', 'hydroperoxy': 'OOH', 'homo': '', 'Epoxy': 'Ep', 'fluoro': 'F', 'chloro': 'Cl', 'methylene': 'My', 'sulfooxy': 'S', 'amino': 'NH2'}

class FattyAcidParserEventHandler(BaseParserEventHandler):
    
    def __init__(self):
        super().__init__()
        
        self.registered_events["lipid_pre_event"] = self.reset_lipid
        self.registered_events["lipid_post_event"] = self.build_lipid
        self.registered_events["fatty_acid_post_event"] = self.set_fatty_acid
        
        self.registered_events["acid_single_type_pre_event"] = self.set_fatty_acyl_type
        self.registered_events["double_bond_position_pre_event"] = self.set_double_bond_information
        self.registered_events["double_bond_position_post_event"] = self.add_double_bond_information
        self.registered_events["db_number_post_event"] = self.set_double_bond_position
        self.registered_events["cistrans_post_event"] = self.set_cistrans
        
        ## lengths
        self.registered_events["functional_length_pre_event"] = self.reset_length
        self.registered_events["db_length_pre_event"] = self.reset_length
        self.registered_events["fatty_length_pre_event"] = self.reset_length
        self.registered_events["functional_length_post_event"] = self.set_functional_length
        self.registered_events["db_length_post_event"] = self.set_db_length
        self.registered_events["fatty_length_post_event"] = self.set_fatty_length
        
        ## numbers
        self.registered_events["notation_specials_pre_event"] = self.special_number
        self.registered_events["notation_last_digit_pre_event"] = self.last_number
        self.registered_events["notation_second_digit_pre_event"] = self.second_number
        
        ## functional groups
        self.registered_events["functional_group_pre_event"] = self.set_functional_group
        self.registered_events["functional_group_post_event"] = self.add_functional_group
        self.registered_events["functional_pos_pre_event"] = self.set_functional_pos
        self.registered_events["functional_group_type_pre_event"] = self.set_functional_type
        
        ## cyclo / epoxy
        self.registered_events["cyclo_pre_event"] = self.set_functional_group
        self.registered_events["cyclo_post_event"] = self.add_cyclo
        self.registered_events["epoxy_pre_event"] = self.set_functional_group
        self.registered_events["epoxy_post_event"] = self.add_epoxy
        
        ## dioic
        self.registered_events["dioic_pre_event"] = self.set_functional_group
        self.registered_events["dioic_post_event"] = self.set_dioic
        self.registered_events["dioic_acid_pre_event"] = self.set_fatty_acyl_type
        
        ## prosta
        self.registered_events["prosta_pre_event"] = self.set_prosta
        self.registered_events["prosta_post_event"] = self.add_cyclo
        self.registered_events["reduction_pre_event"] = self.set_functional_group
        self.registered_events["reduction_post_event"] = self.reduction
        
        ## recursion
        self.registered_events["recursion_description_pre_event"] = self.set_recursion
        self.registered_events["recursion_description_post_event"] = self.add_recursion
        self.registered_events["recursion_pos_pre_event"] = self.set_recursion_pos
        
        
        
        
    def reset_lipid(self, node):
        self.level = LipidLevel.ISOMERIC_SUBSPECIES
        self.lipid = None
        self.headgroup = ""
        self.current_fa = [FattyAcid("FA")]
        self.db_numbers = -1
        self.tmp = {"fa1": {}}
        #self.debug = "full"
        
        
        
    def set_double_bond_position(self, node):
        fa_i = "fa%i" % len(self.current_fa)
        pos = int(node.get_text())
        self.tmp[fa_i]["db_position"] = pos - (sum([1 for p in self.tmp["reduction"] if p < pos]) if "reduction" in self.tmp else 0)
        
        
    def set_recursion_pos(self, node):
        fa_i = "fa%i" % len(self.current_fa)
        self.tmp[fa_i]["recursion_pos"] = int(node.get_text())
        
        
    def set_recursion(self, node):
        self.current_fa.append(FattyAcid("FA"))
        fa_i = "fa%i" % len(self.current_fa)
        self.tmp[fa_i] = {"recursion_pos": 0}
        
        
    def add_recursion(self, node):
        fa_i = "fa%i" % len(self.current_fa)
        pos = self.tmp[fa_i]["recursion_pos"]
        
        fa = self.current_fa.pop()
        fa.position = pos
        
        if self.headgroup not in self.current_fa[-1].functional_groups: self.current_fa[-1].functional_groups[self.headgroup] = []
        self.current_fa[-1].functional_groups[self.headgroup].append(fa)
        
        
        
        
    def set_fatty_acid(self, node):
        if "noyloxy" in self.current_fa[-1].functional_groups and self.headgroup == "FA":
            self.headgroup = "FAHFA"
            
            acyl = AcylAlkylGroup(self.current_fa[-1].functional_groups["noyloxy"][0])
            acyl.position = fa.position
            del self.current_fa[-1].functional_groups["noyloxy"]
            
            if "acyl" not in self.current_fa[-1].functional_groups: self.current_fa[-1].functional_groups["acyl"] = []
            self.current_fa[-1].functional_groups["acyl"].append(acyl)
                
        elif "-1-yl" in self.current_fa[-1].functional_groups and self.headgroup == "cyclopentyl":
            self.current_fa[-1] += self.current_fa[-1].functional_groups["-1-yl"][0]
            del self.current_fa[-1].functional_groups["-1-yl"]
            
        elif "cyclopentyl" in self.current_fa[-1].functional_groups and self.headgroup == "FA":
            self.current_fa[-1] += self.current_fa[-1].functional_groups["cyclopentyl"][0]
            del self.current_fa[-1].functional_groups["cyclopentyl"]
        
        
        #else: 
        #    raise LipidException("Unknown nested lipid name structure")
        
        
    def set_double_bond_information(self, node):
        fa_i = "fa%i" % len(self.current_fa)
        self.tmp[fa_i]["db_position"] = 0
        self.tmp[fa_i]["db_cistrans"] = ""
        
        
    def reduction(self, node):
        self.current_fa[-1].num_carbon -= len(self.tmp["fg_pos"])
        self.tmp["reduction"] = [p for p in self.tmp["fg_pos"]]
        
        
    def add_double_bond_information(self, node):
        fa_i = "fa%i" % len(self.current_fa)
        pos = self.tmp[fa_i]["db_position"]
        cistrans = self.tmp[fa_i]["db_cistrans"].upper()
        if cistrans  in {'R', 'S', 'A', 'B'}: return
        
        del self.tmp[fa_i]["db_position"]
        del self.tmp[fa_i]["db_cistrans"]
        if type(self.current_fa[-1].double_bonds) == int: self.current_fa[-1].double_bonds = {}
        if pos not in self.current_fa[-1].double_bonds or len(self.current_fa[-1].double_bonds[pos]) == 0:
            self.current_fa[-1].double_bonds[pos] = cistrans
        
        
        
    def set_dioic(self, node):
        self.headgroup = "FA"
        
        pos = self.tmp["fg_pos"][1] if len(self.tmp["fg_pos"]) == 2 else self.current_fa[-1].num_carbon
        for fg in {"OH", "oxo"}:
            func_group = get_functional_group(fg)
            func_group.position = pos
            if fg not in self.current_fa[-1].functional_groups: self.current_fa[-1].functional_groups[fg] = []
            self.current_fa[-1].functional_groups[fg].append(func_group)
        
        
        
    def set_cistrans(self, node):
        self.tmp["fa%i" % len(self.current_fa)]["db_cistrans"] = node.get_text()
        
        
        
    def set_fatty_acyl_type(self, node):
        t = node.get_text()
        
        if t == "nol": self.headgroup = "FOH"
        elif t in {"noic acid", "dioic_acid"}: self.headgroup = "FA"
        elif t == "nal": self.headgroup = "FAL"
        elif t == "nyl acetate": self.headgroup = "WE"
        else: self.headgroup = t
        
        
        
    def set_lipid_level(self, level):
        self.level = self.level if self.level.value < level.value else level
        
        
    
    def reset_length(self, node):
        self.tmp["length"] = 0
    
    
    
    def set_functional_length(self, node):
        if self.tmp["length"] != len(self.tmp["fg_pos"]):
            raise LipidException("Length of functional group '%i' does not match with number of its positions '%i'" % (self.tmp["length"], len(self.tmp["fg_pos"])))
    
    
    
    def set_fatty_length(self, node):
        self.current_fa[-1].num_carbon += self.tmp["length"]
        
        
        
    def add_cyclo(self, node):
        
        start = self.tmp["fg_pos"][0]
        end = self.tmp["fg_pos"][1]
        
        cyclo_db = None
        # check double bonds
        if type(self.current_fa[-1].double_bonds) == dict and len(self.current_fa[-1].double_bonds) > 0:
            cyclo_db = {db_pos: val for db_pos, val in self.current_fa[-1].double_bonds.items() if start <= db_pos <= end}
            
            for pos in cyclo_db:
                del self.current_fa[-1].double_bonds[pos]
                
        # check functional_groups
        cyclo_fg = {}
        for fg, fg_list in self.current_fa[-1].functional_groups.items():
            remove_list = []
            for i, func_group in enumerate(fg_list):
                if start <= func_group.position <= end:
                    if fg not in cyclo_fg: cyclo_fg[fg] = []
                    cyclo_fg[fg].append(func_group)
                    remove_list.append(i)
                
            for i in remove_list[::-1]: del fg_list[i]
                
            
        cycle = Cycle(end - start + 1, start = start, end = end, double_bonds = cyclo_db, functional_groups = cyclo_fg)
        
        if "cy" not in self.current_fa[-1].functional_groups: self.current_fa[-1].functional_groups["cy"] = []
        self.current_fa[-1].functional_groups["cy"].append(cycle)
        
        
        
    def add_epoxy(self, node):
        self.tmp["fg_pos"] = self.tmp["fg_pos"][:1]
        self.tmp["fg_type"] = "Epoxy"
        
        
        
    def set_db_length(self, node):
        pass
        #d = self.current_fa[-1].get_double_bonds()
        
        #if len(self.current_fa[-1].double_bonds) != self.tmp["length"]:
        #    raise LipidException("Double bond count does not match with number of double bond positions")
        
    
    
    def special_number(self, node):
        self.tmp["length"] += special_numbers[node.get_text()]
        
        
        
    def last_number(self, node):
        self.tmp["length"] += last_numbers[node.get_text()]
        
        
        
    def second_number(self, node):
        self.tmp["length"] += second_numbers[node.get_text()]
        
        
        
    def set_functional_group(self, node):
        self.tmp["fg_pos"] = []
        self.tmp["fg_type"] = ""
        
        
        
    def set_prosta(self, node):
        minus_pos = (sum([1 for p in self.tmp["reduction"] if p < 8]) if "reduction" in self.tmp else 0)
        self.tmp["fg_pos"] = [8 - minus_pos, 12 - minus_pos]
        self.tmp["fg_type"] = "cy"
        
        
        
    def add_functional_group(self, node):
        t = self.tmp["fg_type"]
        if t not in func_groups: raise LipidException("Unknown functional group: '%s'" % t)
        t = func_groups[t]
        
        if len(t) == 0: return
    
        fg = get_functional_group(t)
        if t not in self.current_fa[-1].functional_groups: self.current_fa[-1].functional_groups[t] = []
        for pos in self.tmp["fg_pos"]:
            fg_insert = fg.copy()
            fg_insert.position = pos - (sum([1 for p in self.tmp["reduction"] if p < pos]) if "reduction" in self.tmp else 0)
            self.current_fa[-1].functional_groups[t].append(fg_insert)
        
        
    def set_functional_pos(self, node):
        self.tmp["fg_pos"].append(int(node.get_text()))
        
        
        
    def set_functional_type(self, node):
        self.tmp["fg_type"] = node.get_text()
        
        
        
    def build_lipid(self, node):
        if type(self.current_fa[-1].double_bonds) == dict and len(self.current_fa[-1].double_bonds) > 0:
            if sum(len(ct) > 0 for p, ct in self.current_fa[-1].double_bonds.items()) != len(self.current_fa[-1].double_bonds):
                self.set_lipid_level(LipidLevel.STRUCTURAL_SUBSPECIES)
        
        
        lipid_level_class = None
        if self.level == LipidLevel.ISOMERIC_SUBSPECIES: lipid_level_class = LipidIsomericSubspecies
        if self.level == LipidLevel.STRUCTURAL_SUBSPECIES: lipid_level_class = LipidStructuralSubspecies
        if self.level == LipidLevel.MOLECULAR_SUBSPECIES: lipid_level_class = LipidMolecularSubspecies
        if self.level == LipidLevel.SPECIES: lipid_level_class = LipidSpecies
        
        headgroup = HeadGroup(self.headgroup)
        
        self.lipid = LipidAdduct()
        self.lipid.lipid = lipid_level_class(headgroup, self.current_fa)
            
        self.content = self.lipid
