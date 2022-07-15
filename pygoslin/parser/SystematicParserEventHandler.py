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

last_numbers = {'un': 1, 'hen': 1, 'do': 2, 'di': 2, 'tri': 3, 'buta': 4, 'but': 4, 'tetra': 4, 'penta': 5, 'pent': 5, 'hexa': 6, 'hex': 6, 'hepta': 7, 'hept': 7, 'octa': 8, 'oct': 8, 'nona': 9, 'non': 9}

second_numbers = {'deca': 10, 'dec': 10, 'eicosa': 20, 'eicos': 20 , 'cosa': 20, 'cos': 20, 'triaconta': 30, 'triacont': 30, 'tetraconta': 40, 'tetracont': 40, 'pentaconta': 50, 'pentacont': 50, 'hexaconta': 60, 'hexacont': 60, 'heptaconta': 70, 'heptacont': 70, 'octaconta': 80, 'octacont': 80, 'nonaconta': 90, 'nonacont': 90}

special_numbers = {'meth': 1, 'etha': 2, 'eth': 2,  'propa': 3, 'isoprop': 3, 'prop': 3, 'propi': 3, 'propio': 3, 'buta': 4, 'but': 4, 'butr': 4, 'valer': 5, 'eicosa': 20, 'eicos': 20, 'icosa': 20, 'icos': 20, 'prosta': 20, 'prost': 20, 'prostan': 20}

func_groups = {'keto': 'oxo', 'ethyl': 'Et', 'hydroxy': "OH", 'phospho': 'Ph', 'oxo': 'oxo', 'bromo': 'Br', 'methyl': 'Me', 'hydroperoxy': 'OOH', 'homo': '', 'Epoxy': 'Ep', 'fluro': 'F', 'fluoro': 'F', 'chloro': 'Cl', 'methylene': 'My', 'sulfooxy': 'Su', 'amino': 'NH2', 'sulfanyl': 'SH', 'methoxy': 'OMe', 'iodo': 'I', 'cyano': 'CN', 'nitro': 'NO2', 'OH': 'OH', 'thio': 'SH', 'mercapto': 'SH', 'carboxy': "COOH", 'acetoxy': 'Ac', 'cysteinyl': 'Cys', 'phenyl': 'Phe', 's-glutathionyl': "SGlu", 's-cysteinyl': "SCys", "butylperoxy": "BOO", "dimethylarsinoyl": "MMAs", "methylsulfanyl": "SMe", "imino": "NH", 's-cysteinylglycinyl': "SCG"}

ate = {'formate': 1, 'acetate': 2, 'butyrate': 4, 'propionate': 3, 'valerate': 5, 'isobutyrate': 4}


class SystematicParserEventHandler(BaseParserEventHandler):
    
    def __init__(self):
        super().__init__()
        
        self.registered_events["lipid_pre_event"] = self.reset_lipid
        self.registered_events["lipid_post_event"] = self.build_lipid
        
        ## FA category
        self.registered_events["FA_post_event"] = self.build_FA
        
        
        ## GL category
        self.registered_events["GL_post_event"] = self.build_GL
        self.registered_events["GP_post_event"] = self.build_GL_GP_SP
        self.registered_events["generic_fa_pre_event"] = self.set_generic_fa
        self.registered_events["generic_fa2_pre_event"] = self.set_generic_fa
        self.registered_events["generic_fa3_pre_event"] = self.set_generic_fa
        self.registered_events["generic_fa_rob_post_event"] = self.add_generic_fa
        self.registered_events["generic_fa_rob_s_post_event"] = self.add_generic_fa
        self.registered_events["generic_pos_pre_event"] = self.set_generic_pos
        self.registered_events["gl_ending_pre_event"] = self.set_gl_class
        self.registered_events["xdg_xmg_ending_pre_event"] = self.set_gl_class
        self.registered_events["gp_ending_pre_event"] = self.set_gp_class
        
        
        
        
        ## fatty acyl
        self.registered_events["fa_pre_event"] = self.set_fa
        self.registered_events["fa_post_event"] = self.build_fa
        self.registered_events["fatty_acid_post_event"] = self.set_fatty_acid
        self.registered_events["acetyl_pre_event"] = self.set_acetyl
        
        self.registered_events["acid_single_type_pre_event"] = self.set_fatty_acyl_type
        self.registered_events["ol_ending_pre_event"] = self.set_fatty_acyl_type
        self.registered_events["double_bond_position_pre_event"] = self.set_double_bond_information
        self.registered_events["double_bond_position_post_event"] = self.add_double_bond_information
        self.registered_events["db_number_post_event"] = self.set_double_bond_position
        self.registered_events["cistrans_post_event"] = self.set_cistrans
        self.registered_events["acid_type_double_post_event"] = self.check_db
        
        ## lengths
        self.registered_events["functional_length_pre_event"] = self.reset_length
        self.registered_events["fatty_length_pre_event"] = self.reset_length
        self.registered_events["functional_length_post_event"] = self.set_functional_length
        self.registered_events["fatty_length_post_event"] = self.set_fatty_length
        
        ## numbers
        self.registered_events["notation_specials_pre_event"] = self.special_number
        self.registered_events["notation_last_digit_pre_event"] = self.last_number
        self.registered_events["notation_second_digit_pre_event"] = self.second_number
        
        ## functional groups
        self.registered_events["functional_group_pre_event"] = self.set_functional_group
        self.registered_events["functional_group_post_event"] = self.add_functional_group
        self.registered_events["functional_pos_pre_event"] = self.set_functional_pos
        self.registered_events["functional_position_pre_event"] = self.set_functional_position
        self.registered_events["functional_group_type_pre_event"] = self.set_functional_type
        
        ## cyclo / epoxy
        self.registered_events["cyclo_position_pre_event"] = self.set_functional_group
        self.registered_events["cyclo_position_post_event"] = self.rearrange_cycle
        self.registered_events["epoxy_pre_event"] = self.set_functional_group
        self.registered_events["epoxy_post_event"] = self.add_epoxy
        self.registered_events["cycle_pre_event"] = self.set_cycle
        self.registered_events["methylene_post_event"] = self.set_methylene
        self.registered_events["furan_pre_event"] = self.set_furan
        
        ## dioic
        self.registered_events["dioic_pre_event"] = self.set_functional_group
        self.registered_events["dioic_post_event"] = self.set_dioic
        self.registered_events["dioic_acid_pre_event"] = self.set_fatty_acyl_type
        self.registered_events["dial_post_event"] = self.set_dial
        
        
        ## prosta
        self.registered_events["prosta_pre_event"] = self.set_prosta
        self.registered_events["prosta_post_event"] = self.add_cyclo
        self.registered_events["reduction_pre_event"] = self.set_functional_group
        self.registered_events["reduction_post_event"] = self.reduction
        self.registered_events["homo_post_event"] = self.homo
        
        
        ## recursion
        self.registered_events["recursion_description_pre_event"] = self.set_recursion
        self.registered_events["recursion_description_post_event"] = self.add_recursion
        self.registered_events["recursion_pos_pre_event"] = self.set_recursion_pos
        self.registered_events["yl_ending_pre_event"] = self.set_yl_ending
        self.registered_events["acetic_acid_post_event"] = self.set_acetic_acid
        self.registered_events["acetic_recursion_pre_event"] = self.set_recursion
        self.registered_events["acetic_recursion_post_event"] = self.add_recursion
        
        self.registered_events["hydroxyl_number_pre_event"] = self.add_hydroxyl
        self.registered_events["ol_pre_event"] = self.setup_hydroxyl
        self.registered_events["ol_post_event"] = self.add_hydroxyls
        self.registered_events["ol_pos_post_event"] = self.set_yl_ending
        
        
        
        ## wax esters
        self.registered_events["wax_ester_pre_event"] = self.set_recursion
        self.registered_events["wax_ester_post_event"] = self.add_wax_ester
        self.registered_events["ate_post_event"] = self.set_ate
        self.registered_events["isoprop_post_event"] = self.set_iso
        self.registered_events["isobut_post_event"] = self.set_iso
        
        ## CoA
        self.registered_events["CoA_post_event"] = self.set_coa
        self.registered_events["methyl_pre_event"] = self.set_methyl
        
        ## CAR
        self.registered_events["CAR_pre_event"] = self.set_car
        self.registered_events["CAR_post_event"] = self.add_car
        
        ## amine
        self.registered_events["ethanolamine_post_event"] = self.add_ethanolamine
        self.registered_events["amine_n_pre_event"] = self.set_recursion
        self.registered_events["amine_n_post_event"] = self.add_amine
        self.registered_events["amine_post_event"] = self.add_amine_name
        
        ## functional group position summary
        self.registered_events["fg_pos_summary_pre_event"] = self.set_functional_group
        self.registered_events["fg_pos_summary_post_event"] = self.add_summary
        self.registered_events["func_stereo_pre_event"] = self.add_func_stereo
        
        
        
    def reset_lipid(self, node):
        self.level = LipidLevel.ISOMERIC_SUBSPECIES
        self.headgroup = ""
        self.fatty_acyl_stack = []
        self.fatty_acyl_list = []
        self.fatty_acyl_dict = {}
        self.tmp = {}
        self.fatty_acyl = None
        
        #self.debug = "full"
        
        
        
    def set_gp_class(self, node):
        gp_class = node.get_text()
        if gp_class == 'sn-glycero-3-phosphocholine': self.headgroup = "PC"
        elif gp_class == 'sn-glycero-3-phosphoethanolamine': self.headgroup = "PE"
        elif gp_class == 'sn-glycero-3-phosphoserine': self.headgroup = "PS"
        
        
        
    def set_gl_class(self, node):
        gl_class = node.get_text()
        if gl_class in {'glycerol', 'rac-glycerol', 'sn-glycerol'}: self.tmp["gl_class"] = ""
        elif gl_class == '3-o-β-d-galactosyl-sn-glycerol': self.tmp["gl_class"] = "MG"
        elif gl_class in {'3-o-(6\'-o-α-d-galactosyl-β-d-galactosyl)-sn-glycerol', '3-o-(α-d-galactosyl1-6)-β-d-galactosyl-sn-glycerol'}: self.tmp["gl_class"] = "DG"
        elif gl_class == '3-(6\'-sulfo-α-d-quinovosyl)-sn-glycerol': self.tmp["gl_class"] = "SQ"
        
        
    def set_fa(self, node):
        self.fatty_acyl_stack.append(FattyAcid("FA"))
        self.tmp["fa1"] = {"fa_type": ""}
        


    def build_FA(self, node):
        self.headgroup = self.tmp["fa1"]["fa_type"]
        self.fatty_acyl_list = self.fatty_acyl_stack
        
        
    
    def build_fa(self, node):        
        if "cyclo_yl" in self.tmp:
            self.tmp["fg_pos"] = [[1, ""], [self.tmp["cyclo_len"], ""]]
            self.add_cyclo(node)
            del self.tmp["cyclo_yl"]
            del self.tmp["cyclo_len"]
            
        
        if "post_adding" in self.tmp:
            def add_position(func_group, pos):
                func_group.position += func_group.position >= pos
                if type(func_group) == Cycle:
                    func_group.start += func_group.start >= pos
                    func_group.end += func_group.end >= pos
                for fg_name, fg_list in func_group.functional_groups.items():
                    for fg in fg_list:
                        add_position(fg, pos)
            curr_fa = self.fatty_acyl_stack[-1]
            curr_fa.num_carbon += len(self.tmp["post_adding"])
            for pos in self.tmp["post_adding"]:
                add_position(curr_fa, pos)
                if type(curr_fa.double_bonds) == dict:
                    curr_fa.double_bonds = {(k + (k >= pos)): v for k, v in curr_fa.double_bonds.items()}
        
        if type(self.fatty_acyl_stack[-1].double_bonds) == dict and len(self.fatty_acyl_stack[-1].double_bonds) > 0:
            if sum(len(ct) > 0 for p, ct in self.fatty_acyl_stack[-1].double_bonds.items()) != len(self.fatty_acyl_stack[-1].double_bonds):
                self.set_lipid_level(LipidLevel.STRUCTURAL_SUBSPECIES)
                
        
        
        
    def set_acetyl(self, node):
        self.fatty_acyl_stack[-1] = FattyAcid("FA", num_carbon = 2)
        
        
        
    def build_GL(self, node):
        self.headgroup = self.tmp["gl_class"]
        del self.tmp["gl_class"]
        del self.tmp["generic_pos"]
        if len(self.fatty_acyl_dict) == 1: self.headgroup += "MG"
        if len(self.fatty_acyl_dict) == 2: self.headgroup += "DG"
        if len(self.fatty_acyl_dict) == 3: self.headgroup += "TG"
        self.build_GL_GP_SP(None)

        
        
        
    def build_GL_GP_SP(self, node):
        lipid_class = get_class(self.headgroup)
        if lipid_class == UNDEFINED_LIPID_CLASS:
            raise LipidParsingException("Unknown lipid class '%s'." % self.headgroup)
        
        self.fatty_acyl_list = [FattyAcid("FA%i" % i) for i in range(1, all_lipids[lipid_class]["max_fa"] + 1)]

        for fa_pos, fa in self.fatty_acyl_dict.items():
            self.fatty_acyl_list[fa_pos - 1] = fa
        
        
        
    def set_generic_fa(self, node):
        self.tmp["generic_pos"] = []
        
        
        
    def set_generic_pos(self, node):
        self.tmp["generic_pos"].append(int(node.get_text()))
        
       
       
    def add_generic_fa(self, node):
        if len(self.fatty_acyl_stack) != 1:
            raise LipidParsingException("Unable to parse lipid.")
                
        fa = self.fatty_acyl_stack.pop()
        
        fa_type = self.tmp["fa1"]["fa_type"]
        if (len(fa_type) == 2 and fa_type[-2:] == "yl") or (len(fa_type) >= 3 and fa_type[-3:] != "oyl"):
            fa.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMANYL
            if type(fa.double_bonds) == dict and 1 in fa.double_bonds:
                del fa.double_bonds[1]
                fa.lipid_FA_bond_type = LipidFaBondType.ETHER_PLASMENYL
                
        for pos in self.tmp["generic_pos"]:
            fa_copy = fa.copy()
            fa_copy.name = "FA%i" % pos
            self.fatty_acyl_dict[pos] = fa_copy
        self.tmp["generic_pos"] = []
        
        
        
    def set_car(self, node):
        self.tmp["fg_pos"] = []
        self.tmp["fg_type"] = ""
        
       
       
    def set_furan(self, node):
        self.tmp["furan"] = True
        
        
        
    def add_ol(self, node):
        self.tmp["fa1"]["fa_type"] = "FOH"
        
        
        
    def homo(self, node):
        self.tmp["post_adding"] = list(p[0] for p in self.tmp["fg_pos"])
        
        
        
    def add_summary(self, node):
        fa_i = "fa%i" % len(self.fatty_acyl_stack)
        self.tmp[fa_i]["fg_pos_summary"] = {k: v.upper() for k, v in self.tmp["fg_pos"]}
        
        
        
    def set_acetic_acid(self, node):
        self.fatty_acyl_stack[-1].num_carbon += 2
        self.tmp["fa1"]["fa_type"] = "FA"
        
        
        
    def add_func_stereo(self, node):
        self.tmp["fg_pos"][-1][1] = node.get_text()
        
        
        
    def set_yl_ending(self, node):
        l = int(node.get_text()) - 1
        if l == 0: return
    
        curr_fa = self.fatty_acyl_stack[-1]
        if l == 1:
            fname = "Me"
            fg = get_functional_group(fname)
        elif l == 2:
            fname = "Et"
            fg = get_functional_group(fname)
        else:
            fa = FattyAcid("FA", num_carbon = l)
            # shift functional groups
            for fg_name, fg_list in curr_fa.functional_groups.items():
                remove_item = []
                for i, func_group in enumerate(fg_list):
                    if func_group.position <= l:
                        remove_item.append(i)
                        if fg_name not in fa.functional_groups: fa.functional_groups[fg_name] = []
                        func_group.position = l + 1 - func_group.position
                        fa.functional_groups[fg_name].append(func_group)
                for i in remove_item[::-1]:
                    del curr_fa.functional_groups[fg_name][i]
            curr_fa.functional_groups = {k: v for k, v in curr_fa.functional_groups.items() if len(v) > 0}
            
            #shift double bonds
            if type(curr_fa.double_bonds) == dict:
                fa.double_bonds = {l + 1 - k: v for k, v in curr_fa.double_bonds.items() if k <= l}
                for k in list(k for k in curr_fa.double_bonds if k <= l): del curr_fa.double_bonds[k]
            
            fname = "cc"
            fg = CarbonChain(fa)
            
        curr_fa.num_carbon -= l
        fg.position = l
        curr_fa.shift_positions(-l)
        if fname not in curr_fa.functional_groups: curr_fa.functional_groups[fname] = []
        curr_fa.functional_groups[fname].append(fg)
        
        
        
    def set_methylene(self, node):
        self.tmp["fg_type"] = "methylene"
        if len(self.tmp["fg_pos"]) > 1:
            if self.tmp["fg_pos"][0][0] < self.tmp["fg_pos"][1][0]: self.tmp["fg_pos"][1][0] += 1
            elif self.tmp["fg_pos"][0][0] > self.tmp["fg_pos"][1][0]: self.tmp["fg_pos"][0][0] += 1
            self.fatty_acyl_stack[-1].num_carbon += 1
            self.tmp["add_methylene"] = True
     
    
    def check_db(self, node):
        fa_i = "fa%i" % len(self.fatty_acyl_stack)
        curr_fa = self.fatty_acyl_stack[-1]
        if "fg_pos_summary" in self.tmp[fa_i]:
            if type(curr_fa.double_bonds) != dict: curr_fa.double_bonds = {}
            for k, v in self.tmp[fa_i]["fg_pos_summary"].items():
                if v in {"E", "Z", ""} and k > 0 and k not in curr_fa.double_bonds: curr_fa.double_bonds[k] = v
                    
        
        
    def add_car(self, node):
        self.tmp["fa1"]["fa_type"] = "CAR"
        
        
        
    def add_hydroxyl(self, node):
        self.tmp["hydroxyl_pos"].append(int(node.get_text()))
    
    
    
    def set_coa(self, node):
        self.tmp["fa1"]["fa_type"] = "CoA"
        
        
        
    def set_methyl(self, node):
        self.fatty_acyl_stack[-1].num_carbon += 1
    
    
    
    def setup_hydroxyl(self, node):
        self.tmp["hydroxyl_pos"] = []
        
        
        
    def set_iso(self, node):
        curr_fa = self.fatty_acyl_stack[-1]
        curr_fa.num_carbon -= 1
        fg = get_functional_group("Me")
        fg.position = 2
        if "Me" not in curr_fa.functional_groups: curr_fa.functional_groups["Me"] = []
        curr_fa.functional_groups["Me"].append(fg)
        
        
        
    def set_ate(self, node):
        self.fatty_acyl_stack[-1].num_carbon += ate[node.get_text()]
        self.tmp["fa1"]["fa_type"] = "WE"
        
        
        
    def add_amine(self, node):
        fa_i = "fa%i" % len(self.fatty_acyl_stack)
        fa = self.fatty_acyl_stack.pop()
        
        fa.name += "1"
        self.fatty_acyl_stack[-1].name += "2"
        fa.lipid_FA_bond_type = LipidFaBondType.AMINE
        self.fatty_acyl_stack.insert(0, fa)
        
        
        
    def add_amine_name(self, node):
        self.tmp["fa1"]["fa_type"] = "NA"
        
        
        
    def add_ethanolamine(self, node):
        self.tmp["fa1"]["fa_type"] = "NAE"
        
        
    
    def add_hydroxyls(self, node):
        if len(self.tmp["hydroxyl_pos"]) > 1:
            fg_oh = get_functional_group("OH")
            for pos in sorted(self.tmp["hydroxyl_pos"], reverse = True)[:-1]:
                fg_insert = fg_oh.copy()
                fg_insert.position = pos
                if "OH" not in self.fatty_acyl_stack[-1].functional_groups: self.fatty_acyl_stack[-1].functional_groups["OH"] = []
                self.fatty_acyl_stack[-1].functional_groups["OH"].append(fg_insert)
        
        
        
    def set_double_bond_position(self, node):
        fa_i = "fa%i" % len(self.fatty_acyl_stack)
        pos = int(node.get_text())
        self.tmp[fa_i]["db_position"] = pos - (sum([1 for p in self.tmp["reduction"] if p < pos]) if "reduction" in self.tmp else 0)
        
        
        
    def set_recursion_pos(self, node):
        fa_i = "fa%i" % len(self.fatty_acyl_stack)
        self.tmp[fa_i]["recursion_pos"] = int(node.get_text())
        
        
        
    def set_recursion(self, node):
        self.tmp["fg_pos"] = []
        self.tmp["fg_type"] = ""
        self.fatty_acyl_stack.append(FattyAcid("FA"))
        fa_i = "fa%i" % len(self.fatty_acyl_stack)
        self.tmp[fa_i] = {"recursion_pos": 0}
        
        
        
    def add_recursion(self, node):
        fa_i = "fa%i" % len(self.fatty_acyl_stack)
        pos = self.tmp[fa_i]["recursion_pos"]
        
        
        fa = self.fatty_acyl_stack.pop()
        fa.position = pos
        curr_fa = self.fatty_acyl_stack[-1]
        
        if "cyclo_yl" in self.tmp:
            fname = "cyclo"
            del self.tmp["cyclo_yl"]
        else:
            fname = self.tmp["fa1"]["fa_type"]
        if fname not in curr_fa.functional_groups: curr_fa.functional_groups[fname] = []
        curr_fa.functional_groups[fname].append(fa)
        
        self.tmp["added_func_group"] = True
        
        
        
    def add_wax_ester(self, node):
        fa_i = "fa%i" % len(self.fatty_acyl_stack)
        fa = self.fatty_acyl_stack.pop()
        
        fa.name += "1"
        fa.lipid_FA_bond_type = LipidFaBondType.AMINE
        self.fatty_acyl_stack[-1].name += "2"
        self.fatty_acyl_stack.insert(0, fa)
        
        
        
    def set_cycle(self, node):
        self.tmp["cyclo"] = True
        
        
        
    def set_fatty_acid(self, node):
        def switch_position(func_group, switch):
            func_group.position = switch - func_group.position
            for fg_name, fg_list in func_group.functional_groups.items():
                for fg in fg_list:
                    switch_position(fg, switch)
            
        curr_fa = self.fatty_acyl_stack[-1]
        
        if "noyloxy" in curr_fa.functional_groups:
            if self.tmp["fa1"]["fa_type"] == "FA": self.tmp["fa1"]["fa_type"] = "FAHFA"
            
            while len(curr_fa.functional_groups["noyloxy"]) > 0:
                fa = curr_fa.functional_groups["noyloxy"].pop()
            
                
                acyl = AcylAlkylGroup(fa)
                acyl.position = fa.position
                
                if "acyl" not in curr_fa.functional_groups: curr_fa.functional_groups["acyl"] = []
                curr_fa.functional_groups["acyl"].append(acyl)
            del curr_fa.functional_groups["noyloxy"]
            
            
        elif "nyloxy" in curr_fa.functional_groups:
            while len(curr_fa.functional_groups["nyloxy"]) > 0:
                fa = curr_fa.functional_groups["nyloxy"].pop()
            
                
                alkyl = AcylAlkylGroup(fa, alkyl = True)
                alkyl.position = fa.position
                
                if "alkyl" not in curr_fa.functional_groups: curr_fa.functional_groups["alkyl"] = []
                curr_fa.functional_groups["alkyl"].append(alkyl)
            del curr_fa.functional_groups["nyloxy"]
            
                
        elif sum([k[-2:] == "yl" for k in curr_fa.functional_groups]) > 0:
            while True:
                try:
                    yl = [k for k in curr_fa.functional_groups if k[-2:] == "yl"][0]
                except Exception:
                    break
            
                while len(curr_fa.functional_groups[yl]) > 0:
                    fa = curr_fa.functional_groups[yl].pop()
                    
                    if "furan" in self.tmp:
                        self.tmp['cyclo_len'] = 4
                        fa.shift_positions(4)
                        for fg, fg_list in fa.functional_groups.items():
                            if fg not in curr_fa.functional_groups: curr_fa.functional_groups[fg] = fg_list
                            else: curr_fa.functional_groups[fg] += fg_list
                            
                        curr_fa.num_carbon += 4 + fa.num_carbon
                        
                        if type(curr_fa.double_bonds) == int: curr_fa.double_bonds = {}
                        if type(fa.double_bonds) == dict:
                            for pos, ez in fa.double_bonds.items():
                                curr_fa.double_bonds[pos + cyclo_len] = ez
                        curr_fa.double_bonds[1] = "E"
                        curr_fa.double_bonds[3] = "E"
                        self.tmp["cyclo_yl"] = True
                        
                        
                    
                    elif "cyclo" in self.tmp:
                        cyclo_len = curr_fa.num_carbon
                        self.tmp['cyclo_len'] = cyclo_len
                        switch_position(curr_fa, 2 + cyclo_len)
                        if type(curr_fa.double_bonds) == dict: curr_fa.double_bonds = {2 + cyclo_len - k: v for k, v in curr_fa.double_bonds.items()}
                        fa.shift_positions(cyclo_len)
                        
                        for fg, fg_list in fa.functional_groups.items():
                            if fg not in curr_fa.functional_groups: curr_fa.functional_groups[fg] = fg_list
                            else: curr_fa.functional_groups[fg] += fg_list
                            
                        
                        curr_fa.num_carbon = cyclo_len + fa.num_carbon
                        
                        if type(curr_fa.double_bonds) == int: curr_fa.double_bonds = {}
                        if type(fa.double_bonds) == dict:
                            for pos, ez in fa.double_bonds.items():
                                curr_fa.double_bonds[pos + cyclo_len] = ez
                        self.tmp["cyclo_yl"] = True
                        
                    else:
                        ## add carbon chains here here
                        ## special chains: i.e. ethyl, methyl
                        fg_name = ""
                        if (fa.double_bonds if type(fa.double_bonds) == int else len(fa.double_bonds)) == 0 and len(fa.functional_groups) == 0:
                            if fa.num_carbon == 1:
                                fg_name = "Me"
                                fg = get_functional_group(fg_name)
                            
                            elif fa.num_carbon == 2:
                                fg_name = "Et"
                                fg = get_functional_group(fg_name)
                                
                            if len(fg_name) > 0:
                                fg.position = fa.position
                                if fg_name not in curr_fa.functional_groups: curr_fa.functional_groups[fg_name] = []
                                curr_fa.functional_groups[fg_name].append(fg)
                            
                        if len(fg_name) == 0:
                            cc = CarbonChain(fa, position = fa.position)
                            if "cc" not in curr_fa.functional_groups: curr_fa.functional_groups["cc"] = []
                            curr_fa.functional_groups["cc"].append(cc)
                            
                if "cyclo" in self.tmp: del self.tmp["cyclo"]
                del curr_fa.functional_groups[yl]
           
           
        if "cyclo" in curr_fa.functional_groups:
            fa = curr_fa.functional_groups["cyclo"][0]
            del curr_fa.functional_groups["cyclo"]
            start_pos, end_pos = curr_fa.num_carbon + 1, curr_fa.num_carbon + (self.tmp["cyclo_len"] if 'cyclo_len' in self.tmp else 5)
            fa.shift_positions(start_pos - 1)
            
            if "cy" in curr_fa.functional_groups:
                for cy in curr_fa.functional_groups["cy"]:
                    cy.shift_positions(start_pos - 1)
            
            for fg, fg_list in fa.functional_groups.items():
                if fg not in curr_fa.functional_groups: curr_fa.functional_groups[fg] = fg_list
                else: curr_fa.functional_groups[fg] += fg_list
            
            curr_fa.num_carbon += fa.num_carbon
            
            if type(curr_fa.double_bonds) == int: curr_fa.double_bonds = {}
            if type(fa.double_bonds) == dict:
                for pos, ez in fa.double_bonds.items():
                    curr_fa.double_bonds[pos + start_pos - 1] = ez
                    
            self.tmp["fg_pos"] = [[start_pos, ""], [end_pos, ""]]
            self.add_cyclo(node)
            
            if "cyclo_len" in self.tmp: del self.tmp["cyclo_len"]
            if "cyclo" in self.tmp: del self.tmp["cyclo"]
            if "furan" in self.tmp: del self.tmp["furan"]
            
            
        elif "cyclo" in self.tmp:
            self.tmp["fg_pos"] = [[1, ""], [curr_fa.num_carbon, ""]]
            self.add_cyclo(node)
            del self.tmp["cyclo"]
        
            
        
    def set_double_bond_information(self, node):
        fa_i = "fa%i" % len(self.fatty_acyl_stack)
        self.tmp[fa_i]["db_position"] = 0
        self.tmp[fa_i]["db_cistrans"] = ""
        
        
        
    def reduction(self, node):
        self.fatty_acyl_stack[-1].num_carbon -= len(self.tmp["fg_pos"])
        self.tmp["reduction"] = [p[0] for p in self.tmp["fg_pos"]]
        
        
        
    def add_double_bond_information(self, node):
        fa_i = "fa%i" % len(self.fatty_acyl_stack)
        pos = self.tmp[fa_i]["db_position"]
        cistrans = self.tmp[fa_i]["db_cistrans"].upper()
        if cistrans == "" and "fg_pos_summary" in self.tmp[fa_i] and pos in self.tmp[fa_i]["fg_pos_summary"]: cistrans = self.tmp[fa_i]["fg_pos_summary"][pos]
        if pos == 0: return
        
        del self.tmp[fa_i]["db_position"]
        del self.tmp[fa_i]["db_cistrans"]
        if type(self.fatty_acyl_stack[-1].double_bonds) == int: self.fatty_acyl_stack[-1].double_bonds = {}
        
        if pos not in self.fatty_acyl_stack[-1].double_bonds or len(self.fatty_acyl_stack[-1].double_bonds[pos]) == 0:
            self.fatty_acyl_stack[-1].double_bonds[pos] = cistrans
        
        
        
    def set_dioic(self, node):
        self.tmp["fa1"]["fa_type"] = "FA"
        
        pos = self.tmp["fg_pos"][1][0] if len(self.tmp["fg_pos"]) == 2 else self.fatty_acyl_stack[-1].num_carbon
        self.fatty_acyl_stack[-1].num_carbon -= 1
        func_group = get_functional_group("COOH")
        func_group.position = pos - 1
        if "COOH" not in self.fatty_acyl_stack[-1].functional_groups: self.fatty_acyl_stack[-1].functional_groups["COOH"] = []
        self.fatty_acyl_stack[-1].functional_groups["COOH"].append(func_group)
        
        
        
    def set_cistrans(self, node):
        self.tmp["fa%i" % len(self.fatty_acyl_stack)]["db_cistrans"] = node.get_text()
        
        
        
    def set_fatty_acyl_type(self, node):
        t = node.get_text()
        
        if t[-2:] == "ol": self.tmp["fa1"]["fa_type"] = "FOH"
        elif t in {"noic acid", "nic acid", "dioic_acid"}: self.tmp["fa1"]["fa_type"] = "FA"
        elif t in {"nal", "dial"}: self.tmp["fa1"]["fa_type"] = "FAL"
        elif t in {"acetate", "noate", "nate"}: self.tmp["fa1"]["fa_type"] = "WE"
        elif t == "ne":
            self.tmp["fa1"]["fa_type"] = "HC"
            self.fatty_acyl_stack[-1].lipid_FA_bond_type = LipidFaBondType.AMINE
            
        else: self.tmp["fa1"]["fa_type"] = t
        
        
        
    def set_lipid_level(self, level):
        self.level = self.level if self.level.value < level.value else level
        
        
        
    def set_dial(self, node):
        curr_fa = self.fatty_acyl_stack[-1]
        pos = curr_fa.num_carbon
        fg = get_functional_group("oxo")
        fg.position = pos
        if "oxo" not in curr_fa.functional_groups: curr_fa.functional_groups["oxo"] = []
        curr_fa.functional_groups["oxo"].append(fg)
        
        
    
    def reset_length(self, node):
        self.tmp["length"] = 0
    
    
    
    def set_functional_length(self, node):
        if self.tmp["length"] != len(self.tmp["fg_pos"]):
            raise LipidException("Length of functional group '%i' does not match with number of its positions '%i'" % (self.tmp["length"], len(self.tmp["fg_pos"])))
    
    
    
    def set_fatty_length(self, node):
        self.fatty_acyl_stack[-1].num_carbon += self.tmp["length"]
        
        
        
    def rearrange_cycle(self, node):
        if "post_adding" in self.tmp:
            self.fatty_acyl_stack[-1].num_carbon += len(self.tmp["post_adding"])
            del self.tmp["post_adding"]
            
        curr_fa = self.fatty_acyl_stack[-1]
        start = self.tmp["fg_pos"][0][0]
        if "cy" in curr_fa.functional_groups:
            for cy in curr_fa.functional_groups["cy"]:
                shift = start - cy.position
                if shift == 0: continue
                cy.rearrange_functional_groups(curr_fa, shift)
        
        
        
    def add_cyclo(self, node):
        start = self.tmp["fg_pos"][0][0]
        end = self.tmp["fg_pos"][1][0]
        
            
        cyclo_db = None
        # check double bonds
        if type(self.fatty_acyl_stack[-1].double_bonds) == dict and len(self.fatty_acyl_stack[-1].double_bonds) > 0:
            cyclo_db = {db_pos: val for db_pos, val in self.fatty_acyl_stack[-1].double_bonds.items() if start <= db_pos <= end}
            
            for pos in cyclo_db:
                del self.fatty_acyl_stack[-1].double_bonds[pos]
                
        # check functional_groups
        cyclo_fg, remove_list, curr_fa = {}, set(), self.fatty_acyl_stack[-1]
        
        if "noyloxy" in curr_fa.functional_groups:
            remove_item = []
            for i, func_group in enumerate(curr_fa.functional_groups["noyloxy"]):
                if start <= func_group.position <= end:
                    cc = CarbonChain(func_group, position = func_group.position)
                    
                    if "cc" not in curr_fa.functional_groups: curr_fa.functional_groups["cc"] = []
                    curr_fa.functional_groups["cc"].append(cc)
                    remove_item.append(i)
                    
            for i in remove_item[::-1]: del curr_fa.functional_groups["noyloxy"][i]
            if len(curr_fa.functional_groups["noyloxy"]) == 0: remove_list.add("noyloxy")
            
        for fg, fg_list in curr_fa.functional_groups.items():
            remove_item = []
            
            for i, func_group in enumerate(fg_list):
                if start <= func_group.position <= end:
                    if fg not in cyclo_fg: cyclo_fg[fg] = []
                    cyclo_fg[fg].append(func_group)
                    remove_item.append(i)
                    
            for i in remove_item[::-1]: del curr_fa.functional_groups[fg][i]
            if len(fg_list) == 0: remove_list.add(fg)
            
        for fg in remove_list: del curr_fa.functional_groups[fg]
        cy = end - start + 1
        
        bridge_chain = []
        if "furan" in self.tmp:
            bridge_chain.append(Element.O)
            cy += 1

        cycle = Cycle(cy, start = start, end = end, double_bonds = cyclo_db, functional_groups = cyclo_fg, bridge_chain = bridge_chain)
        if "cy" not in self.fatty_acyl_stack[-1].functional_groups: self.fatty_acyl_stack[-1].functional_groups["cy"] = []
        self.fatty_acyl_stack[-1].functional_groups["cy"].append(cycle)
        
        
        
    def add_epoxy(self, node):
        self.tmp["fg_pos"] = self.tmp["fg_pos"][:1]
        self.tmp["fg_type"] = "Epoxy"
        
    
    
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
        self.tmp["fg_pos"] = [[8 - minus_pos, ""], [12 - minus_pos, ""]]
        self.tmp["fg_type"] = "cy"
        
        
        
    def add_functional_group(self, node):
        if "added_func_group" in self.tmp: 
            del self.tmp["added_func_group"]
            return
        
        elif "add_methylene" in self.tmp: 
            del self.tmp["add_methylene"]
            self.add_cyclo(node)
            return
        
        t = self.tmp["fg_type"]
        
        if t != "acetoxy":
            if t not in func_groups: raise LipidException("Unknown functional group: '%s'" % t)
            t = func_groups[t]
            if len(t) == 0: return
            fg = get_functional_group(t)
        else:
            fg = AcylAlkylGroup(FattyAcid("O", num_carbon = 2))
            
            
        
        if t not in self.fatty_acyl_stack[-1].functional_groups: self.fatty_acyl_stack[-1].functional_groups[t] = []
        for pos in self.tmp["fg_pos"]:
            fg_insert = fg.copy()
            fg_insert.position = pos[0] - (sum([1 for p in self.tmp["reduction"] if p < pos[0]]) if "reduction" in self.tmp else 0)
            self.fatty_acyl_stack[-1].functional_groups[t].append(fg_insert)
        
        
        
    def set_functional_position(self, node):
        self.tmp["fg_pos"].append([0, ""])
        
        
        
    def set_functional_pos(self, node):
        self.tmp["fg_pos"][-1][0] = int(node.get_text())
        
        
        
    def set_functional_type(self, node):
        self.tmp["fg_type"] = node.get_text()
        
        
        
    def build_lipid(self, node):
        lipid_level_class = None
        if self.level == LipidLevel.ISOMERIC_SUBSPECIES: lipid_level_class = LipidIsomericSubspecies
        if self.level == LipidLevel.STRUCTURAL_SUBSPECIES: lipid_level_class = LipidStructuralSubspecies
        if self.level == LipidLevel.MOLECULAR_SUBSPECIES: lipid_level_class = LipidMolecularSubspecies
        if self.level == LipidLevel.SPECIES: lipid_level_class = LipidSpecies
        
        headgroup = HeadGroup(self.headgroup)
        
        self.content = LipidAdduct()
        self.content.lipid = lipid_level_class(headgroup, self.fatty_acyl_list)
            
