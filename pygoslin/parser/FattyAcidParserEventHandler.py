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

last_numbers = {'un': 1, 'do': 2, 'di': 2, 'tri': 3, 'buta': 4, 'but': 4, 'tetra': 4, 'penta': 5, 'pent': 5, 'hexa': 6, 'hex': 6, 'hepta': 7, 'hept': 7, 'octa': 8, 'oct': 8, 'nona': 9, 'non': 9}
second_numbers = {'deca': 10, 'dec': 10, 'cosa': 20, 'cos': 20, 'triaconta': 30, 'triacont': 30, 'tetraconta': 40, 'tetracont': 40}
special_numbers = {'buta': 4, 'deca': 10, 'eicosa': 20, 'heneicosa': 21, 'triaconta': 30, 'tetraconta': 40}

class FattyAcidParserEventHandler(BaseParserEventHandler):
    
    def __init__(self):
        super().__init__()
        
        self.registered_events["lipid_pre_event"] = self.reset_lipid
        self.registered_events["lipid_post_event"] = self.build_lipid
        
        self.registered_events["acid_single_type_pre_event"] = self.set_fatty_acyl_type
        self.registered_events["double_bond_position_pre_event"] = self.set_double_bond_information
        self.registered_events["double_bond_position_post_event"] = self.add_double_bond_information
        self.registered_events["db_number_post_event"] = self.set_double_bond_position
        self.registered_events["cistrans_post_event"] = self.set_cistrans
        
        ## lengths
        self.registered_events["methly_length_pre_event"] = self.reset_length
        self.registered_events["db_length_pre_event"] = self.reset_length
        self.registered_events["fatty_length_pre_event"] = self.reset_length
        self.registered_events["methly_length_post_event"] = self.set_methyl_length
        self.registered_events["db_length_post_event"] = self.set_db_length
        self.registered_events["fatty_length_post_event"] = self.set_fatty_length
        
        ## numbers
        self.registered_events["notation_specials_pre_event"] = self.special_number
        self.registered_events["notation_last_digit_pre_event"] = self.last_number
        self.registered_events["notation_second_digit_pre_event"] = self.second_number
        
        
    def reset_lipid(self, node):
        self.level = LipidLevel.ISOMERIC_SUBSPECIES
        self.lipid = None
        self.headgroup = ""
        self.current_fa = [FattyAcid("FA")]
        self.db_numbers = -1
        self.tmp = {"fa1": {}}
        
        
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
        
    def set_fatty_acyl_type(self, node):
        t = node.get_text()
        
        if t == "nol": self.headgroup = "FOH"
        elif t == "noic acid": self.headgroup = "FA"
        elif t == "nal": self.headgroup = "FAL"
        
        
    
    def reset_length(self, node):
        self.tmp["length"] = 0
    
    
    def set_methyl_length(self, node):
        pass 
    
    
    def set_fatty_length(self, node):
        self.current_fa[-1].num_carbon = self.tmp["length"]
        
        
    def set_db_length(self, node):
        if len(self.current_fa[-1].double_bonds) != self.tmp["length"]:
            raise LipidException("Double bond count does not match with number of double bond positions")
        
    
    def special_number(self, node):
        self.tmp["length"] = special_numbers[node.get_text()]
        
    def last_number(self, node):
        self.tmp["length"] += last_numbers[node.get_text()]
        
    def second_number(self, node):
        self.tmp["length"] += second_numbers[node.get_text()]
        
        
    def build_lipid(self, node):
        lipid_level_class = None
        if self.level == LipidLevel.ISOMERIC_SUBSPECIES: lipid_level_class = LipidIsomericSubspecies
        if self.level == LipidLevel.STRUCTURAL_SUBSPECIES: lipid_level_class = LipidStructuralSubspecies
        if self.level == LipidLevel.MOLECULAR_SUBSPECIES: lipid_level_class = LipidMolecularSubspecies
        if self.level == LipidLevel.SPECIES: lipid_level_class = LipidSpecies
        
        headgroup = HeadGroup(self.headgroup)
        
        self.lipid = LipidAdduct()
        self.lipid.lipid = lipid_level_class(headgroup, self.current_fa)
            
        self.content = self.lipid
