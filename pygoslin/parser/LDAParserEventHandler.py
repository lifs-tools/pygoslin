"""
MIT License

Copyright (c) 2026 Dominik Kopczynski   -   dominik.kopczynski {at} univie.ac.at
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

from pygoslin.parser.GoslinParserEventHandler import GoslinParserEventHandler
from pygoslin.domain.FunctionalGroup import get_functional_group
from pygoslin.domain.LipidLevel import LipidLevel


class LDAParserEventHandler(GoslinParserEventHandler):
    """
    Event handler for LDA.g4, the sphingolipid names written by
    LipidDataAnalyzer (LDA), e.g. 'Cer d18:1/n16:0' or 'HexCer t18:0/h24:0'.
    It reuses the Goslin handler for building the lipid and only
    registers the events of the LDA grammar.
    """

    # default LDA encodings, see lcbHydroxylationEncoding.txt and
    # faHydroxylationEncoding.txt in https://github.com/ThallingerLab/LDA2
    # LDA head group -> Goslin class name
    HEADGROUPS = {"SphBase": "SPB", "Cer1P": "CerP"}
    LCB_HYDROXYLS = {"m": 1, "d": 2, "t": 3, "q": 4, "p": 5, "h": 6}
    FA_HYDROXYLS = {"n": 0, "h": 1, "i": 2, "k": 3}

    def __init__(self):
        super().__init__()
        self.registered_events = {}

        self.registered_events["lipid_pre_event"] = self.reset_lipid
        self.registered_events["lipid_post_event"] = self.build_lipid

        self.registered_events["hg_lsl_pre_event"] = self.set_lda_head_group_name
        self.registered_events["hg_dsl_pre_event"] = self.set_lda_head_group_name

        self.registered_events["sl_species_pre_event"] = self.set_species_level
        self.registered_events["sl_unsorted_pre_event"] = self.set_molecular_subspecies_level

        self.registered_events["lcb_pre_event"] = self.new_lcb
        self.registered_events["lcb_post_event"] = self.clean_lcb
        self.registered_events["fa_pre_event"] = self.new_fa
        self.registered_events["fa_post_event"] = self.append_fa

        self.registered_events["lcb_hydroxyl_pre_event"] = self.add_lcb_hydroxyl
        self.registered_events["fa_hydroxyl_pre_event"] = self.add_fa_hydroxyl
        self.registered_events["db_count_pre_event"] = self.add_double_bonds
        self.registered_events["carbon_pre_event"] = self.add_carbon


    def set_lda_head_group_name(self, node):
        head_group = node.get_text()
        self.head_group = LDAParserEventHandler.HEADGROUPS.get(head_group, head_group)


    def add_hydroxyl_count(self, num_h):
        if num_h <= 0: return
        functional_group = get_functional_group("OH").copy()
        functional_group.count = num_h
        if "OH" not in self.current_fa.functional_groups: self.current_fa.functional_groups["OH"] = []
        self.current_fa.functional_groups["OH"].append(functional_group)


    def add_lcb_hydroxyl(self, node):
        num_h = LDAParserEventHandler.LCB_HYDROXYLS[node.get_text()]
        # the hydroxy group at C1 of a regular long chain base carries the head group
        if self.sp_regular_lcb(): num_h -= 1
        self.add_hydroxyl_count(num_h)


    def add_fa_hydroxyl(self, node):
        self.add_hydroxyl_count(LDAParserEventHandler.FA_HYDROXYLS[node.get_text()])
