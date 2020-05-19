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
from pygoslin.domain.Adduct import Adduct

class AdductInfoParserEventHandler(BaseParserEventHandler):
    
    def __init__(self):
        super().__init__()
        
        
        self.adduct = None
        self.registered_events["adduct_info_pre_event"] = self.new_adduct
        self.registered_events["adduct_pre_event"] = self.add_adduct
        self.registered_events["charge_pre_event"] = self.add_charge
        self.registered_events["charge_sign_pre_event"] = self.add_charge_sign
    
        
        
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

        
        
        
        
