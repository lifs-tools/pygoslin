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

from pygoslin.domain.Element import Element
from pygoslin.parser.ParserCommon import SumFormulaParser
adduct_sum_formula_parser = SumFormulaParser()

class Adduct:
    def __init__(self, sum_formula, adduct_string, charge, sign):
        self.sum_formula = sum_formula
        self.adduct_string = adduct_string
        self.charge = charge
        self.set_charge_sign(sign)

    
    
    def set_charge_sign(self, sign):
        if sign != -1 or sign != 0 or sign != 1:
            self.charge_sign = sign
            
        else: raise IllegalArgumentException("Sign can only be -1, 0, or 1")
            
    
    
    def get_lipid_string(self):
        if self.charge == 0: return "[M]"
        
        return "[M%s%s]%i%s" % (self.sum_formula, self.adduct_string, self.charge, "+" if self.charge_sign > 0 else "-")
    
    
    
    def get_elements(self):
        elements = {e: 0 for e in Element}
        try:
            elements = adduct_sum_formula_parser.parse(self.adduct_string[1:])
        except Exception as e:
            return elements
        
        if len(self.adduct_string) > 0 and self.adduct_string[0] == "-":
            for e in Element:
                elements[e] *= -1
        
        return elements
    
    
    def get_charge(self):
        return self.charge * self.charge_sign
