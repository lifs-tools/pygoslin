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

from pygoslin.domain.Element import Element, element_order, heavy_shortcut, heavy_to_regular
from pygoslin.parser.ParserCommon import SumFormulaParser
from pygoslin.domain.LipidExceptions import *


class Adduct:
    adduct_sum_formula_parser = SumFormulaParser()

    adducts = {"+H": {Element.H: 1},
               "+2H": {Element.H: 2},
               "+3H": {Element.H: 3},
               "+4H": {Element.H: 4},
               "-H": {Element.H: -1},
               "-2H": {Element.H: -2},
               "-3H": {Element.H: -3},
               "-4H": {Element.H: -4},
               "+H-H2O": {Element.H: -1, Element.O: -1},
               "+NH4": {Element.N: 1, Element.H: 4},
               "+Cl": {Element.Cl: 1},
               "+HCOO": {Element.H: 1, Element.C: 1, Element.O: 2},
               "+CH3COO": {Element.H: 3, Element.C: 2, Element.O: 2}
               }
    
    adduct_charges = {"+H": 1, "+2H": 2, "+3H": 3, "+4H": 4, "-H": -1,
                      "-2H": -2, "-3H": -3, "-4H": -4, "+H-H2O": 1,
                      "+NH4": 1, "+Cl": -1, "+HCOO": -1, "+CH3COO": -1
                      }
    
    def __init__(self, sum_formula, adduct_string, charge = 0, sign = 1):
        self.sum_formula = sum_formula
        self.adduct_string = adduct_string
        self.charge = charge
        self.set_charge_sign(sign)
        self.heavy_elements = {e: 0 for e in Element}
    
    
    def set_charge_sign(self, sign):
        if sign in {-1, 1}:
            self.charge_sign = sign
    
    
    def get_heavy_isotope_string(self):
        return "".join(["%s%i" % (heavy_shortcut[e], self.heavy_elements[e]) if self.heavy_elements[e] > 1 else heavy_shortcut[e] for e in element_order if self.heavy_elements[e] > 0])
            
    
    
    def get_lipid_string(self):
        if self.charge == 0: return "[M%s]" % self.get_heavy_isotope_string()
        
        return "[M%s%s%s]%i%s" % (self.sum_formula, self.get_heavy_isotope_string(), self.adduct_string, self.charge, "+" if self.charge_sign > 0 else "-")
    
    
    
    def get_elements(self):
        elements = {e: 0 for e in Element}
        
        for e, num in self.heavy_elements.items():
            if num > 0:
                elements[heavy_to_regular[e]] -= num
                elements[e] += num
        
        if len(self.adduct_string) > 0:
            if self.adduct_string in Adduct.adducts:
                if Adduct.adduct_charges[self.adduct_string] != self.get_charge():
                    raise ConstraintViolationException("Provided charge '%i' in contradiction to adduct '%s' charge '%i'." % (self.get_charge(), self.adduct_string, Adduct.adduct_charges[self.adduct_string]))
                    
                for k, v in Adduct.adducts[self.adduct_string].items():
                    elements[k] += v
                
            else:
                raise ConstraintViolationException("Adduct '%s' is unknown." % self.adduct_string)
        
        return elements
    
    
    def get_charge(self):
        return self.charge * self.charge_sign
