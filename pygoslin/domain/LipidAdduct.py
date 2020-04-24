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


from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.Element import *

class LipidAdduct:

    def __init__(self):
        self.lipid = None
        self.adduct = None
        self.fragment = None
        self.sum_formula = None
        
        
        
    def get_lipid_string(self, level = None):
        lipid_name = []
        
        if self.lipid != None: lipid_name.append(self.lipid.get_lipid_string(level))
        else: return ""
        
        if level == None or level not in {LipidLevel.CLASS, LipidLevel.CATEGORY}:
            if self.adduct != None: lipid_name.append(self.adduct.get_lipid_string())
        
        return "".join(lipid_name)
        
        
        
    def get_lipid_fragment_string(self, level = None):
        lipid_name = []
        
        if self.lipid != None: lipid_name.append(self.lipid.get_lipid_string(level))
        else: return ""
        
        if level == None or level not in {LipidLevel.CLASS, LipidLevel.CATEGORY}:
            if self.adduct != None: lipid_name.append(self.adduct.get_lipid_string())
        
            if self.fragment != None:
                lipid_name.append(" - ")
                lipid_name.append(self.fragment.get_lipid_string())
        
        return "".join(lipid_name)
    
    
    
    def get_elements(self):
        elements = {e: 0 for e in Element}
        
        if self.lipid != None:
            lipid_elements = self.lipid.get_elements()
            for e in Element: elements[e] += lipid_elements[e]
                
        if self.adduct != None:
            adduct_elements = self.adduct.get_elements()
            for e in Element:
                elements[e] += adduct_elements[e]
        return elements
    
    
    
    def get_mass(self):
        elements = self.get_elements()
        charge = charge = self.adduct.get_charge() if self.adduct != None else 0
        mass = sum([element_masses[e] * elements[e] for e in Element])
        if charge != 0: mass = (mass - charge * electron_rest_mass) / abs(charge)
            
        return mass
    
    
    
    
    def get_sum_formula(self):
        return compute_sum_formula(self.get_elements())
