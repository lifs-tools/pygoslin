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


from pygoslin.domain.Element import *
from pygoslin.domain.LipidLevel import LipidLevel

class FunctionalGroup:
    def __init__(self, name, position = -1, count = 1, stereochemistry = None, elements = None, functional_groups = None):
        self.name = name
        self.position = position
        self.count = count
        self.stereochemistry = stereochemistry
        self.elements = {e: 0 for e in Element} if elements == None else {k: v for k, v in elements.items()}
        self.functional_groups = functional_groups if functional_groups != None else {}
        
        
    def clone(self, fg):
        self.name = fg.name
        self.elements = {k: v for k, v in fg.elements.items()}
        
    def copy(self):
        return FunctionalGroup(self.name, position = self.position, count = self.count, stereochemistry = self.stereochemistry, elements = {k: v for k, v in self.elements.items()})
        
        
    def get_elements(self):
        elements = {k: v * self.count for k, v in self.elements.items()}
        for fg, fg_list in self.functional_groups.items():
            for func_group in fg_list:
                func_group.compute_elements()
                for k, v in func_group.get_elements().items():
                    if k not in elements: elements[k] = 0
                    elements[k] += v * func_group.count
        return elements
    
    
    def compute_elements(self):
        for fg, fg_list in self.functional_groups.items():
            for func_group in fg_list:
                func_group.compute_elements()
        
        
    def to_string(self, level):
        fg_string = ""
        if level == LipidLevel.ISOMERIC_SUBSPECIES:
            if str.isnumeric(self.name[0]): fg_string = "%i(%s)" % (self.position, self.name) if self.position > -1 else self.name
            else: fg_string = "%i%s" % (self.position, self.name) if self.position > -1 else self.name
            
        else:
            fg_string = "(%s)%i" % (self.name, self.count) if self.count > 1 else self.name
        if self.stereochemistry != None and self.stereochemistry != "" and level == LipidLevel.ISOMERIC_SUBSPECIES: fg_string += "[%s]" % self.stereochemistry
                
        return fg_string
    
    
    def get_num_oxygens(self):
        num_oxygens = self.elements[Element.O] * self.count if Element.O in self.elements else 0
        for fg, fg_list in self.functional_groups.items():
            for func_group in fg_list:
                func_group.compute_elements()
                if Element.O in func_group.elements:
                    num_oxygens += func_group.elements[Element.O] * func_group.count
        return num_oxygens
    
    
    
    def __iadd__(self, fg):
        for k, v in fg.elements.items():
            if k not in self.elements: self.elements[k] = 0
            self.elements[k] += v * fg.count
        
        for fg, fg_list in fg.functional_groups.items():
            for func_group in fg_list:
                func_group.compute_elements()
                for k, v in func_group.get_elements().items():
                    if k not in self.elements: self.elements[k] = 0
                    self.elements[k] += v * func_group.count
        
        return self
    
    
    
_known_functional_groups = {"Et": FunctionalGroup("Et", elements = {Element.C: 2, Element.H: 5}),
                           "Me": FunctionalGroup("Me", elements = {Element.C: 1, Element.H: 3}),
                           "Br": FunctionalGroup("Br", elements = {Element.Br: 1}),
                           "Cl": FunctionalGroup("Cl", elements = {Element.Cl: 1}),
                           "F": FunctionalGroup("F", elements = {Element.F: 1}),
                           "I": FunctionalGroup("I", elements = {Element.I: 1}),
                           "NO2": FunctionalGroup("NO2", elements = {Element.N: 1, Element.O: 2}),
                           "Ep": FunctionalGroup("Ep", elements = {Element.O: 1}),
                           "OO": FunctionalGroup("OO", elements = {Element.O: 2, Element.H: 1}),
                           "OMe": FunctionalGroup("OMe", elements = {Element.O: 1, Element.C: 1, Element.H: 3}),
                           "oxy": FunctionalGroup("oxy", elements = {}),
                           "NH2": FunctionalGroup("NH2", elements = {Element.N: 1, Element.H: 2}),
                           "OOH": FunctionalGroup("OOH", elements = {Element.O: 2, Element.H: 1}),
                           "SH": FunctionalGroup("SH", elements = {Element.S: 1, Element.H: 1}),
                           "OH": FunctionalGroup("OH", elements = {Element.O: 1, Element.H: 1}),
                           "oxo": FunctionalGroup("oxo", elements = {Element.O: 1}),
                           "CN": FunctionalGroup("CN", elements = {Element.C: 1, Element.N: 1}),
                           "P": FunctionalGroup("P", elements = {Element.P: 1}),
                           "S": FunctionalGroup("S", elements = {Element.S: 1}),
                           "COOH": FunctionalGroup("COOH", elements = {Element.C: 1, Element.O: 2, Element.H: 1}),
                           "G": FunctionalGroup("G", elements = {}),
                           "T": FunctionalGroup("T", elements = {Element.S: 1, Element.O: 3, Element.H: 1}),
                           "H": FunctionalGroup("H", elements = {Element.H: 1}),
                           "COG": FunctionalGroup("COG", elements = {}),
                           "COT": FunctionalGroup("COT", elements = {}),
                           "Hex": FunctionalGroup("Hex", elements = {}),
                           "Gal": FunctionalGroup("Gal", elements = {}),
                           "Glc": FunctionalGroup("Glc", elements = {}),
                           "Man": FunctionalGroup("Man", elements = {}),
                           "Neu": FunctionalGroup("Neu", elements = {}),
                           "HexNAc": FunctionalGroup("NexNAc", elements = {}),
                           "GalNAc": FunctionalGroup("GalNAc", elements = {}),
                           "GlcNAc": FunctionalGroup("GlcNAc", elements = {}),
                           "OGlcNAc": FunctionalGroup("OGlcNAc", elements = {}),
                           "OGlc": FunctionalGroup("OGlc", elements = {}),
                           "NeuAc": FunctionalGroup("NeuAc", elements = {}),
                           "NeuGc": FunctionalGroup("NeuGc", elements = {}),
                           "Kdn": FunctionalGroup("Kdn", elements = {}),
                           "GlcA": FunctionalGroup("GlcA", elements = {}),
                           "Xyl": FunctionalGroup("Xyl", elements = {}),
                           "Fuc": FunctionalGroup("Fuc", elements = {}),
                           "NeuAc2": FunctionalGroup("NeuAc2", elements = {}),
                           "SHex": FunctionalGroup("SHex", elements = {}),
                           "S(3')Hex": FunctionalGroup("S(3')Hex", elements = {}),
                           "S(3′)Hex": FunctionalGroup("S(3')Hex", elements = {}),
                           "NAc": FunctionalGroup("NAc", elements = {}),
                           "Nac": FunctionalGroup("Nac", elements = {}),
                           "SGal": FunctionalGroup("SGal", elements = {}),
                           "S(3')Gal": FunctionalGroup("S(3')Gal", elements = {}),
                           "S(3′)Gal": FunctionalGroup("S(3')Gal", elements = {}),
                           "HexA": FunctionalGroup("HexA", elements = {}),
                           "O": FunctionalGroup("O", elements = {Element.O: 1})}



def get_functional_group(name):
    if name in _known_functional_groups:
        return _known_functional_groups[name].copy()
    raise Exception("Name '%s' not registered in functional group list" % name)


    
class HeadGroupDecorator(FunctionalGroup):
    def __init__(self, name, position = -1, count = 1, elements = None):
        super().__init__(name, position = position, count = count, elements = elements)
        
        
    def to_string(self, level):
        return "%s-" % self.name
    
    

    
class AcylAlkylGroup(FunctionalGroup):
    def __init__(self, fa, position = -1, count = 1, alkyl = False):
        super().__init__("O", position = position, count = count)
        if fa != None: self.funcional_groups["alkyl" if alkyl else "acyl"] = [fa]
        self.alkyl = alkyl

    def to_string(self, level):
        acyl_alkyl_string = []
        if level == LipidLevel.ISOMERIC_SUBSPECIES: acyl_alkyl_string.append("%i" % self.position)
        acyl_alkyl_string.append("O(")
        if self.alkyl: acyl_alkyl_string.append("FA ")
        fa = self.functional_groups["alkyl" if self.alkyl else "acyl"][0]
        acyl_alkyl_string.append(fa.to_string(level))
        acyl_alkyl_string.append(")")
        
        return "".join(acyl_alkyl_string)
        
        
        
