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
    def __init__(self, name, position = -1, count = 1, double_bonds = 0, stereochemistry = None, elements = None, functional_groups = None):
        self.name = name
        self.position = position
        self.count = count
        self.stereochemistry = stereochemistry
        self.ring_stereo = ""
        self.double_bonds = double_bonds
        self.elements = {e: 0 for e in Element} if elements == None else {k: v for k, v in elements.items()}
        self.functional_groups = functional_groups if functional_groups != None else {}
        
        
    def copy(self):
        functional_group = FunctionalGroup(self.name, position = self.position, count = self.count, double_bonds = self.double_bonds, stereochemistry = self.stereochemistry, elements = {k: v for k, v in self.elements.items()})
        for fg, fg_list in self.functional_groups.items():
            if fg not in functional_group.functional_groups: functional_group.functional_groups[fg] = []
            for func_group in fg_list:
                functional_group.functional_groups[fg].append(func_group)
        functional_group.ring_stereo = self.ring_stereo
        return functional_group
        
        
        
    def get_elements(self):
        self.compute_elements()
        elements = {e: 0 for e in Element}
        for k, v in self.elements.items(): elements[k] = v
        for k, v in self.get_functional_group_elements().items(): elements[k] += v
        return elements
    
    
    
    
    def get_functional_group_elements(self):
        elements = {e: 0 for e in Element}
        
        for fg, fg_list in self.functional_groups.items():
            for func_group in fg_list:
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
            if str.isnumeric(self.name[0]): fg_string = "%i%s(%s)" % (self.position, self.ring_stereo, self.name) if self.position > -1 else self.name
            else: fg_string = "%i%s%s" % (self.position, self.ring_stereo, self.name) if self.position > -1 else self.name
            
        else:
            fg_string = "(%s)%i" % (self.name, self.count) if self.count > 1 else self.name
        if self.stereochemistry != None and self.stereochemistry != "" and level == LipidLevel.ISOMERIC_SUBSPECIES: fg_string += "[%s]" % self.stereochemistry
                
        return fg_string
    
    
    def get_double_bonds(self):
        db = self.count * (self.double_bonds if type(self.double_bonds) == int else len(self.double_bonds))
        for fg, fg_list in self.functional_groups.items():
            for func_group in fg_list:
                db += func_group.get_double_bonds()
                
        return db
    
    
    def __iadd__(self, fgroup):
        
        for e in Element:
            if e not in self.elements: self.elements[e] = 0
        
        for k, v in fgroup.get_elements().items():
            self.elements[k] += v * fgroup.count
                 
        
        return self

    
class HeadgroupDecorator(FunctionalGroup):
    def __init__(self, name, position = -1, count = 1, elements = None, suffix = False, level = None):
        super().__init__(name, position = position, count = count, elements = elements)
        self.suffix = suffix
        self.lowest_visible_level = level
        
    def copy(self):
        return HeadgroupDecorator(self.name, position = self.position, count = self.count, elements = self.elements, suffix = self.suffix, level = self.lowest_visible_level)
        
        
        
        
    def to_string(self, level):
        if not self.suffix: return self.name
    
        decorator_string = ""
        if self.lowest_visible_level == None or self.lowest_visible_level.value <= level.value:
            
            if "decorator_alkyl" in self.functional_groups and len(self.functional_groups["decorator_alkyl"]) > 0:
                decorator_string = self.functional_groups["decorator_alkyl"][0].to_string(level) if level != LipidLevel.SPECIES else "Alk"
                
            elif "decorator_acyl" in self.functional_groups and len(self.functional_groups["decorator_acyl"]) > 0:
                decorator_string = "FA %s" % self.functional_groups["decorator_acyl"][0].to_string(level) if level != LipidLevel.SPECIES else "FA"
                
            else:
                decorator_string = self.name
                
            decorator_string = "(%s)" % decorator_string
            
        return decorator_string
    
    



    
class AcylAlkylGroup(FunctionalGroup):
    def __init__(self, fa, position = -1, count = 1, alkyl = False):
        
        super().__init__("O", position = position, count = count)
        self.alkyl = alkyl
        if fa != None: self.functional_groups["alkyl" if self.alkyl else "acyl"] = [fa]
        self.double_bonds = 1
        self.elements[Element.O] = 0 if self.alkyl else 1
        self.elements[Element.H] = -1
        
        

    def to_string(self, level):
        acyl_alkyl_string = []
        if level == LipidLevel.ISOMERIC_SUBSPECIES: acyl_alkyl_string.append("%i" % self.position)
        acyl_alkyl_string.append("O(")
        if not self.alkyl: acyl_alkyl_string.append("FA ")
        fa = self.functional_groups["alkyl" if self.alkyl else "acyl"][0]
        acyl_alkyl_string.append(fa.to_string(level))
        acyl_alkyl_string.append(")")
        
        return "".join(acyl_alkyl_string)
    
    
    
_known_functional_groups = {"OH": FunctionalGroup("OH", elements = {Element.O: 1}), # hydroxyl
                           "Me": FunctionalGroup("Me", elements = {Element.C: 1, Element.H: 2}), # methyl
                           "dMe": FunctionalGroup("dMe", elements = {Element.C: 1}), # methylen
                           "oxo": FunctionalGroup("oxo", elements = {Element.O: 1, Element.H: -2}, double_bonds = 1), # keto
                           "COOH": FunctionalGroup("COOH", elements = {Element.C: 1, Element.O: 2}), # carboxyl
                           "Ep": FunctionalGroup("Ep", elements = {Element.O: 1, Element.H: -2}), # epoxy
                           "OO": FunctionalGroup("OO", elements = {Element.O: 2}),  # peroxy
                           "OMe": FunctionalGroup("OMe", elements = {Element.O: 1, Element.C: 1, Element.H: 2}), # methoxy
                           "oxy": FunctionalGroup("oxy", elements = {Element.O: 1}), # Alkoxy / ether
                           "Et": FunctionalGroup("Et", elements = {Element.C: 2, Element.H: 4}), # ethyl
                           "Cl": FunctionalGroup("Cl", elements = {Element.Cl: 1, Element.H: -1}),
                           "F": FunctionalGroup("F", elements = {Element.F: 1, Element.H: -1}),
                           "Br": FunctionalGroup("Br", elements = {Element.Br: 1, Element.H: -1}),
                           "I": FunctionalGroup("I", elements = {Element.I: 1, Element.H: -1}),
                           "NH2": FunctionalGroup("NH2", elements = {Element.N: 1, Element.H: 1}),
                           "NO2": FunctionalGroup("NO2", elements = {Element.N: 1, Element.O: 2, Element.H: -1}),
                           "OOH": FunctionalGroup("OOH", elements = {Element.O: 2}),
                           "SH": FunctionalGroup("SH", elements = {Element.S: 1}),
                           "CN": FunctionalGroup("CN", elements = {Element.C: 1, Element.N: 1, Element.H: -1}),
                           "P": FunctionalGroup("P", elements = {Element.P: 1, Element.O: 4, Element.H: 1}),
                           "S": FunctionalGroup("S", elements = {Element.S: 1, Element.O: 4}),
                           "T": FunctionalGroup("T", elements = {Element.S: 1, Element.O: 3, Element.H: 1}),
                           "G": FunctionalGroup("G", elements = {Element.N: 1, Element.H: 1}),
                           "Hex": HeadgroupDecorator("Hex", elements = {Element.C: 6, Element.H: 10, Element.O: 5}),
                           "Gal": HeadgroupDecorator("Gal", elements = {Element.C: 6, Element.H: 10, Element.O: 5}),
                           "Glc": HeadgroupDecorator("Glc", elements = {Element.C: 6, Element.H: 10, Element.O: 5}),
                           "NeuAc": HeadgroupDecorator("NeuAc", elements = {Element.O: 8, Element.N: 1, Element.C: 11, Element.H: 17}),
                           "SGal": HeadgroupDecorator("SGal", elements = {Element.O: 8, Element.H: 10, Element.C: 6, Element.S: 1}),
                           "S(3')Gal": HeadgroupDecorator("S(3')Gal", elements = {Element.O: 8, Element.H: 10, Element.C: 6, Element.S: 1}),
                           "S(3′)Gal": HeadgroupDecorator("S(3')Gal", elements = {Element.O: 8, Element.H: 10, Element.C: 6, Element.S: 1}),
                           "SHex": HeadgroupDecorator("SHex", elements = {Element.O: 8, Element.H: 10, Element.C: 6, Element.S: 1}),
                           "S(3')Hex": HeadgroupDecorator("S(3')Hex", elements = {Element.O: 8, Element.H: 10, Element.C: 6, Element.S: 1}),
                           "S(3′)Hex": HeadgroupDecorator("S(3')Hex", elements = {Element.O: 8, Element.H: 10, Element.C: 6, Element.S: 1}),
                           "GlcA": HeadgroupDecorator("GlcA", elements = {Element.C: 6, Element.O: 6, Element.H: 8}),
                           "HexA": HeadgroupDecorator("HexA", elements = {Element.C: 6, Element.O: 6, Element.H: 8}),
                           "HexNAc": HeadgroupDecorator("NexNAc", elements = {Element.C: 8, Element.H: 13, Element.O: 5, Element.N: 1}),
                           "GalNAc": HeadgroupDecorator("GalNAc", elements = {Element.C: 8, Element.H: 13, Element.O: 5, Element.N: 1}),
                           "GlcNAc": HeadgroupDecorator("GlcNAc", elements = {Element.C: 8, Element.H: 13, Element.O: 5, Element.N: 1}),
                           "Man": HeadgroupDecorator("Man", elements = {Element.C: 6, Element.H: 9, Element.O: 5}),
                           "Neu": HeadgroupDecorator("Neu", elements = {Element.C: 9, Element.H: 14, Element.O: 7, Element.N: 1}),
                           "NeuGc": HeadgroupDecorator("NeuGc", elements = {Element.C: 11, Element.H: 17, Element.N: 1, Element.O: 9}),
                           "NAc": HeadgroupDecorator("NAc", elements = {Element.C: 5, Element.H: 7, Element.N: 1, Element.O: 2, Element.S: 1}),
                           "Nac": HeadgroupDecorator("Nac", elements = {Element.C: 5, Element.H: 7, Element.N: 1, Element.O: 2, Element.S: 1}),
                           "Fuc": HeadgroupDecorator("Fuc", elements = {Element.C: 6, Element.H: 10, Element.O: 4}),
                           "Kdn": HeadgroupDecorator("Kdn", elements = {Element.C: 9, Element.H: 14, Element.O: 8}),
                           "Xyl": HeadgroupDecorator("Xyl", elements = {Element.C: 29, Element.H: 52, Element.O: 26}),
                           "COG": HeadgroupDecorator("COG", elements = {Element.C: 18, Element.H: 19, Element.N: 5, Element.O: 1}),
                           "COT": HeadgroupDecorator("COT", elements = {Element.C: 6, Element.H: 14, Element.N: 2, Element.O: 2}),
                           
                           "H": FunctionalGroup("H", elements = {Element.H: 1}),
                           "OGlcNAc": HeadgroupDecorator("OGlcNAc", elements = {}),
                           "OGlc": HeadgroupDecorator("OGlc", elements = {Element.C: 6, Element.H: 10, Element.O: 5}),
                           "NeuAc2": HeadgroupDecorator("NeuAc2", elements = {}),
                           "O": FunctionalGroup("O", elements = {Element.O: 1})
                           }



def get_functional_group(name):
    if name in _known_functional_groups:
        return _known_functional_groups[name].copy()
    raise Exception("Name '%s' not registered in functional group list" % name)



        
        
        
