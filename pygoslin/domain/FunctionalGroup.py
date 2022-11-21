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

import pygoslin
from pygoslin.domain.Element import *
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.parser.ParserCommon import SumFormulaParser, Parser
from os import path

class FunctionalGroup:
    def __init__(self, name, position = -1, count = 1, double_bonds = 0, stereochemistry = None, elements = None, functional_groups = None, is_atomic = False, stereo_bound = False):
        self.name = name
        self.position = position
        self.count = count
        self.stereochemistry = stereochemistry
        self.stereo_bound = stereo_bound
        self.ring_stereo = ""
        self.double_bonds = double_bonds
        self.is_atomic = is_atomic
        self.elements = {e: 0 for e in Element} if elements == None else {k: v for k, v in elements.items()}
        self.functional_groups = functional_groups if functional_groups != None else {}
        
        
    def copy(self):
        functional_group = FunctionalGroup(self.name, position = self.position, count = self.count, double_bonds = self.double_bonds, stereochemistry = self.stereochemistry, elements = {k: v for k, v in self.elements.items()}, is_atomic = self.is_atomic, stereo_bound = self.stereo_bound)
        for fg, fg_list in self.functional_groups.items():
            if fg not in functional_group.functional_groups: functional_group.functional_groups[fg] = []
            for func_group in fg_list:
                functional_group.functional_groups[fg].append(func_group.copy())
        functional_group.ring_stereo = self.ring_stereo
        return functional_group
        
        
        
    def get_elements(self):
        self.compute_elements()
        elements = {e: 0 for e in Element}
        for k, v in self.elements.items(): elements[k] = v
        for k, v in self.get_functional_group_elements().items(): elements[k] += v
        return elements
    
    
    
    def shift_positions(self, shift):
        self.position += shift
        for fg_name, fg_list in self.functional_groups.items():
            for fg in fg_list:
                fg.shift_positions(shift)
                
                
                
    def stereo_information_missing(self):
        missing = self.stereo_bound and self.stereochemistry == ""
        for fg_name, fg_list in self.functional_groups.items():
            for fg in fg_list:
                missing |= fg.stereo_information_missing()
        return missing
    
    
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
        
        
        
    def add_functional_group(self, func_group):
        if func_group.name not in self.functional_groups: self.functional_groups[func_group.name] = []
        self.functional_groups[func_group.name].append(func_group)
        
        
    def to_string(self, level):
        fg_string = ""
        if level in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE}:
            if str.isnumeric(self.name[0]): fg_string = "%i%s(%s)" % (self.position, self.ring_stereo, self.name) if self.position > -1 else self.name
            else: fg_string = "%i%s%s" % (self.position, self.ring_stereo, self.name) if self.position > -1 else self.name
            
        else:
            fg_string = "(%s)%i" % (self.name, self.count) if self.count > 1 else self.name
        if self.stereochemistry != None and self.stereochemistry != "" and level == LipidLevel.COMPLETE_STRUCTURE: fg_string += "[%s]" % self.stereochemistry
                
        return fg_string
    
    
    
    def get_total_functional_group_count(self, fg_name):
        return sum([fg_item.count for fg_item in self.functional_groups[fg_name]]) if fg_name in self.functional_groups else 0
    
    
    
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
        if not self.suffix: return self.name + ("%i" % self.count if self.count > 1 else "")
    
        decorator_string = ""
        if self.lowest_visible_level == None or self.lowest_visible_level.value <= level.value:
            if "decorator_alkyl" in self.functional_groups:
                if len(self.functional_groups["decorator_alkyl"]) > 0:
                    decorator_string = self.functional_groups["decorator_alkyl"][0].to_string(level) if level.value > LipidLevel.SPECIES.value else "Alk"
                else:
                    decorator_string = "Alk"
                
            elif "decorator_acyl" in self.functional_groups:
                if len(self.functional_groups["decorator_acyl"]) > 0:
                    decorator_string = "FA %s" % self.functional_groups["decorator_acyl"][0].to_string(level) if level.value > LipidLevel.SPECIES.value else "FA"
                else:
                    decorator_string = "FA"
                
            else:
                decorator_string = self.name
                
            decorator_string = "(%s)" % decorator_string
            
        return decorator_string
    
    



    
class AcylAlkylGroup(FunctionalGroup):
    def __init__(self, fa, position = -1, count = 1, alkyl = False, N_bond = False):
        
        super().__init__("O", position = position, count = count)
        self.alkyl = alkyl
        if fa != None: self.functional_groups["alkyl" if self.alkyl else "acyl"] = [fa]
        self.double_bonds = int(not self.alkyl)
        self.set_N_bond_type(N_bond)
        
        if self.N_bond:
            self.elements[Element.H] = 2 if self.alkyl else 0
            self.elements[Element.O] = -1 if self.alkyl else 0
            self.elements[Element.N] = 1
            
        else:
            self.elements[Element.H] = 1 if self.alkyl else -1
            self.elements[Element.O] = 0 if self.alkyl else 1
        
        
    def copy(self):
        return AcylAlkylGroup(self.functional_groups["alkyl" if self.alkyl else "acyl"][0].copy(), alkyl = self.alkyl, position = self.position, count = self.count, N_bond = self.N_bond)
        
        
        
    def set_N_bond_type(self, N_bond):
        self.N_bond = N_bond
        
        if self.N_bond:
            self.elements[Element.H] = 2 if self.alkyl else 0
            self.elements[Element.O] = -1 if self.alkyl else 0
            self.elements[Element.N] = 1
            
        else:
            self.elements[Element.H] = 1 if self.alkyl else -1
            self.elements[Element.O] = 0 if self.alkyl else 1
        
        

    def to_string(self, level):
        acyl_alkyl_string = []
        if level in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE}: acyl_alkyl_string.append("%i" % self.position)
        acyl_alkyl_string.append("%s(" % ("N" if self.N_bond else "O"))
        if not self.alkyl: acyl_alkyl_string.append("FA ")
        fa = self.functional_groups["alkyl" if self.alkyl else "acyl"][0]
        acyl_alkyl_string.append(fa.to_string(level))
        acyl_alkyl_string.append(")")
        
        return "".join(acyl_alkyl_string)
    
    



    
class CarbonChain(FunctionalGroup):
    def __init__(self, fa, position = -1, count = 1):
        super().__init__("cc", position = position, count = count)
        if fa != None: self.functional_groups["cc"] = [fa]
        
        self.elements[Element.H] = 1
        self.elements[Element.O] = -1
        
        
    def copy(self):
        return CarbonChain(self.functional_groups["cc"][0].copy(), position = self.position, count = self.count)
        
        

    def to_string(self, level):
        return "%s(%s)" % ((str(self.position) if level in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE} else ""), self.functional_groups["cc"][0].to_string(level))
    
    
    
_known_functional_groups, s = {}, SumFormulaParser()
fg_dir_name = path.dirname(pygoslin.__file__)
fg_file_name = path.join(fg_dir_name, "data", "goslin", "functional-groups.csv")
with open(fg_file_name, mode = "rt", encoding= "utf-8") as fg_infile:
    fg_infile.readline()
    for line in fg_infile:
        line = line.strip()
        row = Parser.split_string(line.strip(), ",", '"', True)
        e = s.parse(row[2]) if len(row[2]) > 0 else {e: 0 for e in Element}
        row.append(row[1])
        for i in range(7, len(row)):
            key = row[i]
            if len(key) == 0: continue
        
            if row[0] == "FG":
                _known_functional_groups[key] = FunctionalGroup(row[1], elements = e, double_bonds = int(row[3]), is_atomic = row[4] == "1", stereo_bound = row[5] == "1")
                
            elif row[0] == "HGD":
                _known_functional_groups[key] = HeadgroupDecorator(row[1], elements = e)
        
    

def get_functional_group(name):
    if name in _known_functional_groups:
        return _known_functional_groups[name].copy()
    raise Exception("Name '%s' not registered in functional group list" % name)



        
        
        
