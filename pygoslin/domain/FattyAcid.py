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


from pygoslin.domain.FunctionalGroup import *
from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.Element import Element
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidLevel import LipidLevel

class FattyAcid(FunctionalGroup):

    def __init__(self, name, num_carbon = 0, double_bonds = 0, functional_groups = None, lipid_FA_bond_type = LipidFaBondType.ESTER, lcb = False, position = 0):
        super().__init__(name, double_bonds = double_bonds)
        self.position = position
        self.num_carbon = num_carbon
        self.lipid_FA_bond_type = lipid_FA_bond_type
        self.lcb = lcb
        
        num_double_bonds = len(self.double_bonds) if type(self.double_bonds) != int else self.double_bonds
        if num_carbon < 0 or num_carbon == 1:
            raise ConstraintViolationException("FattyAcid must have at least 2 carbons!")
        
        if num_double_bonds < 0:
            raise ConstraintViolationException("FattyAcid must have at least 0 double bonds!")
        
        if position < 0:
            raise ConstraintViolationException("FattyAcid must be at least 0 at position 0!")
        
        
    def copy(self):
        fa = FattyAcid(self.name)
        fa.position = self.position
        fa.num_carbon = self.num_carbon
        fa.lipid_FA_bond_type = self.lipid_FA_bond_type
        fa.lcb = self.lcb
        fa.double_bonds = {key: value for key, value in self.double_bonds.items()} if type(self.double_bonds) != int else self.double_bonds
        fa.functional_groups = {}
        for fg, fg_list in self.functional_groups.items():
            fa.functional_groups[fg] = []
            for fg_item in fg_list:
                func_group = FunctionalGroup("")
                func_group.clone(fg_item)
                fa.function_groups[fg].append(func_group.copy())
        return fa

        
    def clone(self, fa):
        self.name = fa.name
        self.position = fa.position
        self.num_carbon = fa.num_carbon
        self.lipid_FA_bond_type = fa.lipid_FA_bond_type
        self.lcb = fa.lcb
        self.double_bonds = {key: value for key, value in fa.double_bonds.items()} if type(fa.double_bonds) != int else fa.double_bonds
        self.functional_groups = {}
        for fg, fg_list in fa.functional_groups.items():
            self.functional_groups[fg] = []
            for fg_item in fg_list:
                func_group = FunctionalGroup("")
                func_group.clone(fg_item)
                self.function_groups[fg].append(func_group.copy())
                
                
                
    def get_double_bonds(self):
        return super().get_double_bonds() + (self.lipid_FA_bond_type == LipidFaBondType.ETHER_PLASMENYL)
        
                
                
        
    def to_string(self, level):
        fa_string = [self.lipid_FA_bond_type.prefix()]
        num_carbon = self.num_carbon
        double_bonds = self.double_bonds
        num_oxygen = 0
        
        if num_carbon == 0 and double_bonds == 0 and level not in {LipidLevel.ISOMERIC_SUBSPECIES, LipidLevel.STRUCTURAL_SUBSPECIES}:
            return ""
        
        if level == LipidLevel.MOLECULAR_SUBSPECIES:
            num_carbon = self.get_elements()[Element.C]
            double_bonds = self.get_double_bonds() - (self.lipid_FA_bond_type == LipidFaBondType.ETHER_PLASMENYL)
            num_oxygen = self.get_functional_group_elements()[Element.O]
            

        fa_string.append("%i" % num_carbon)
        
        if type(double_bonds) != int:
            fa_string.append(":%i" % len(double_bonds))
            if level == LipidLevel.ISOMERIC_SUBSPECIES:
                db_positions = ["%i%s" % (k, double_bonds[k]) for k in sorted(double_bonds.keys())]
                db_pos = "(%s)" % ",".join(db_positions) if len (double_bonds) > 0 else ""
                fa_string.append(db_pos)
            elif level == LipidLevel.STRUCTURAL_SUBSPECIES:
                db_positions = ["%i" % k for k in sorted(double_bonds.keys())]
                db_pos = "(%s)" % ",".join(db_positions) if len (double_bonds) > 0 else ""
                fa_string.append(db_pos)
            
        else:
            fa_string.append(":%i" % double_bonds)
        
        if level == LipidLevel.ISOMERIC_SUBSPECIES:
            for fg in sorted(self.functional_groups.keys(), key = lambda x: x.lower()):
                fg_list = self.functional_groups[fg]
                fg_summary = ",".join([func_group.to_string(level) for func_group in sorted(fg_list, key = lambda x: x.position)])
                if len(fg_summary) > 0: fa_string.append(";%s" % fg_summary)
        
        elif level == LipidLevel.STRUCTURAL_SUBSPECIES:
            for fg in sorted(self.functional_groups.keys()):
                fg_list = self.functional_groups[fg]
                if len(fg_list) > 0:
                    
                    if fg in {"acyl", "alkyl", "cy", "cc"}:
                        fg_summary = ",".join([func_group.to_string(level) for func_group in fg_list])
                        if len(fg_summary) > 0: fa_string.append(";%s" % fg_summary)
                    
                    else:
                        fg_count = sum([func_group.count for func_group in fg_list])
                        if fg_count > 1: fa_string.append(";(%s)%i" % (fg, fg_count) if not fg_list[0].is_atomic else ";%s%i" % (fg, fg_count))
                        else: fa_string.append(";%s" % fg)
        
        else:
            elements = self.get_functional_group_elements()
            additional_elements = ";".join(("%s%i" % (element_shortcut[e], elements[e]) if elements[e] > 1 else "%s" % element_shortcut[e]) for e in element_order[2:] if e in elements and elements[e] > 0)
            
            if len(additional_elements) > 0:
                fa_string.append(";")
                fa_string.append(additional_elements)
        
        return "".join(fa_string)
    
    


    def compute_elements(self):
        self.elements = {e: 0 for e in Element}
        num_double_bonds = len(self.double_bonds) if type(self.double_bonds) != int else self.double_bonds
        if self.lipid_FA_bond_type == LipidFaBondType.ETHER_PLASMENYL: num_double_bonds += 1
        
        if self.num_carbon == 0 and num_double_bonds == 0:
            self.elements[Element.H] = 1
            return
        
        if not self.lcb:
            
            self.elements[Element.C] = self.num_carbon # carbon
            if self.lipid_FA_bond_type == LipidFaBondType.ESTER:
                self.elements[Element.H] = (2 * self.num_carbon - 1 - 2 * num_double_bonds) # hydrogen
                self.elements[Element.O] = 1 # oxygen
            
            elif self.lipid_FA_bond_type == LipidFaBondType.ETHER_PLASMENYL:
                self.elements[Element.H] = (2 * self.num_carbon - 1 - 2 * num_double_bonds + 2) # hydrogen
            
            elif self.lipid_FA_bond_type == LipidFaBondType.ETHER_PLASMANYL:
                self.elements[Element.H] = ((self.num_carbon + 1) * 2 - 1 - 2 * num_double_bonds) # hydrogen
            
            elif self.lipid_FA_bond_type == LipidFaBondType.AMINE:
                self.elements[Element.H] = (2 * self.num_carbon + 1 - 2 * num_double_bonds) # hydrogen
                
            else:
                raise LipidException("Mass cannot be computed for fatty acyl chain with bond type: %s" % self.lipid_FA_bond_type)
                
        else:
            # long chain base
            self.elements[Element.C] = self.num_carbon # carbon
            self.elements[Element.H] = (2 * (self.num_carbon - num_double_bonds) + 1) # hydrogen
            self.elements[Element.N] = 1 # nitrogen
            self.elements[Element.O] = 1 # oxygen
            
            
