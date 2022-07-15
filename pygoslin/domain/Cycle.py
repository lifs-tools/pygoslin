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


from pygoslin.domain.FunctionalGroup import *
from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.Element import Element
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidLevel import LipidLevel

class Cycle(FunctionalGroup):
    def __init__(self, cycle, start = None, end = None, double_bonds = None, functional_groups = None, bridge_chain = None):
        super().__init__("cy", functional_groups = functional_groups)
        self.count = 1
        self.cycle = cycle
        self.position = start
        self.start = start
        self.end = end
        self.double_bonds = double_bonds if double_bonds != None else 0
        self.elements[Element.H] = -2
        self.bridge_chain = bridge_chain if bridge_chain != None else []
        
        if type(start) != type(end):
            raise ConstraintViolationException("Cycle data start and end values not of same type!")
        
            
                
        
    def copy(self):
        cy = Cycle(cycle = self.cycle, start = self.start, end = self.end, bridge_chain = list(self.bridge_chain))
        cy.double_bonds = {k: v for k, v in self.double_bonds.items()} if type(self.double_bonds) == dict else self.double_bonds
        for fg, fg_list in self.functional_groups.items():
            cy.functional_groups[fg] = [func_group.copy() for func_group in fg_list]
        return cy
                
    
    
    def get_double_bonds(self):
        return super().get_double_bonds() + 1
    
    
    
    def rearrange_functional_groups(self, parent, shift):
        ## put everything back into parent
        if type(parent.double_bonds) != dict: parent.double_bonds = {}
        if type(self.double_bonds) == dict and len(self.double_bonds) > 0:
            for k, v in self.double_bonds.items(): parent.double_bonds[k] = v
            self.double_bonds = {}
        
        fgroup = parent.functional_groups
        for fg, fg_list in self.functional_groups.items():
            if fg not in fgroup: fgroup[fg] = []
            fgroup[fg] += fg_list
        self.functional_groups = {}
            
        # shift the cycle
        self.shift_positions(shift)
            
        ## take back what's mine
        # check double bonds
        if type(parent.double_bonds) == dict and len(parent.double_bonds) > 0:
            self.double_bonds = {db_pos: val for db_pos, val in parent.double_bonds.items() if self.start <= db_pos <= self.end}
            
            for pos in self.double_bonds:
                del parent.double_bonds[pos]
        
        # check functional groups
        remove_list = set()
        for fg, fg_list in fgroup.items():
            remove_item = []
            
            for i, func_group in enumerate(fg_list):
                if self.start <= func_group.position <= self.end and func_group != self:
                    if fg not in self.functional_groups: self.functional_groups[fg] = []
                    self.functional_groups[fg].append(func_group)
                    remove_item.append(i)
                    
            for i in remove_item[::-1]: del fgroup[fg][i]
            if len(fg_list) == 0: remove_list.add(fg)
            
        for fg in remove_list: del fgroup[fg]
        
        
    
    def shift_positions(self, shift):
        super().shift_positions(shift)
        self.start += shift
        self.end += shift
        if type(self.double_bonds) == dict:
            self.double_bonds = {k + shift: v for k, v in self.double_bonds.items()}
        
    
    
    def compute_elements(self):
        self.elements = {e: 0 for e in element_order}
        self.elements[Element.H] = -2
        if self.double_bonds != None:
            self.elements[Element.H] -= 2 * (self.double_bonds if type(self.double_bonds) == int else len(self.double_bonds))
            
        for chain_element in self.bridge_chain:
            if chain_element == Element.C:
                self.elements[Element.C] += 1
                self.elements[Element.H] += 2
                
            if chain_element == Element.N:
                self.elements[Element.N] += 1
                self.elements[Element.H] += 1
                
            if chain_element == Element.As:
                self.elements[Element.As] += 1
                self.elements[Element.H] += 1
                
            if chain_element == Element.P:
                self.elements[Element.P] += 1
                self.elements[Element.H] += 1
                
            if chain_element == Element.O:
                self.elements[Element.O] += 1
                
            if chain_element == Element.S:
                self.elements[Element.S] += 1
            
        # add all implicit carbon chain elements
        if self.start != None and self.end != None:
            n = max(0, self.cycle - (self.end - self.start + 1 + len(self.bridge_chain)))
            self.elements[Element.C] += n
            self.elements[Element.H] += 2 * n
            
        
        
    def to_string(self, level):
        cycle_string = ["["]
        if self.start != None and level == LipidLevel.FULL_STRUCTURE:
            cycle_string.append("%i-%i" % (self.start, self.end))
        
        if level in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE, LipidLevel.STRUCTURE_DEFINED} and len(self.bridge_chain) > 0:
            cycle_string.append("".join(element_shortcut[e] for e in self.bridge_chain))
        cycle_string.append("cy%i" % self.cycle)    
          
        if level in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE, LipidLevel.STRUCTURE_DEFINED}:
            if self.double_bonds != None:
                if type(self.double_bonds) != int:
                    cycle_string.append(":%i" % len(self.double_bonds))
                    db_positions = ["%i%s" % (k, self.double_bonds[k]) for k in sorted(self.double_bonds.keys())] if level in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE} else ["%i" % k for k in sorted(self.double_bonds.keys())]
                    db_pos = "(%s)" % ",".join(db_positions) if len (self.double_bonds) > 0 else ""
                    cycle_string.append(db_pos)
                else:
                    cycle_string.append(":%i" % self.double_bonds)
            
        
        if level in {LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE}:
            for fg in sorted(self.functional_groups.keys(), key = lambda x: x.lower()):
                fg_list = self.functional_groups[fg]
                fg_summary = ",".join([func_group.to_string(level) for func_group in sorted(fg_list, key = lambda x: x.position)])
                if len(fg_summary) > 0: cycle_string.append(";%s" % fg_summary)
        
        elif level == LipidLevel.STRUCTURE_DEFINED:
            for fg in sorted(self.functional_groups.keys()):
                fg_list = self.functional_groups[fg]
                if len(fg_list) > 0:
                    if len(fg_list) == 1:
                        fg_summary = ",".join([func_group.to_string(level) for func_group in fg_list])
                        if len(fg_summary) > 0: cycle_string.append(";%s" % fg_summary)
                
                    else:
                        fg_count = sum([func_group.count for func_group in fg_list])
                        if fg_count > 1: cycle_string.append(";(%s)%i" % (fg, fg_count))
                        else: cycle_string.append(";%s" % fg)
                    
        cycle_string.append("]")
        if level == LipidLevel.COMPLETE_STRUCTURE and self.stereochemistry != None: cycle_string.append("[%s]" % self.stereochemistry)
        
        return "".join(cycle_string)
        
        
