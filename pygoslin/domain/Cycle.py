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

class Cycle(FunctionalGroup):
    def __init__(self, cycle, start = None, end = None, double_bonds = None, functional_groups = None):
        super().__init__("cy", functional_groups = functional_groups)
        self.count = 1
        self.cycle = cycle
        self.position = start
        self.start = start
        self.end = end
        self.double_bonds = double_bonds if double_bonds != None else 0
        self.elements[Element.H] = -2
        
        if type(start) != type(end):
            raise ConstraintViolationException("Cycle data start and end values not of same type!")
        
        if type(start) == int:
            if end - start + 1 != cycle:
                raise ConstraintViolationException("Cycle data start (%i) and end (%i) values do not correspond to count (%i)!" % (start, end, cycle))
                
        
    def clone(self, cyc):
        pass
    
    
    def get_double_bonds(self):
        return super().get_double_bonds() + 1
    
    
    def rearrange_functional_groups(self, parent):
        ## put everything back into parent
        if type(parent.double_bonds) != dict: parent.double_bonds = {}
        if type(self.double_bonds) == dict and len(self.double_bonds) > 0:
            for k, v in self.double_bonds.items(): parent.double_bonds[k] = v
        
        for fg, fg_list in self.functional_groups.items():
            if fg not in fgroup: fgroup[fg] = []
            fgroup[fg] += fg_list
            
            
        ## take back what's mine
        remove_list = set()
        for fg, fg_list in parent.functional_groups.items():
            remove_item = []
            
            for i, func_group in enumerate(fg_list):
                if self.start <= func_group.position <= self.end:
                    if fg not in self.functional_groups: self.functional_groups[fg] = []
                    self.functional_groups[fg].append(func_group)
                    remove_item.append(i)
                    
            for i in remove_item[::-1]: del parent.functional_groups[fg][i]
            if len(fg_list) == 0: remove_list.add(fg)
            
        for fg in remove_list: del parent.functional_groups[fg]
        
        
    
    def shift_positions(self, shift):
        super().shift_positions(shift)
        self.start += shift
        self.end += shift
        if type(self.double_bonds) == dict:
            self.double_bonds = {k + shift: v for k, v in self.double_bonds.items()}
        
    
    def compute_elements(self):
        self.elements[Element.H] = -2
        if self.double_bonds != None:
            self.elements[Element.H] -= 2 * (self.double_bonds if type(self.double_bonds) == int else len(self.double_bonds))
        
        
    def to_string(self, level):
        cycle_string = ["["]
        if self.start != None and level == LipidLevel.ISOMERIC_SUBSPECIES:
            cycle_string.append("%i-%i" % (self.start, self.end))
        cycle_string.append("cy%i" % self.cycle)    
            
        if self.double_bonds != None:
            if type(self.double_bonds) != int:
                cycle_string.append(":%i" % len(self.double_bonds))
                db_positions = ["%i%s" % (k, self.double_bonds[k]) for k in sorted(self.double_bonds.keys())]
                db_pos = "(%s)" % ",".join(db_positions) if len (self.double_bonds) > 0 else ""
                cycle_string.append(db_pos)
            else:
                cycle_string.append(":%i" % self.double_bonds)
        
        
        if level == LipidLevel.ISOMERIC_SUBSPECIES:
            for fg in sorted(self.functional_groups.keys(), key = lambda x: x.lower()):
                fg_list = self.functional_groups[fg]
                fg_summary = ",".join([func_group.to_string(level) for func_group in sorted(fg_list, key = lambda x: x.position)])
                if len(fg_summary) > 0: cycle_string.append(";%s" % fg_summary)
        
        elif level == LipidLevel.STRUCTURAL_SUBSPECIES:
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
        if self.stereochemistry != None: cycle_string.append("[%s]" % self.stereochemistry)
        
        return "".join(cycle_string)
        
        
