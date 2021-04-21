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
        super().__init__("cy", functional_groups)
        self.count = 1
        self.cycle = cycle
        self.start = start
        self.end = end
        self.double_bonds = double_bonds
        self.elements[Element.H] = -2
        
        if type(start) != type(end):
            raise ConstraintViolationException("Cycle data start and end values not of same type!")
        
        if type(start) == int:
            if end - start + 1 != cycle:
                raise ConstraintViolationException("Cycle data start (%i) and end (%i) values do not correspond to count (%i)!" % (start, end, cycle))
                
        
    def clone(self, cyc):
        pass
        
        
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
            for fg in sorted(self.functional_groups.keys()):
                fg_list = self.functional_groups[fg]
                fg_summary = ",".join([func_group.to_string(level) for func_group in fg_list])
                if len(fg_summary) > 0: cycle_string.append(";%s" % fg_summary)
        
        elif level == LipidLevel.STRUCTURAL_SUBSPECIES:
            for fg in sorted(self.functional_groups.keys()):
                fg_list = self.functional_groups[fg]
                if len(fg_list) > 0:
                    if len(fg_list) == 1:
                        fg_summary = ",".join([func_group.to_string(level) for func_group in fg_list])
                        if len(fg_summary) > 0: cycle_string.append(";%s" % fg_summary)
                
                    else:
                        fg_count = sum([func_group.cycle for func_group in fg_list])
                        if fg_count > 1: cycle_string.append(";(%s)%i" % (fg, fg_count))
                        else: cycle_string.append(";%s" % fg)
                    
        cycle_string.append("]")
        if self.stereochemistry != None: cycle_string.append("[%s]" % self.stereochemistry)
        
        return "".join(cycle_string)
        
        
