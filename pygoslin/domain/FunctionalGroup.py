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
    def __init__(self, name, position = -1, count = 1, elements = None):
        self.name = name
        self.position = position
        self.count = count
        self.elements = {e: 0 for e in Element} if elements == None else elements
        
        
    def clone(self, fg):
        self.name = fg.name
        self.elements = {k: v for k, v in fg.elements.items()}
        
        
    def get_elements(self):
        return self.elements
        
        
    def to_string(self, level):
        fg_string = ""
        if level == LipidLevel.ISOMERIC_SUBSPECIES:
            if str.isnumeric(self.name[0]): fg_string = "%i(%s)" % (self.position, self.name)
            else: fg_string = "%i%s" % (self.position, self.name)
            
        else:
            fg_string = "(%s)%i" % (self.name, self.count) if self.count > 1 else self.name
                
        return fg_string
    
    
    def get_num_oxygens(self):
        return self.elements[Element.O]
    
    
    def __iadd__(self, fg):
        for k, v in fg.elements.items():
            if k not in self.elements: self.elements[k] = 0
            self.elements[k] += v
        return self
