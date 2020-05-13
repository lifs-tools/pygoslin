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


from pygoslin.parser.BaseParserEventHandler import BaseParserEventHandler
from pygoslin.domain.Element import Element, element_positions
from pygoslin.domain.LipidExceptions import LipidException


class SumFormulaParserEventHandler(BaseParserEventHandler):
    
    def __init__(self):
        super().__init__()
        
        self.content = {e: 0 for e in Element}
        self.element = Element.H
        self.count = 0
        
        self.registered_events["molecule_pre_event"] = self.reset_parser
        self.registered_events["element_group_post_event"] = self.element_group_post_event
        self.registered_events["element_pre_event"] = self.element_pre_event
        self.registered_events["single_element_pre_event"] = self.single_element_group_pre_event
        self.registered_events["count_pre_event"] = self.count_pre_event
            
        
        
    def reset_parser(self, node):
        self.content = {e: 0 for e in Element}
    
    
    def element_group_post_event(self, node):
        self.content[self.element] += self.count
    
    
    def element_pre_event(self, node):
        
        parsed_element = node.get_text()
        if parsed_element in element_positions:
            self.element = element_positions[parsed_element]
            
        else:
            raise LipidException("Error: element '%s' is unknown" % parsed_element)
    
    
    
    
    def single_element_group_pre_event(self, node):
        
        parsed_element = node.get_text()
        if parsed_element in element_positions:
            self.element = element_positions[parsed_element]
            self.content[self.element] += 1
            
        else:
            raise LipidException("Error: element '%s' is unknown" % parsed_element)
    
    
    
    
    def count_pre_event(self, node):
        self.count = int(node.get_text())
