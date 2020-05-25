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


from pygoslin.parser.GoslinParserEventHandler import GoslinParserEventHandler
from pygoslin.domain.Fragment import Fragment

class GoslinFragmentParserEventHandler(GoslinParserEventHandler):
    
    def __init__(self):
        self.fragment = None
        
        super().__init__()
        self.reset_lipid(None)
    
        self.registered_events["lipid_post_event"] = self.build_lipid_fragment
        self.registered_events["fragment_name_pre_event"] = self.add_fragment
    
           
        
        
    def build_lipid_fragment(self, node):
        super().build_lipid(node)
        self.lipid.fragment = self.fragment
        

        
    def add_fragment(self, node):
        self.fragment = Fragment(node.get_text())
        
