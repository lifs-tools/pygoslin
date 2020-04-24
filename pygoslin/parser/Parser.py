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



import pygoslin
from os import path
from pygoslin.parser.GoslinParserEventHandler import GoslinParserEventHandler
from pygoslin.parser.GoslinFragmentParserEventHandler import GoslinFragmentParserEventHandler
from pygoslin.parser.LipidMapsParserEventHandler import LipidMapsParserEventHandler
from pygoslin.parser.SwissLipidsParserEventHandler import SwissLipidsParserEventHandler
from pygoslin.parser.HmdbParserEventHandler import HmdbParserEventHandler
from pygoslin.parser.ParserCommon import Parser
           

        
    
    
class GoslinParser(Parser):
    def __init__(self):
        self.event_handler = GoslinParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        super().__init__(self.event_handler, dir_name + "/data/goslin/Goslin.g4", Parser.DEFAULT_QUOTE)
        
        
        
class GoslinFragmentParser(Parser):
    def __init__(self):
        self.event_handler = GoslinFragmentParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        super().__init__(self.event_handler, dir_name + "/data/goslin/GoslinFragments.g4", Parser.DEFAULT_QUOTE)
        
        
class LipidMapsParser(Parser):
    def __init__(self):
        self.event_handler = LipidMapsParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        super().__init__(self.event_handler, dir_name + "/data/goslin/LipidMaps.g4", Parser.DEFAULT_QUOTE)
        
        
class SwissLipidsParser(Parser):
    def __init__(self):
        self.event_handler = SwissLipidsParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        super().__init__(self.event_handler, dir_name + "/data/goslin/SwissLipids.g4", Parser.DEFAULT_QUOTE)
        
        
        
        
class HmdbParser(Parser):
    def __init__(self):
        self.event_handler = HmdbParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        super().__init__(self.event_handler, dir_name + "/data/goslin/HMDB.g4", Parser.DEFAULT_QUOTE)
        
        
class LipidParser:
    def __init__(self):
        self.parser_list = [GoslinParser(), LipidMapsParser(), SwissLipidsParser(), HmdbParser(), GoslinFragmentParser()]
        
    def parse(self, lipid_name):
        for parser in self.parser_list:
            lipid = parser.parse(lipid_name, False)
            if lipid != None:
                return lipid
            
        raise LipidException("Lipid not found")

#init_tables()
