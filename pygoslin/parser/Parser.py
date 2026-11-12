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
from os import path
from pygoslin.parser.GoslinParserEventHandler import GoslinParserEventHandler
from pygoslin.parser.LipidMapsParserEventHandler import LipidMapsParserEventHandler
from pygoslin.parser.SwissLipidsParserEventHandler import SwissLipidsParserEventHandler
from pygoslin.parser.HmdbParserEventHandler import HmdbParserEventHandler
from pygoslin.parser.ShorthandParserEventHandler import ShorthandParserEventHandler
from pygoslin.parser.FattyAcidParserEventHandler import FattyAcidParserEventHandler
from pygoslin.parser.LDAParserEventHandler import LDAParserEventHandler
#from pygoslin.parser.SystematicParserEventHandler import SystematicParserEventHandler

from pygoslin.parser.ParserCommon import Parser
from pygoslin.domain.LipidExceptions import LipidException, LipidParsingException
from pygoslin.domain.LipidClass import get_class, all_lipids, UNDEFINED_LIPID_CLASS
from pygoslin.domain.LipidCategory import LipidCategory
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.HeadGroup import HeadGroup
from pygoslin.domain.LipidSpecies import LipidSpecies
from pygoslin.domain.LipidAdduct import LipidAdduct

        
    
    
class GoslinParser(Parser):
    def __init__(self):
        self.event_handler = GoslinParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        file_name = path.join(dir_name, "data", "goslin", "Goslin.g4")
        super().__init__(self.event_handler, file_name, Parser.DEFAULT_QUOTE)
        
        

class LDAParser(Parser):
    def __init__(self):
        self.event_handler = LDAParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        file_name = path.join(dir_name, "data", "goslin", "LDA.g4")
        super().__init__(self.event_handler, file_name, Parser.DEFAULT_QUOTE)



class ShorthandParser(Parser):
    def __init__(self):
        self.event_handler = ShorthandParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        file_name = path.join(dir_name, "data", "goslin", "Shorthand2020.g4")
        super().__init__(self.event_handler, file_name, Parser.DEFAULT_QUOTE)
        
        
"""
class SystematicParser(Parser):
    def __init__(self):
        self.event_handler = SystematicParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        file_name = path.join(dir_name, "data", "goslin", "Systematic.g4")
        super().__init__(self.event_handler, file_name, Parser.DEFAULT_QUOTE)
        
    def parse(self, lipid_name):
        return super().parse(lipid_name.lower())
"""
        
        
class ClassLevelParser:
    """
    Parses names that consist of a lipid class or a lipid category only,
    e.g. 'PC', 'TG', 'HexCer', 'CE' (lipid class level) or 'GP', 'SP'
    (lipid category level). Class names and their synonyms come from
    lipid-list.csv, category names from LipidCategory. A category name that
    is also the primary name of a class (e.g. 'FA', 'SL') is read as that
    class; a category name that is only a class synonym (e.g. 'ST', a synonym
    of 'ST 27:1;O') is read as the category. Names with fatty acyl information
    are not handled here; they are parsed by the grammar based parsers.
    """
    CATEGORIES = {c.name: c for c in LipidCategory if c != LipidCategory.UNDEFINED}

    def __init__(self):
        self.parser_event_handler = self
        self.content = None
        self.error_message = ""
        self.word_in_grammar = False

    def parse(self, lipid_name, raise_error = True):
        self.content, self.error_message, self.word_in_grammar = None, "", False
        name = lipid_name.strip(" ")
        lipid_class = get_class(name)
        is_category = name in ClassLevelParser.CATEGORIES
        if lipid_class != UNDEFINED_LIPID_CLASS and all_lipids[lipid_class]["category"] != LipidCategory.UNDEFINED and (not is_category or all_lipids[lipid_class]["name"] == name):
            headgroup = HeadGroup(name)
            level = LipidLevel.CLASS
        elif is_category:
            headgroup = HeadGroup(name)
            headgroup.lipid_class = UNDEFINED_LIPID_CLASS
            headgroup.lipid_category = ClassLevelParser.CATEGORIES[name]
            level = LipidLevel.CATEGORY
        else:
            self.error_message = "Lipid '%s' is neither a lipid class nor a lipid category" % name
            if raise_error: raise LipidParsingException(self.error_message)
            return None

        lipid = LipidSpecies(headgroup)
        lipid.info.level = level
        self.content = LipidAdduct()
        self.content.lipid = lipid
        self.word_in_grammar = True
        return self.content



class FattyAcidParser(Parser):
    def __init__(self):
        self.event_handler = FattyAcidParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        file_name = path.join(dir_name, "data", "goslin", "FattyAcids.g4")
        super().__init__(self.event_handler, file_name, Parser.DEFAULT_QUOTE)
        
    def parse(self, lipid_name, raise_error = True):
        return super().parse(lipid_name.lower(), raise_error = raise_error)
        
        
class LipidMapsParser(Parser):
    def __init__(self):
        self.event_handler = LipidMapsParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        file_name = path.join(dir_name, "data", "goslin", "LipidMaps.g4")
        super().__init__(self.event_handler, file_name, Parser.DEFAULT_QUOTE)
        
        
class SwissLipidsParser(Parser):
    def __init__(self):
        self.event_handler = SwissLipidsParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        file_name = path.join(dir_name, "data", "goslin", "SwissLipids.g4")
        super().__init__(self.event_handler, file_name, Parser.DEFAULT_QUOTE)
        
        
        
        
class HmdbParser(Parser):
    def __init__(self):
        self.event_handler = HmdbParserEventHandler()
        dir_name = path.dirname(pygoslin.__file__)
        file_name = path.join(dir_name, "data", "goslin", "HMDB.g4")
        super().__init__(self.event_handler, file_name, Parser.DEFAULT_QUOTE)
        
        
        
class LipidParser:
    def __init__(self):
        self.parser_list = [ShorthandParser(), GoslinParser(), FattyAcidParser(), LipidMapsParser(), SwissLipidsParser(), HmdbParser(), LDAParser(), ClassLevelParser()]
        self.parser_event_handler = None
        """
        self.trivial_names = {}
        dir_name = path.dirname(pygoslin.__file__)
        file_name = path.join(dir_name, "data", "goslin", "trivial-names.csv")
        
        with open(file_name, "rt", encoding= "utf-8") as trivial_infile:
            trivial_infile.readline()
            for line in trivial_infile:
                line = line.strip()
                if len(line) < 1: continue
                row = Parser.split_string(line.strip(), ",", '"')
                if len(row) < 2: continue
            
                value = row[0].strip('"')
                for key in row[1:]:
                    key = key.strip('"').lower()
                    if len(key) < 1: continue
                
                    if key in self.trivial_names:
                        print("Error: key %s in trivial names list occurs twice." % key)
                        exit()
                    self.trivial_names[key] = value
        """
                
        
    def parse(self, lipid_name):
        # check if lipid name has a trivial name and translate it via predefined table
        #lipid_name_lower = lipid_name.lower()
        #if lipid_name_lower in self.trivial_names:
        #    lipid_name = self.trivial_names[lipid_name_lower]
        
        self.parser_event_handler = None
        # go through all parsers
        for parser in self.parser_list:
            lipid = parser.parse(lipid_name, raise_error = False)
            if lipid != None:
                self.parser_event_handler = parser.parser_event_handler
                return lipid
            
        raise LipidException("Lipid not found")

