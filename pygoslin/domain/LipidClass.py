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
from pygoslin.domain.LipidCategory import *
from pygoslin.parser.ParserCommon import SumFormulaParser
from pygoslin.domain.Element import Element
from enum import Enum
from os import path
import pygoslin
import csv

UNDEFINED_LIPID_CLASS = 0

def get_category(name):
    global class_string_to_category
    return class_string_to_category[name] if name in class_string_to_category else LipidCategory.UNDEFINED
    
    
    
def get_class(name):
    global class_string_to_class
    return class_string_to_class[name] if name in class_string_to_class else UNDEFINED_LIPID_CLASS
    
    
    
    
class_string_to_class = {}
class_string_to_category = {}
all_lipids = [{"name": "UNDEFINED",
                "category": LipidCategory.UNDEFINED,
                "description": "Undefined lipid class",
                "max_fa": 0,
                "poss_fa": set(),
                "synonyms": []}]


all_lipids_dir_name = path.dirname(pygoslin.__file__)
with open(all_lipids_dir_name + "/data/goslin/lipid-list.csv", mode = "rt") as infile:
    line = infile.readline() # skip title row
    lipid_reader = csv.reader(infile, delimiter=',', quotechar='"')
    sum_formula_parser = SumFormulaParser()
    
    i = 1
    for row in lipid_reader:
        if len(row) > 0 and len(row[0]) > 0:
            
            lipid_dict = {}
            
            lipid_dict["name"] = row[0]
            lipid_dict["category"] = category_string_to_category[row[1]]
            lipid_dict["description"] = row[2]
            lipid_dict["max_fa"] = int(row[3])
            lipid_dict["poss_fa"] = set(int(pos) for pos in row[4].split("|"))
            row[5] = row[5].strip(" ")
            lipid_dict["elements"] = sum_formula_parser.parse(row[5]) if len(row[5]) > 0 else {e: 0 for e in Element}
            lipid_dict["synonyms"] = [r for r in row[6:] if len(r) > 0]
            
            
            all_lipids.append(lipid_dict)
            class_string_to_class[lipid_dict["name"]] = i
            class_string_to_category[lipid_dict["name"]] = lipid_dict["category"]
            for class_name in lipid_dict["synonyms"]:
                if len(class_name) > 0:
                    class_string_to_class[class_name] = i
                    class_string_to_category[class_name] = lipid_dict["category"]
            i += 1
        

        
