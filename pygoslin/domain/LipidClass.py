from pygoslin.domain.LipidCategory import *
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
all_lipids = [("UNDEFINED", LipidCategory.UNDEFINED, "Undefined lipid class")]

all_lipids_dir_name = path.dirname(pygoslin.__file__)
with open(all_lipids_dir_name + "/data/goslin/lipid-list.csv", mode = "rt") as infile:
    line = infile.readline() # skip title row
    lipid_reader = csv.reader(infile, delimiter=',', quotechar='"')
    i = 1
    for row in lipid_reader:
        if len(row) > 0 and len(row[0]) > 0:
            all_lipids.append(tuple(r for r in row if len(r) > 0))
            class_string_to_class[row[0]] = i
            class_string_to_category[row[0]] = category_string_to_category[row[1]]
            for class_name in row[3:]:
                if len(class_name) > 0:
                    class_string_to_class[class_name] = i
                    class_string_to_category[class_name] = category_string_to_category[row[1]]
            i += 1
        

        
