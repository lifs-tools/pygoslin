from pygoslin.domain.LipidCategory import *
from enum import Enum
from os import path
import pygoslin
import csv

UNDEFINED_LIPID_CLASS = -1

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
            lipid_data = list(r for r in row if len(r) > 0)
            lipid_data[1] = category_string_to_category[lipid_data[1]]
            lipid_data[3] = int(lipid_data[3])
            all_lipids.append(tuple(lipid_data))
            class_string_to_class[lipid_data[0]] = i
            class_string_to_category[lipid_data[0]] = lipid_data[1]
            for class_name in lipid_data[4:]:
                if len(class_name) > 0:
                    class_string_to_class[class_name] = i
                    class_string_to_category[class_name] = lipid_data[1]
            i += 1
        

        
