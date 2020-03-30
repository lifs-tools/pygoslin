from pygoslin.domain.LipidExceptions import RuntimeException
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.LipidCategory import LipidCategory
from pygoslin.domain.LipidClass import *

class LipidSpecies:
    def __init__(self, head_group, lipid_category = None, lipid_class = None, lipid_species_info = None):
        self.head_group = head_group.strip(" ")
        self.lipid_category = lipid_category if lipid_category != None else get_category(self.head_group)
        
        self.lipid_class = lipid_class if lipid_class != None else get_class(self.head_group)
        self.info = lipid_species_info
        self.use_head_group = False
        
        
    def clone(self, fa):
        self.head_group = fa.head_group
        self.lipid_category = fa.lipid_category
        self.lipid_class = fa.lipid_class
        self.info = LipidSpeciesInfo(fa)
        self.use_head_group = fa.use_head_group
            


    def get_lipid_string(self, level = None):
        
        if level == None:
            if self.info != None:
                level = self.info.level
            else:
                raise RuntimeException("LipidSpecies does not know how to create a lipid string for level %s" % (level if level != None else " unknown"))
        
        if level == LipidLevel.CATEGORY:
            return self.lipid_category.name
        
        elif level == LipidLevel.CLASS:
            return all_lipids[self.lipid_class][0] if not self.use_head_group else self.head_group
        
        elif level == LipidLevel.SPECIES:
            lipid_string = [all_lipids[self.lipid_class][0] if not self.use_head_group else self.head_group]
            if self.info != None and self.info.num_carbon > 0:
                
                lipid_string += " " if all_lipids[self.lipid_class][1] != LipidCategory.ST else "/"
                lipid_string += self.info.to_string()
            return "".join(lipid_string)
        
        else:
            raise RuntimeException("LipidSpecies does not know how to create a lipid string for level %s" + level)
