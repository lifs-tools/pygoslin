from pygoslin.domain.LipidLevel import LipidLevel

class LipidAdduct:

    def __init__(self):
        self.lipid = None
        self.adduct = None
        self.fragment = None
        self.sum_formula = None
        
        
    def get_lipid_string(self, level = None):
        lipid_name = []
        
        if self.lipid != None: lipid_name.append(self.lipid.get_lipid_string(level))
        else: return ""
        
        if self.adduct != None: lipid_name.append(self.adduct.get_lipid_string())
        
        if self.fragment != None: lipid_name.append(self.fragment.get_lipid_string())
        
        return "".join(lipid_name)