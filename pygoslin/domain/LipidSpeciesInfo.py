
class LipidSpeciesInfo:
    
    def __init__(self, fa = None):
        self.level = None
        self.num_carbon = fa.num_carbon if fa != None else 0
        self.num_hydroxyl = fa.num_hydroxyl if fa != None else 0
        self.num_double_bonds = fa.num_double_bonds if fa != None else 0
        self.lipid_FA_bond_type = fa.lipid_FA_bond_type if fa != None else None