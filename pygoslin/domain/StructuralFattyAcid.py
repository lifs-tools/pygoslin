from pygoslin.domain.MolecularFattyAcid import MolecularFattyAcid

class StructuralFattyAcid(MolecularFattyAcid):
    def __init__(self, name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position):
        super().__init__(name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position)
        
    def clone(self, fa):
        super().clone(fa)
        
    def to_string(self, special_case = False):
        return super().to_string(special_case)