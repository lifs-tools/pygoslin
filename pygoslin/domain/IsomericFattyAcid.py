from pygoslin.domain.StructuralFattyAcid import StructuralFattyAcid

class IsomericFattyAcid(StructuralFattyAcid):

    def __init__(self, name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position, double_bond_positions):
        super().__init__(name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position)
        self.double_bond_positions = {key: double_bond_positions[key] for key in double_bond_positions}

    def clone(self, fa):
        self.double_bond_positions = {key: double_bond_positions[key] for key in fa.double_bond_positions}
        super().clone(fa)