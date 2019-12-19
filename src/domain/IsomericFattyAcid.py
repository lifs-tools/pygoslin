from domain.StructuralFattyAcid import StructuralFattyAcid

class IsomericFattyAcid(StructuralFattyAcid):

    def __init__(self, name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, double_bond_positions, position):
        super().__init__(self, name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position)
        self.double_bond_positions = {key: double_bond_positions[key] for key in double_bond_positions}
