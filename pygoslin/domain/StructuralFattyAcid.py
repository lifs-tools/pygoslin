from pygoslin.domain.MolecularFattyAcid import MolecularFattyAcid

class StructuralFattyAcid(MolecularFattyAcid):
    def __init__(self, name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position):
        super().__init__(name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position)
        