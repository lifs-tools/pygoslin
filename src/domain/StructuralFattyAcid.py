from .MolecularFattyAcid import MolecularFattyAcid

class StructuralFattyAcid(MolecularFattyAcid):
    
    def __init__(self, name, position, num_carbon, num_hydroxy, num_double_bonds, lipid_FA_bond_type, lcb):
        super().__init__(name, position, num_carbon, num_hydroxy, num_double_bonds, lipid_FA_bond_type, lcb)
        