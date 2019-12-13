from .StructuralFattyAcid import StructuralFattyAcid

class IsomericFattyAcid(StructuralFattyAcid):
    
    def __init__(self, name, position, num_carbon, num_hydroxy, lipid_FA_bond_type, lcb, double_bond_positions):
        super().__init__(self, name, position, nCarbon, nHydroxy, len(double_bond_positions), lipidFaBondType, lcb)
        self.double_bond_positions = {key: double_bond_positions[key] for key in double_bond_positions}
