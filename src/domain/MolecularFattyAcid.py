from .FattyAcid import FattyAcid
from .LipidExceptions import ConstraintViolationException

class MolecularFattyAcid(FattyAcid):
    
    def __init__(self, name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position = -1):
        super().__init__(name, num_carbon, num_hydroxyl, lipid_FA_bond_type, lcb, position)
        if num_double_bonds < 0:
            raise ConstraintViolationException("MolecularFattyAcid must have at least 0 double bonds!")
            
        self.num_double_bonds = num_double_bonds
    