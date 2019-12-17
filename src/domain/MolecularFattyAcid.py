from .FattyAcid import FattyAcid
from .LipidExceptions import ConstraintViolationException

class MolecularFattyAcid(FattyAcid):
        
    
    def __init__(self, name, num_carbon, num_hydroxy, num_double_bonds, lipid_FA_bond_type, lcb, position = -1):
        super().__init__(name, position, num_carbon, num_hydroxy, lipid_FA_bond_type, lcb)
        if num_double_bonds < 0:
            raise ConstraintViolationException("MolecularFattyAcid must have at least 0 double bonds!")
            
        self.num_double_bonds = num_double_bonds
    