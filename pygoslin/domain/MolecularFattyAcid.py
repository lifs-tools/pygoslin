from pygoslin.domain.FattyAcid import FattyAcid
from pygoslin.domain.LipidExceptions import ConstraintViolationException

class MolecularFattyAcid(FattyAcid):
    
    def __init__(self, name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position = -1):
        super().__init__(name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position)
        if num_double_bonds < 0:
            raise ConstraintViolationException("MolecularFattyAcid must have at least 0 double bonds!")
            
    
    def clone(self, fa):
        super().clone(fa)
        
    
    def to_string(self, special_case = False):
        return super().to_string(special_case)