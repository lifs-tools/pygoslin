
class FattyAcid:

    def __init__(self, name, position, num_carbon, num_hydroxy, lipid_FA_bond_type, lcb):
        self.name = name
        self.position = position
        self.num_carbon = carbon
        self.num_hydroxy = num_hydroxy
        self.lipid_FA_bond_type = lipid_FA_bond_type
        self.lcb = lcb
        
        if num_carbon < 2:
            raise ConstraintViolationException("FattyAcid must have at least 2 carbons!")
        
        if position < -1:
            raise ConstraintViolationException("FattyAcid position must be greater or equal to -1 (undefined) or greater or equal to 0 (0 = first position)!")
        
        if nHydroxy < 0:
            raise ConstraintViolationException("FattyAcid must have at least 0 hydroxy groups!")
        