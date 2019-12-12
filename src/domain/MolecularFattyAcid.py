
class MolecularFattyAcid(FattyAcid):

    def __init__(self, name, num_carbon, num_hydroxy, num_double_bonds, lipid_FA_bond_type, lcb):
        this(name, -1, num_carbon, num_hydroxy, num_double_bonds, lipid_FA_bond_type, lcb)
        
    
    def __init__(self, name, position, num_carbon, num_hydroxy, num_double_bonds, lipid_FA_bond_type, lcb):
        super(name, position, num_carbon, nHydroxy, lipid_FA_bond_type, lcb)
        if num_double_bonds < 0:
            raise ConstraintViolationException("MolecularFattyAcid must have at least 0 double bonds!")
            
        self.num_double_bonds = num_double_bonds
