from pygoslin.domain.StructuralFattyAcid import StructuralFattyAcid

class IsomericFattyAcid(StructuralFattyAcid):

    def __init__(self, name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position, double_bond_positions):
        super().__init__(name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position)
        self.double_bond_positions = {key: double_bond_positions[key] for key in double_bond_positions}

    def clone(self, fa):
        self.double_bond_positions = {key: double_bond_positions[key] for key in fa.double_bond_positions}
        super().clone(fa)
        
    
    def to_string(self, special_case = False):
        
        suffix = self.lipid_FA_bond_type.suffix()
        dbp = self.double_bond_positions
        db_positions = ["%i%s" % (k, dbp[k]) for k in sorted(dbp.keys())]
        db_pos = "(%s)" % ",".join(db_positions) if len (dbp) > 0 else ""
        
        return "%s%i:%i%s%s%s" % ("O-" if special_case and len(suffix) > 0 else "", self.num_carbon, self.num_double_bonds, db_pos, ";" + str(self.num_hydroxyl) if self.num_hydroxyl > 0 else "", suffix)