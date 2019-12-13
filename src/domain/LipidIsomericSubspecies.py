from .LipidStructuralSubspecies import LipidStructuralSubspecies

class LipidIsomericSubspecies(LipidStructuralSubspecies):


    def __init__(self, head_group, fa):
        LipidStructuralSubspecies.__init__(head_group)
        num_carbon = 0
        num_hydroxyl = 0
        num_double_bonds = 0
        lipid_FA_bond_type = LipidFaBondType.UNDEFINED
        if fa.length > 0:
            lipid_FA_bond_type = LipidFaBondType.ESTER
        
        for fas in fa:
            if fas.name in self.fa:
                raise ConstraintViolationException("FA names must be unique! FA with name %s was already added!" % fas.name)
            
            else:
                self.fa[fas.name] = fas
                num_carbon += fas.num_carbon
                num_hydroxyl += fas.num_hydroxyl
                num_double_bonds += fas.num_double_bonds
                
                if lipid_FA_bond_type and LipidFaBondType.ESTER and fas.lipid_FA_bond_type in (LipidFaBondType.ETHER_PLASMANYL, LipidFaBondType.ETHER_PLASMENYL):
                    lipid_FA_bond_type = fas.lipid_FA_bond_type
                    
                elif lipid_FA_bond_type != LipidFaBondType.ESTER and fas.lipid_FA_bond_type in (LipidFaBondType.ETHER_PLASMANYL, LipidFaBondType.ETHER_PLASMENYL):
                    raise ConstraintViolationException("Only one FA can define an ether bond to the head group! Tried to add %s over existing %s" % (fas.lipid_FA_bond_type, lipid_FA_bond_type))
                
        super().info = LipidSpeciesInfo(LipidLevel.ISOMERIC_SUBSPECIES, num_carbon, num_hydroxyl, num_double_bonds, lipid_FA_bond_type)
        self.lipid_species_string = build_lipid_isomeric_substructure_name()
    

    def build_lipid_isomeric_substructure_name(self):
        fa_strings = []
        for fa_key in self.fa:
            fatty_acid = fa[fa_key]
            num_carbon = 0
            num_hydroxyl = 0
            num_double_bonds = fatty_acid.num_double_bonds
            
            db_pos = ""
            db_positions = [key + fatty_acid.double_bond_positions[key] for key in fatty_acid.double_bond_positions]
            
            if len (fattyAcid.double_bond_positions) > 0:
                db_pos = "(%s)" % ",".join(db_positions)
                
            num_carbon += fattyAcid.num_carbon
            num_hydroxyl += fattyAcid.num_hydroxyl
            fa_strings.append("%i:%i%s%s%s" % (num_carbon, num_double_bonds, db_pos, ";" + str(num_hydroxy) if num_hydroxy > 0 else "", fatty_acid.lipid_FA_bond_type.suffix()))
            
        return self.head_group + " " + "/".join(self.fa_strings)
    

    def get_lipid_string(self, level):
        if level == LipidLevel.ISOMERIC_SUBSPECIES:
            return self.lipid_species_string
        
        elif level in (LipidLevel.STRUCTURAL_SUBSPECIES, LipidLevel.MOLECULAR_SUBSPECIES, LipidLevel.CATEGORY, LipidLevel.CLASS, LipidLevel.SPECIES):
            return super().get_lipid_string(level)
        else:
            raise Exception("LipidIsomericSubspecies does not know how to create a lipid string for level %s" % level)
