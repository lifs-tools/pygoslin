from .LipidMolecularSubspecies import LipidMolecularSubspecies

class LipidStructuralSubspecies(LipidMolecularSubspecies):


    def __init__(self, head_group, fa):
        super().__init__(headGroup)
        num_carbon = 0
        num_hydroxyl = 0
        num_double_bonds = 0
        lipid_FA_bond_type = LipidFaBondType.ESTER
        for fas in fa:
            if fas.name in self.fa:
                raise ConstraintViolationException("FA names must be unique! FA with name %s was already added!" % fas.name)
            else:
                self.fa[fas.name] = fas
                num_carbon += fas.num_carbon
                num_hydroxyl += fas.num_hydroxyl
                num_double_bonds += fas.getNDoubleBonds();
                if lipid_FA_bond_type == LipidFaBondType.ESTER and fas.lipid_FA_bond_type in (LipidFaBondType.ETHER_PLASMANYL, LipidFaBondType.ETHER_PLASMENYL):
                    lipid_FA_bond_type = fas.lipid_FA_bond_type
                    
                elif lipid_FA_bond_type != LipidFaBondType.ESTER and fas.lipid_FA_bond_type in (LipidFaBondType.ETHER_PLASMANYL, LipidFaBondType.ETHER_PLASMENYL):
                    raise ConstraintViolationException("Only one FA can define an ether bond to the head group! Tried to add %s over existing %s" % (fas.lipid_FA_bond_type, lipid_FA_bond_type))
                
        super.info = LipidSpeciesInfo(LipidLevel.STRUCTURAL_SUBSPECIES, num_carbon, num_hydroxyl, num_double_bonds, lipid_FA_bond_type)
        self.lipid_species_string = super().build_lipid_subspecies_name("/")
    

    
    def get_lipid_string(self, level):
        if level == LipidLevel.STRUCTURAL_SUBSPECIES:
            return self.lipid_species_string
        elif level in (LipidLevel.MOLECULAR_SUBSPECIES, LipidLevel.CATEGORY, LipidLevel.CLASS, LipidLevel.SPECIES):
            return super().get_lipid_string(level)
            
        else:
            raise RuntimeException("LipidStructuralSubspecies does not know how to create a lipid string for level %s")