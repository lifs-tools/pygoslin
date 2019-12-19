from pygoslin.domain.LipidMolecularSubspecies import LipidMolecularSubspecies
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidSpeciesInfo import LipidSpeciesInfo
from pygoslin.domain.LipidLevel import LipidLevel

class LipidStructuralSubspecies(LipidMolecularSubspecies):


    def __init__(self, head_group, fa = []):
        super().__init__(head_group)
        num_carbon = 0
        num_hydroxyl = 0
        num_double_bonds = 0
        lipid_FA_bond_type = LipidFaBondType.ESTER
        for fas in fa:
            if fas.name in self.fa:
                raise ConstraintViolationException("FA names must be unique! FA with name %s was already added!" % fas.name)
            else:
                self.fa[fas.name] = fas
                self.fa_list.append(fas)
                num_carbon += fas.num_carbon
                num_hydroxyl += fas.num_hydroxyl
                num_double_bonds += fas.num_double_bonds
                if lipid_FA_bond_type == LipidFaBondType.ESTER and fas.lipid_FA_bond_type in (LipidFaBondType.ETHER_PLASMANYL, LipidFaBondType.ETHER_PLASMENYL):
                    lipid_FA_bond_type = fas.lipid_FA_bond_type
                    
                elif lipid_FA_bond_type != LipidFaBondType.ESTER and fas.lipid_FA_bond_type in (LipidFaBondType.ETHER_PLASMANYL, LipidFaBondType.ETHER_PLASMENYL):
                    raise ConstraintViolationException("Only one FA can define an ether bond to the head group! Tried to add %s over existing %s" % (fas.lipid_FA_bond_type, lipid_FA_bond_type))
                
        self.info = LipidSpeciesInfo()
        self.info.level = LipidLevel.STRUCTURAL_SUBSPECIES
        self.info.num_carbon = num_carbon
        self.info.num_hydroxyl = num_hydroxyl
        self.info.num_double_bonds = num_double_bonds
        self.info.lipid_FA_bond_type = lipid_FA_bond_type
    

    
    def get_lipid_string(self, level = None):
        if level == None or level == LipidLevel.STRUCTURAL_SUBSPECIES:
            return super().build_lipid_subspecies_name("/")
        
        elif level in (LipidLevel.MOLECULAR_SUBSPECIES, LipidLevel.CATEGORY, LipidLevel.CLASS, LipidLevel.SPECIES):
            return super().get_lipid_string(level)
            
        else:
            raise RuntimeException("LipidStructuralSubspecies does not know how to create a lipid string for level %s")