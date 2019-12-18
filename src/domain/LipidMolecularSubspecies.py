from .LipidSpecies import LipidSpecies
from .LipidFaBondType import LipidFaBondType
from .LipidExceptions import ConstraintViolationException
from .LipidSpeciesInfo import LipidSpeciesInfo
from .LipidLevel import LipidLevel

class LipidMolecularSubspecies(LipidSpecies):


    def __init__(self, head_group, fa = []):
        super().__init__(head_group)
        self.fa = {}
        self.fa_list = []
        num_carbon = 0
        num_hydroxyl = 0
        num_double_bonds = 0
        lipid_FA_bond_type = LipidFaBondType.ESTER;
        for fas in fa:
            if fas.position != -1:
                raise ConstraintViolationException("MolecularFattyAcid %s must have position set to -1! Was: %i"  % (fas.name, fas.position))
            
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
#                    num_double_bonds += lipid_FA_bond_type.doubleBondCorrection();
#                    log.debug("Correcting double bond count to {} due to ether bond.", num_double_bonds);
                
                elif lipid_FA_bond_type != LipidFaBondType.ESTER and fas.lipid_FA_bond_type in (LipidFaBondType.ETHER_PLASMANYL, LipidFaBondType.ETHER_PLASMENYL):
                    raise ConstraintViolationException("Only one FA can define an ether bond to the head group! Tried to add %s over existing %s" % (fas.lipid_FA_bond_type, lipid_FA_bond_type))
                
        self.info = LipidSpeciesInfo()
        self.info.level = LipidLevel.MOLECULAR_SUBSPECIES
        self.info.num_carbon = num_carbon
        self.info.num_hydroxyl = num_hydroxyl
        self.info.num_double_bonds = num_double_bonds
        self.info.lipid_FA_bond_type = lipid_FA_bond_type
        self.lipid_species_string = self.build_lipid_subspecies_name("_")
    

    def build_lipid_subspecies_name(self, fa_separator):
        fa_strings = []
        for fatty_acid in self.fa_list:
            num_double_bonds = fatty_acid.num_double_bonds
            num_carbon = fatty_acid.num_carbon
            num_hydroxyl = fatty_acid.num_hydroxyl
            fa_strings.append("%i:%i%s%s" % (num_carbon, num_double_bonds, ";" + str(num_hydroxyl) if num_hydroxyl > 0 else "", fatty_acid.lipid_FA_bond_type.suffix()))
            
        return self.head_group + " " + fa_separator.join(fa_strings)
    
    
    
    def get_lipid_string(self, level = None):
        
        if level == None or level == LipidLevel.MOLECULAR_SUBSPECIES:
            return self.lipid_species_string
        
        elif level in (LipidLevel.CATEGORY, LipidLevel.CLASS, LipidLevel.SPECIES):
            return super().get_lipid_string(level)
        else:
            raise Exception("LipidMolecularSubspecies does not know how to create a lipid string for level %s" % level)
    
    
