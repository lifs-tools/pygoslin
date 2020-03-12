from pygoslin.domain.LipidSpecies import LipidSpecies
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidExceptions import ConstraintViolationException
from pygoslin.domain.LipidSpeciesInfo import LipidSpeciesInfo
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.LipidClass import *

class LipidMolecularSubspecies(LipidSpecies):


    def __init__(self, head_group, fa = []):
        super().__init__(head_group)
        self.fa = {}
        self.fa_list = []
        self.special_cases = {class_string_to_class["PC"], class_string_to_class["PE"], class_string_to_class["LPC"], class_string_to_class["LPE"]}
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
    

    def build_lipid_subspecies_name(self, fa_separator):
        special_case = self.lipid_class in self.special_cases
        
        fa_string = " " + fa_separator.join(fatty_acid.to_string(special_case) for fatty_acid in self.fa_list)
        if fa_string == " ": fa_string = ""
        
        return (all_lipids[self.lipid_class][0] if not self.use_head_group else self.head_group) + fa_string
    
    
    
    def get_lipid_string(self, level = None):
        if level == None or level == LipidLevel.MOLECULAR_SUBSPECIES:
            return self.build_lipid_subspecies_name("_")
        
        elif level in (LipidLevel.CATEGORY, LipidLevel.CLASS, LipidLevel.SPECIES):
            return super().get_lipid_string(level)
        else:
            raise Exception("LipidMolecularSubspecies does not know how to create a lipid string for level %s" % level)
    
    
