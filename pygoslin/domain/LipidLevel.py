from enum import Enum

class LipidLevel(Enum):
    # Undefined / non-inferable lipid level
    UNDEFINED = 0
    
    # Mediators, Glycerolipids, Glycerophospholipids, Sphingolipids, Steroids, Prenols
    CATEGORY = 1
    
    # Glyerophospholipids -> Glycerophosphoinositols (PI)
    CLASS = 2
    
    # Phosphatidylinositol (16:0) or PI(16:0)
    SPECIES = 3
    
    # Phosphatidylinositol (8:0-8:0) or PI(8:0-8:0)
    MOLECULAR_SUBSPECIES = 4

    # Phosphatidylinositol (8:0/8:0) or PI(8:0/8:0)
    STRUCTURAL_SUBSPECIES = 5
    
    """
    1,2-dioctanoyl-sn-glycero-3-phospho-1D-myo-inositol
    PE(P-18:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z))
    Phosphatidylethanolamine (P-18:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z))
    """
    ISOMERIC_SUBSPECIES = 6
