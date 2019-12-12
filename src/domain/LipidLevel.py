from enum import Enum

class LipidLevel(Enum):
    # Undefined / non-inferable lipid level
    UNDEFINED = auto()
    
    # Mediators, Glycerolipids, Glycerophospholipids, Sphingolipids, Steroids, Prenols
    CATEGORY = auto()
    
    # Glyerophospholipids -> Glycerophosphoinositols (PI)
    CLASS = auto()
    
    # Phosphatidylinositol (16:0) or PI(16:0)
    SPECIES = auto()
    
    # Phosphatidylinositol (8:0-8:0) or PI(8:0-8:0)
    MOLECULAR_SUBSPECIES = auto()

    # Phosphatidylinositol (8:0/8:0) or PI(8:0/8:0)
    STRUCTURAL_SUBSPECIES = auto()
    
    """
    1,2-dioctanoyl-sn-glycero-3-phospho-1D-myo-inositol
    PE(P-18:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z))
    Phosphatidylethanolamine (P-18:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z))
    """
    ISOMERIC_SUBSPECIES = auto()
}
