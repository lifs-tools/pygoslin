from .LipidCategory import LipidCategory
from enum import Enum

class LipidClass(Enum):

    UNDEFINED = (LipidCategory.UNDEFINED, "UNDEFINED", "Undefined lipid class")
    
    ## Fatty acyls [FA] Fatty acids and conjugates [FA01]
    FA = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "FA")
    
    ## Glycerolipids [GL] Monoradylglycerols [GL01]
    MG = (LipidCategory.GL, "Monoradylglycerols [GL01]", "MG", "MAG")
    
    ## Diradylglycerols [GL02]
    DG = (LipidCategory.GL, "Diradylglycerols [GL02]", "DG", "DAG")
    
    ## Triradylglycerols [GL03]
    TG = (LipidCategory.GL, "Triradylglycerols [GL03]", "TG", "TAG")
    
    # TODO: there are some newer categories in LipidMaps, like Glycosylmono/di-radylglycerols, SQMG and SQDG */
    ## Glycerophospholipids [GP]
    BMP = (LipidCategory.GP, "Monoacylglycerophosphomonoradylglycerols [GP0410]", "BMP")
    CL = (LipidCategory.GP, "Glycerophosphoglycerophosphoglycerols [GP12]", "CL")
    MLCL = (LipidCategory.GP, "Glycerophosphoglycerophosphoglycerols [GP12]", "CL")
    PA = (LipidCategory.GP, "Glycerophosphates [GP10]", "PA")
    LPA = (LipidCategory.GP, "Glycerophosphates [GP10]", "LPA")
    PC = (LipidCategory.GP, "Glycerophosphocholines [GP01]", "PC")
    PC_O = (LipidCategory.GP, "Glycerophosphocholines [GP01]", "PC O")
    LPC = (LipidCategory.GP, "Glycerophosphocholines [GP01]", "LPC")
    LPC_O = (LipidCategory.GP, "Glycerophosphocholines [GP01]", "LPC O")
    PE = (LipidCategory.GP, "Glycerophosphoethanolamines [GP02]", "PE")
    PE_O = (LipidCategory.GP, "Glycerophosphoethanolamines [GP02]", "PE O")
    LPE = (LipidCategory.GP, "Glycerophosphoethanolamines [GP02]", "LPE")
    LPE_O = (LipidCategory.GP, "Glycerophosphoethanolamines [GP02]", "LPE O")
    PG = (LipidCategory.GP, "Glycerophosphoglycerols [GP04]", "PG")
    LPG = (LipidCategory.GP, "Glycerophosphoglycerols [GP04]", "LPG")
    PGP = (LipidCategory.GP, "Glycerophosphoglycerophosphates [GP05]", "PGP")
    PI = (LipidCategory.GP, "Glycerophosphoinositols [GP06]", "PI")
    LPI = (LipidCategory.GP, "Glycerophosphoinositols [GP06]", "LPI")
    PIP = (LipidCategory.GP, "Glycerophosphoinositol monophosphates [GP07]", "PIP")
    PIP_3p = (LipidCategory.GP, "Glycerophosphoinositol monophosphates [GP07]", "PIP[3']")
    PIP_4p = (LipidCategory.GP, "Glycerophosphoinositol monophosphates [GP07]", "PIP[4']")
    PIP_5p = (LipidCategory.GP, "Glycerophosphoinositol monophosphates [GP07]", "PIP[5']")
    PIP2 = (LipidCategory.GP, "Glycerophosphoinositol bisphosphates [GP08]", "PIP2")
    PIP2_3p_4p = (LipidCategory.GP, "Glycerophosphoinositol bisphosphates [GP08]", "PIP2[3',4']")
    PIP2_3p_5p = (LipidCategory.GP, "Glycerophosphoinositol bisphosphates [GP08]", "PIP2[3',5']")
    PIP3 = (LipidCategory.GP, "Glycerophosphoinositol trisphosphates [GP09]", "PIP3")
    PS = (LipidCategory.GP, "Glycerophosphoserines [GP03]", "PS")
    
    ## Sphingolipids
    CER = (LipidCategory.SP, "Ceramides [SP02]", "Cer")
    C1P = (LipidCategory.SP, "Ceramide-1-phosphates [SP0205]", "C1P")
    SPH = (LipidCategory.SP, "Sphingoid bases [SP01]", "SPH")
    S1P = (LipidCategory.SP, "Sphingoid bases [SP01]", "S1P")
    SM = (LipidCategory.SP, "Phosphosphingolipids [SP03]", "SM")
    HEXCER = (LipidCategory.SP, "Neutral glycosphingolipids [SP05]", "HexCer")
    GLCCER = (LipidCategory.SP, "Neutral glycosphingolipids [SP05]", "GlcCer")
    GALCER = (LipidCategory.SP, "Neutral glycosphingolipids [SP05]", "GalCer")
    HEX2CER = (LipidCategory.SP, "Neutral glycosphingolipids [SP05]", "Hex2Cer")
    HEX3CER = (LipidCategory.SP, "Neutral glycosphingolipids [SP05]", "Hex3Cer")
    LACCER = (LipidCategory.SP, "Neutral glycosphingolipids [SP05]", "LacCer")
    
    ## Sterol lipids
    ST = (LipidCategory.ST, "Sterols [ST01]", "ST")
    SE = (LipidCategory.ST, "Steryl esters [ST0102]", "SE")
    # FC = (LipidCategory.ST, "Cholesterol [LMST01010001]", "FC")
    CH = (LipidCategory.ST, "Cholesterol [LMST01010001]", "FC", "Ch", "Cholesterol")
    CHE = (LipidCategory.ST, "Cholesteryl esters [ST0102]", "ChE", "CE")
    


    def get_category(name):
        class_to_category = {key: lipid_class.value[0] for lipid_class in LipidClass for key in lipid_class.value[2:]}
        return class_to_category[name]
    
    def get_class(name):
        class_to_class = {key: lipid_class for lipid_class in LipidClass for key in lipid_class.value[2:]}
        return class_to_class[name]
        