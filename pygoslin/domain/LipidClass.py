from pygoslin.domain.LipidCategory import LipidCategory
from enum import Enum

class LipidClass(Enum):

    UNDEFINED = (LipidCategory.UNDEFINED, "UNDEFINED", "Undefined lipid class")
    
    ## Fatty acyls [FA] Fatty acids and conjugates [FA01]
    FA1 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "FA")
    FA2 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "10-HDoHE")
    FA3 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "11-HDoHE")
    FA4 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "11-HETE")
    FA5 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "11,12-DHET")
    FA6 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "11(12)-EET")
    FA7 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]",  "12-HEPE")
    FA8 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "12-HETE")
    FA9 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "12-HHTrE")
    FA10 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "12-OxoETE")
    FA11 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "12(13)-EpOME")
    FA12 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "13-HODE")
    FA13 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "13-HOTrE")
    FA14 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "14,15-DHET")
    FA15 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "14(15)-EET")
    FA16 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "14(15)-EpETE")
    FA17 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "15-HEPE")
    FA18 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "15-HETE")
    FA19 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "15d-PGJ2")
    FA20 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "16-HDoHE")
    FA21 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "16-HETE")
    FA22 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "18-HEPE")
    FA23 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "5-HEPE")
    FA24 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "5-HETE")
    FA25 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "5-HpETE")
    FA26 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "5-OxoETE")
    FA27 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "5,12-DiHETE")
    FA28 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "5,6-DiHETE")
    FA29 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "5,6,15-LXA4")
    FA30 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "5(6)-EET")
    FA31 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "8-HDoHE")
    FA32 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "8-HETE")
    FA33 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "8,9-DHET")
    FA34 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "8(9)-EET")
    FA35 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "9-HEPE")
    FA36 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "9-HETE")
    FA37 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "9-HODE")
    FA38 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "9-HOTrE")
    FA39 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "9(10)-EpOME")
    FA40 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "alpha-LA")
    FA41 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "DHA")
    FA42 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "EPA")
    FA43 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "Linoleic acid")
    FA44 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "LTB4")
    FA45 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "LTC4")
    FA46 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "LTD4")
    FA47 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "Maresin 1")
    FA48 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "Palmitic acid")
    FA49 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "PGB2")
    FA50 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "PGD2")
    FA51 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "PGE2")
    FA52 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "PGF2alpha")
    FA53 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "PGI2")
    FA54 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "Resolvin D1")
    FA55 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "Resolvin D2")
    FA56 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "Resolvin D3")
    FA57 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "Resolvin D5")
    FA58 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "tetranor-12-HETE")
    FA59 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "TXB1")
    FA60 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "TXB2")
    FA61 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "TXB3")
    FA62 = (LipidCategory.FA, "Fatty acids and conjugates [FA01]", "AA", "Arachidonic acid", "Arachidonic acid")
    
    
    
    ## Glycerolipids [GL]
    MG = (LipidCategory.GL, "Monoradylglycerols [GL01]", "MAG", "MG") ## Glycerolipids [GL] Monoradylglycerols [GL01]
    DG = (LipidCategory.GL, "Diradylglycerols [GL02]", "DAG", "DG") ## Diradylglycerols [GL02]
    TG = (LipidCategory.GL, "Triradylglycerols [GL03]", "TAG", "TG") ## Triradylglycerols [GL03]
    MGDG = (LipidCategory.GL, "Glycosyldiradylglycerols [GL05]", "MGDG")
    DGDG = (LipidCategory.GL, "Glycosyldiradylglycerols [GL05]", "DGDG")
    SQMG = (LipidCategory.GL, "Glycosylmonoradylglycerols [GL04]", "SQMG")
    SQDG = (LipidCategory.GL, "Glycosyldiradylglycerols [GL05]", "SQDG")
    
    
    
    # TODO: there are some newer categories in LipidMaps, like Glycosylmono/di-radylglycerols, SQMG and SQDG */
    ## Glycerophospholipids [GP]
    BMP = (LipidCategory.GP, "Monoacylglycerophosphomonoradylglycerols [GP0410]", "BMP", "LBPA")
    CDPDAG = (LipidCategory.GP, "CDP-Glycerols [GP13]", "CDPDAG", "CDP-DG")
    CL = (LipidCategory.GP, "Glycerophosphoglycerophosphoglycerols [GP12]", "CL")
    MLCL = (LipidCategory.GP, "Glycerophosphoglycerophosphoglycerols [GP12]", "MLCL")
    DLCL = (LipidCategory.GP, "Glycerophosphoglycerophosphoglycerols [GP12]", "DLCL")
    PA = (LipidCategory.GP, "Glycerophosphates [GP10]", "PA")
    LPA = (LipidCategory.GP, "Glycerophosphates [GP10]", "LPA")
    PC = (LipidCategory.GP, "Glycerophosphocholines [GP01]", "PC")
    LPC = (LipidCategory.GP, "Glycerophosphocholines [GP01]", "LPC", "LysoPC")
    PE = (LipidCategory.GP, "Glycerophosphoethanolamines [GP02]", "PE")
    PET = (LipidCategory.GP, "Glycerophosphoethanolamines [GP02]", "PEt")
    LPE = (LipidCategory.GP, "Glycerophosphoethanolamines [GP02]", "LPE", "LysoPE")
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
    PIP2_4p_5p = (LipidCategory.GP, "Glycerophosphoinositol bisphosphates [GP08]", "PIP2[4',5']")
    PIP3 = (LipidCategory.GP, "Glycerophosphoinositol trisphosphates [GP09]", "PIP3")
    PIP3_3p_4p_5p = (LipidCategory.GP, "Glycerophosphoinositol trisphosphates [GP09]", "PIP3[3',4',5']")
    PS = (LipidCategory.GP, "Glycerophosphoserines [GP03]", "PS")
    LPS = (LipidCategory.GP, "Glycerophosphoserines [GP03]", "LPS")
    PIM1 = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "PIM1")
    PIM2 = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "PIM2")
    PIM3 = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "PIM3")
    PIM4 = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "PIM4")
    PIM5 = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "PIM5")
    PIM6 = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "PIM6")
    GLCDG = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "Glc-DG")
    PENME2 = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "PE-NMe2")
    AC2SGL = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "AC2SGL")
    PENME = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "PE-NMe")
    PT = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "PT")
    GLCGP = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "Glc-GP")
    NAPE = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "NAPE")
    LPIM1 = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "LPIM1")
    LPIM2 = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "LPIM2")
    LPIM3 = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "LPIM3")
    LPIM4 = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "LPIM4")
    LPIM5 = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "LPIM5")
    LPIM6 = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "LPIM6")
    CPA = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "CPA")
    SLBPA = (LipidCategory.GP, "Glycerophosphoglycerols [GP04]", "SLBPA")
    GL_6_AC_GlC_GP = (LipidCategory.GP, "Glycosylglycerophospholipids [GP14]", "6-Ac-Glc-GP")
    PNC = (LipidCategory.GP, "Glycerophosphonocholines [GP16]", "PnC")
    PNE = (LipidCategory.GP, "Glycerophosphoinositolglycans [GP15]", "PnE")
    
    ## Sphingolipids
    CER = (LipidCategory.SP, "Ceramides [SP02]", "Cer")
    CERP = (LipidCategory.SP, "Ceramides [SP02]", "CerP", "C1P")
    SM = (LipidCategory.SP, "Phosphosphingolipids [SP03]", "SM")
    HEXCER = (LipidCategory.SP, "Neutral glycosphingolipids [SP05]", "HexCer", "GalCer", "GlcCer")
    HEX2CER = (LipidCategory.SP, "Neutral glycosphingolipids [SP05]", "Hex2Cer", "LacCer")
    HEX3CER = (LipidCategory.SP, "Neutral glycosphingolipids [SP05]", "Hex3Cer", "GB3")
    FMC5 = (LipidCategory.SP, "Neutral glycosphingolipids [SP05]", "FMC-5")
    FMC6 = (LipidCategory.SP, "Neutral glycosphingolipids [SP05]", "FMC-6")
    SHEXCER = (LipidCategory.SP, "Acidic glycosphingolipids [SP06]", "SHexCer", "(3'-sulfo)Galbeta-Cer")
    LCB = (LipidCategory.SP, "Sphingoid bases [SP01]", "LCB", "Sphingosine", "So", "Sphingosine-1-phosphate", "SPH")
    LCBP = (LipidCategory.SP, "Sphingoid bases [SP01]", "LCBP", "Sphinganine", "Sa", "Sphingosine-1-phosphate", "S1P", "SPH-P")
    LHexCer = (LipidCategory.SP, "Hexosylsphingosine", "LHexCer", "HexSph")
    EPC = (LipidCategory.SP, "Phosphosphingolipids [SP03]", "EPC", "PE-Cer")
    GB4 = (LipidCategory.SP, "Neutral glycosphingolipids [SP05]", "GB4")
    GD3 = (LipidCategory.SP, "Acidic glycosphingolipids [SP06]", "GD3")
    GM3 = (LipidCategory.SP, "Acidic glycosphingolipids [SP06]", "GM3")
    GM4 = (LipidCategory.SP, "Acidic glycosphingolipids [SP06]", "GM4")
    IPC = (LipidCategory.SP, "Phosphosphingolipids [SP03]", "IPC", "PI-Cer")
    LSM = (LipidCategory.SP, "Ceramides [SP02]", "LSM", "SPC")
    MIP2C = (LipidCategory.SP, "Phosphosphingolipids [SP03]", "M(IP)2C")
    MIPC = (LipidCategory.SP, "Phosphosphingolipids [SP03]", "MIPC")

    
    
    ## Sterol lipids
    ST = (LipidCategory.ST, "Sterols [ST01]", "ST")
    SE = (LipidCategory.ST, "Steryl esters [ST0102]", "SE")
    CH = (LipidCategory.ST, "Cholesterol [LMST01010001]", "CH", "FC", "Cholesterol")
    CHE = (LipidCategory.ST, "Cholesteryl esters [ST0102]", "ChE", "CE", "Cholesteryl ester", "Cholesterol ester")
    
    
    
    ## 	Saccharolipids
    DAT = (LipidCategory.SL, "Acyltrehaloses [SL03]", "DAT")
    AC2SGL = (LipidCategory.SL, "Acyltrehaloses [SL03]", "AC2SGL")
    PAT16 = (LipidCategory.SL, "Acyltrehaloses [SL03]", "PAT16")
    PAT18 = (LipidCategory.SL, "Acyltrehaloses [SL03]", "PAT18")
    


    def get_category(name):
        class_to_category = {key: lipid_class.value[0] for lipid_class in LipidClass for key in lipid_class.value[2:]}
        return class_to_category[name] if name in class_to_category else LipidCategory.UNDEFINED
    
    def get_class(name):
        class_to_class = {key: lipid_class for lipid_class in LipidClass for key in lipid_class.value[2:]}
        return class_to_class[name] if name in class_to_class else LipidClass.UNDEFINED
        