
from enum import Enum

class LipidCategory(Enum):

    UNDEFINED = 0
    GL = 1 # SLM:000117142 Glycerolipids
    GP = 2 # SLM:000001193 Glycerophospholipids
    SP = 3 # SLM:000000525 Sphingolipids
    ST = 4 # SLM:000500463 Steroids and derivatives
    FA = 5 # SLM:000390054 Fatty acyls and derivatives
    SL = 6 # Saccharolipids
    
category_string_to_category = {"GL": LipidCategory.GL,
                      "GP": LipidCategory.GP,
                      "SP": LipidCategory.SP,
                      "ST": LipidCategory.ST,
                      "FA": LipidCategory.FA,
                      "SL": LipidCategory.SL,
                      "UNDEFINED": LipidCategory.UNDEFINED
                      }
    