
from enum import Enum

class LipidCategory(Enum):

    UNDEFINED = 0

    GL = 1 # SLM:000117142 Glycerolipids
    GP = 2 # SLM:000001193 Glycerophospholipids
    SP = 3 # SLM:000000525 Sphingolipids
    ST = 4 # SLM:000500463 Steroids and derivatives
    FA = 5 # SLM:000390054 Fatty acyls and derivatives
