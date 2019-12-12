
from enum import Enum

class LipidCategory(Enum):

    UNDEFINED = auto()

    GL = auto() # SLM:000117142 Glycerolipids
    GP = auto() # SLM:000001193 Glycerophospholipids
    SP = auto() # SLM:000000525 Sphingolipids
    ST = auto() # SLM:000500463 Steroids and derivatives
    FA = auto() # SLM:000390054 Fatty acyls and derivatives
