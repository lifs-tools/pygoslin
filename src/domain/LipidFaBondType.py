from enum import Enum

class LipidFaBondType(Enum):
    UNDEFINED = auto()
    ESTER = auto()
    ETHER_PLASMANYL = auto()
    ETHER_PLASMENYL = auto()

    def suffix(self):
        if self == ETHER_PLASMANYL: return "a"
        elif self == ETHER_PLASMENYL: return "p"
        else: return ""
    

    def double_bond_correction(self):
        return 1 if self == ETHER_PLASMENYL else 0