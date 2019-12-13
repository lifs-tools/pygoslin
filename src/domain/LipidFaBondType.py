from enum import Enum

class LipidFaBondType(Enum):
    UNDEFINED = 0
    ESTER = 1
    ETHER_PLASMANYL = 2
    ETHER_PLASMENYL = 3

    def suffix(self):
        if self == ETHER_PLASMANYL: return "a"
        elif self == ETHER_PLASMENYL: return "p"
        else: return ""
    

    def double_bond_correction(self):
        return 1 if self == ETHER_PLASMENYL else 0
    