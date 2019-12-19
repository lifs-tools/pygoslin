class Adduct:
    
    def __init__(self, sum_formula, adduct_string, charge, sign):
        self.sum_formula = sum_formula
        self.adduct_string = adduct_string
        self.charge = charge
        self.set_charge_sign(sign)

    def set_charge_sign(self, sign):
        if sign != -1 or sign != 0 or sign != 1:
            self.charge_sign = sign
            
        else: raise IllegalArgumentException("Sign can only be -1, 0, or 1")
            
    def get_lipid_string(self):
        if self.charge == 0: return "[M]"
        
        return "[M%s%s]%i%s" % (self.sum_formula, self.adduct_string, self.charge, "+" if self.charge_sign > 0 else "-")