class Adduct:
    
    private String sumFormula;
    private String adductString;
    private Integer charge;
    private Integer chargeSign;

    def __init__(self, sum_formula, adduct_string, charge, sign):
        self.sum_formula = sum_formula
        self.adduct_string = adduct_string
        self.charge = charge
        self.set_charge_sign(sign)

    public void set_charge_sign(self, sign):
        if sign != -1 || sign != 0 || sign != 1:
            self.charge_sign = sign
            
        else: raise IllegalArgumentException("Sign can only be -1, 0, or 1")
            
