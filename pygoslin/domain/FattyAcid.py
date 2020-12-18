"""
MIT License

Copyright (c) 2020 Dominik Kopczynski   -   dominik.kopczynski {at} isas.de
                   Nils Hoffmann  -  nils.hoffmann {at} isas.de

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.Element import Element
from pygoslin.domain.LipidFaBondType import LipidFaBondType

class FattyAcid:

    def __init__(self, name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position, double_bond_positions = None):
        self.name = name
        self.position = position
        self.num_carbon = num_carbon
        self.num_double_bonds = num_double_bonds
        self.num_hydroxyl = num_hydroxyl
        self.lipid_FA_bond_type = lipid_FA_bond_type
        self.lcb = lcb
        self.double_bond_positions = {key: double_bond_positions[key] for key in double_bond_positions} if double_bond_positions != None else {}
        
        if num_carbon < 2:
            raise ConstraintViolationException("FattyAcid must have at least 2 carbons!")
        
        if num_double_bonds < 0:
            raise ConstraintViolationException("FattyAcid must have at least 0 double bonds!")
        
        if num_hydroxyl < 0:
            raise ConstraintViolationException("FattyAcid must have at least 0 hydroxy groups!")
        
        if position < 0:
            raise ConstraintViolationException("FattyAcid must be at least 0 at position 0!")
        
        
    def clone(self, fa):
        self.name = fa.name
        self.position = fa.position
        self.num_carbon = fa.num_carbon
        self.num_double_bonds = fa.num_double_bonds
        self.num_hydroxyl = fa.num_hydroxyl
        self.lipid_FA_bond_type = fa.lipid_FA_bond_type
        self.lcb = fa.lcb
        self.double_bond_positions = {key: value for key, value in fa.double_bond_positions.items()}
        
        
    def to_string(self, special_case = False):
        if self.num_carbon == 0: return None
        
        suffix = self.lipid_FA_bond_type.suffix()
        dbp = self.double_bond_positions
        db_positions = ["%i%s" % (k, dbp[k]) for k in sorted(dbp.keys())]
        db_pos = "(%s)" % ",".join(db_positions) if len (dbp) > 0 else ""
        
        return "%s%i:%i%s%s%s" % ("O-" if special_case and len(suffix) > 0 else "", self.num_carbon, self.num_double_bonds, db_pos, ";" + str(self.num_hydroxyl) if self.num_hydroxyl > 0 else "", suffix)


    def assemble_fa_smiles(self, fa_structure):
        
        components = []
        last_db = -1
        for i, token in enumerate(fa_structure):
            if token == "C":
                components.append("C")
            elif token == "C=":
                if i - 1 == last_db:
                    components.append("C=")
                else:
                    components.append("\\C=")
            elif token == "=C":
                components.append("C\\")
                last_db = i
            elif token == "/C=":
                if i - 1 == last_db:
                    components.append("C=")
                else:
                    components.append("/C=")
            elif token == "=C/":
                components.append("C/")
                last_db = i
            elif token == "\\C=":
                if i - 1 == last_db:
                    components.append("C=")
                else:
                    components.append("\\C=")
            elif token == "=C\\":
                components.append("C\\")
                last_db = i
            elif token == "C=O":
                components.append("C(=O)")
            elif token == "C-OH":
                components.append("C(O)")
                
        return "".join(components)





    def fa_smiles(self):
        
        if self.num_carbon == 0: return "H"
        
        # check if defined double bond positions are possible
        double_bond_positions = {}
        not_isomeric = True
        if self.num_double_bonds > 0 and len(self.double_bond_positions) > 0:
            n_db = len(self.double_bond_positions) + (1 if self.lipid_FA_bond_type == LipidFaBondType.ETHER_PLASMENYL else 0)
            if self.num_double_bonds != n_db:
                raise ConstraintViolationException("Double bond number (%i) is unequal to number of provided double bond positions (%i)" % (self.num_double_bonds, n_db))
            
            not_isomeric = False
            dbpk = sorted(list(k for k in self.double_bond_positions))
            for i in range(1, len(dbpk)):
                if dbpk[i] - dbpk[i - 1] <= 1:
                    raise ConstraintViolationException("Impossible double bond positions %s" % str(dbpk))
                elif dbpk[i] - dbpk[i - 1] == 2 and self.double_bond_positions[dbpk[i - 1]] == "E" and self.double_bond_positions[dbpk[i]] == "Z":
                    raise ConstraintViolationException("Impossible double bond positions %i%s and %i%s" % (dbpk[i - 1], self.double_bond_positions[dbpk[i - 1]], dbpk[i], self.double_bond_positions[dbpk[i]]))
        
            if dbpk[-1] == self.num_carbon:
                raise ConstraintViolationException("Impossible double bond position %i with fatty acyl length %i" % (dbpk[-1], self.num_carbon))
        
            double_bond_positions = self.double_bond_positions
        
        
        fa_structure = ["C"] * self.num_carbon
        remaining_C = self.num_carbon
        remaining_OH = self.num_hydroxyl
        
        if self.lcb:
            remaining_C -= 3
            fa_structure = fa_structure[:-3]
            remaining_OH -= 2
            if not_isomeric: double_bond_positions = {i * 2 + 1: "E" for i in range(self.num_double_bonds)} 
        
        else:
            if self.lipid_FA_bond_type == LipidFaBondType.ESTER:
                fa_structure[0] = "C=O"
                remaining_C -= 1
                if not_isomeric: double_bond_positions = {i * 2 + 2: "E" for i in range(self.num_double_bonds)} 
            elif self.lipid_FA_bond_type == LipidFaBondType.ETHER_PLASMANYL:
                remaining_C -= 1
                if not_isomeric: double_bond_positions = {i * 2 + 1: "E" for i in range(self.num_double_bonds)}
            elif self.lipid_FA_bond_type == LipidFaBondType.ETHER_PLASMENYL:
                if not_isomeric: double_bond_positions = {i * 2 + 1: "E" for i in range(self.num_double_bonds)}
                else: double_bond_positions[1] = "E"
        
        
        
        if remaining_C < self.num_double_bonds * 2:
            raise ConstraintViolationException("Impossible composition with double bond (%i) hydroxylation (%i) numbers on fatty acyl length %i" % (self.num_double_bonds, self.num_hydroxyl, self.num_carbon))
        
        
        for dbp, ez in double_bond_positions.items():
            if ez == "E":
                fa_structure[dbp - 1] = "\\C="
                fa_structure[dbp] = "=C\\"
            elif ez == "Z":
                fa_structure[dbp - 1] = "\\C="
                fa_structure[dbp] = "=C/"
            else:
                raise ConstraintViolationException("Error unknown double bond E/Z sign '%s'" % ez)       
    
            
                
        if remaining_OH > 0:
            nh = remaining_OH
            for i, c in enumerate(fa_structure):
                if c == "C": 
                    fa_structure[i] = "C-OH"
                    nh -= 1
                if nh == 0: break
            
            if nh > 0:
                raise ConstraintViolationException("Impossible composition with double bond (%i) hydroxylation (%i) numbers on fatty acyl length %i" % (self.num_double_bonds, self.num_hydroxyl, self.num_carbon))      
                
                
        return self.assemble_fa_smiles(fa_structure)
                
                
                

    def get_elements(self):
        elements = {e: 0 for e in Element}
        
        if not self.lcb:
        
            if self.num_carbon > 0 or self.num_double_bonds > 0:
                
                elements[Element.C] = self.num_carbon # carbon
                if self.lipid_FA_bond_type == LipidFaBondType.ESTER:
                    elements[Element.H] = (2 * self.num_carbon - 1 - 2 * self.num_double_bonds) # hydrogen
                    elements[Element.O] = (1 + self.num_hydroxyl) # oxygen
                
                elif self.lipid_FA_bond_type == LipidFaBondType.ETHER_PLASMENYL:
                    elements[Element.H] = (2 * self.num_carbon - 1 - 2 * self.num_double_bonds + 2) # hydrogen
                    elements[Element.O] = self.num_hydroxyl # oxygen
                
                elif self.lipid_FA_bond_type == LipidFaBondType.ETHER_PLASMANYL:
                    elements[Element.H] = ((self.num_carbon + 1) * 2 - 1 - 2 * self.num_double_bonds) # hydrogen
                    elements[Element.O] = self.num_hydroxyl # oxygen
                    
                else:
                    raise LipidException("Mass cannot be computed for fatty acyl chain with bond type: %s" % self.lipid_FA_bond_type)
                
                
                
        else:
            # long chain base
            elements[Element.C] = self.num_carbon # carbon
            elements[Element.H] = (2 * (self.num_carbon - self.num_double_bonds) + 1) # hydrogen
            elements[Element.O] = self.num_hydroxyl # oxygen
            elements[Element.N] = 1 # nitrogen
        return elements
