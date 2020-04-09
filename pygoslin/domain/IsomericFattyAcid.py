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


from pygoslin.domain.StructuralFattyAcid import StructuralFattyAcid

class IsomericFattyAcid(StructuralFattyAcid):

    def __init__(self, name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position, double_bond_positions):
        super().__init__(name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position)
        self.double_bond_positions = {key: double_bond_positions[key] for key in double_bond_positions}

    def clone(self, fa):
        self.double_bond_positions = {key: double_bond_positions[key] for key in fa.double_bond_positions}
        super().clone(fa)
        
    
    def to_string(self, special_case = False):
        
        suffix = self.lipid_FA_bond_type.suffix()
        dbp = self.double_bond_positions
        db_positions = ["%i%s" % (k, dbp[k]) for k in sorted(dbp.keys())]
        db_pos = "(%s)" % ",".join(db_positions) if len (dbp) > 0 else ""
        
        return "%s%i:%i%s%s%s" % ("O-" if special_case and len(suffix) > 0 else "", self.num_carbon, self.num_double_bonds, db_pos, ";" + str(self.num_hydroxyl) if self.num_hydroxyl > 0 else "", suffix)
