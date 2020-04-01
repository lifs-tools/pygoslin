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

class FattyAcid:

    def __init__(self, name, num_carbon, num_double_bonds, num_hydroxyl, lipid_FA_bond_type, lcb, position):
        self.name = name
        self.position = position
        self.num_carbon = num_carbon
        self.num_double_bonds = num_double_bonds
        self.num_hydroxyl = num_hydroxyl
        self.lipid_FA_bond_type = lipid_FA_bond_type
        self.lcb = lcb
        
        if num_carbon < 2:
            raise ConstraintViolationException("FattyAcid must have at least 2 carbons!")
        
        if position < -1:
            raise ConstraintViolationException("FattyAcid position must be greater or equal to -1 (undefined) or greater or equal to 0 (0 = first position)!")
        
        if num_hydroxyl < 0:
            raise ConstraintViolationException("FattyAcid must have at least 0 hydroxy groups!")
        
    def clone(self, fa):
        self.name = fa.name
        self.position = fa.position
        self.num_carbon = fa.num_carbon
        self.num_double_bonds = fa.num_double_bonds
        self.num_hydroxyl = fa.num_hydroxyl
        self.lipid_FA_bond_type = fa.lipid_FA_bond_type
        self.lcb = fa.lcb
        
    def to_string(self, special_case = False):
        suffix = self.lipid_FA_bond_type.suffix()
        return "%s%i:%i%s%s" % ("O-" if special_case and len(suffix) > 0 else "", self.num_carbon, self.num_double_bonds, ";" + str(self.num_hydroxyl) if self.num_hydroxyl > 0 else "", suffix)