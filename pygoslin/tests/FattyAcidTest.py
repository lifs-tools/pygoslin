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


import unittest

from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.FattyAcid import FattyAcid
from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.Element import Element

class FattyAcidTest(unittest.TestCase):
    
    def test_instanceZero(self):
        instanceZero = FattyAcid("FA1", 2, 0, 0, LipidFaBondType.UNDEFINED, False, 0)
        assert 0 == instanceZero.num_double_bonds
        
        
    def test_instanceOne(self):
        instanceOne = FattyAcid("FA1", 2, 1, 0, LipidFaBondType.UNDEFINED, False, 0)
        assert 1 == instanceOne.num_double_bonds
            

    @unittest.expectedFailure
    def test_wrong_bond_type(self):
        instanceZero = FattyAcid("FA1", 2, -1, 0, LipidFaBondType.UNDEFINED, False, 0)
            
            
    def test_name(self):
        instance = FattyAcid("FAX", 2, 0, 0, LipidFaBondType.UNDEFINED, False, 0)
        assert "FAX" == instance.name

        
    def test_position(self):
        instance = FattyAcid("FAX", 2, 0, 0, LipidFaBondType.UNDEFINED, False, 1);
        assert 1 == instance.position
        

    @unittest.expectedFailure
    def test_wrong_position(self):
        instanceZero = FattyAcid("FA1", 20, 2, 0, LipidFaBondType.UNDEFINED, False, -2)


    def test_carbon(self):
        instance = FattyAcid("FAX", 20, 2, 0, LipidFaBondType.ESTER, False, 1)
        assert 20 == instance.num_carbon
        elements = instance.get_elements()
        assert elements[Element.C] == 20
        assert elements[Element.H] == 35
        assert elements[Element.O] == 1
        assert elements[Element.N] == 0


    @unittest.expectedFailure
    def test_wrong_db(self):
        instance = FattyAcid("FAX", 1, 0, 0, LipidFaBondType.UNDEFINED, False, 1)
        

    def test_hydroxyl(self):
        instance = FattyAcid("FAX", 2, 0, 1, LipidFaBondType.UNDEFINED, False, 1)
        assert 1 == instance.num_hydroxyl


    @unittest.expectedFailure
    def test_wrong_hydroxyl(self):
        instance = FattyAcid("FAX", 2, 0, -1, LipidFaBondType.UNDEFINED, False, 1)
    
if __name__ == '__main__':
    unittest.main()
