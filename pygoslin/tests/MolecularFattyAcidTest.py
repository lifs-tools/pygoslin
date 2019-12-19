import unittest

from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.MolecularFattyAcid import MolecularFattyAcid
from pygoslin.domain.LipidExceptions import *

class MolecularFattyAcidTest(unittest.TestCase):
    
    def test_instanceZero(self):
        instanceZero = MolecularFattyAcid("FA1", 2, 0, 0, LipidFaBondType.UNDEFINED, False, 0)
        assert 0 == instanceZero.num_double_bonds
        
        
    def test_instanceOne(self):
        instanceOne = MolecularFattyAcid("FA1", 2, 1, 0, LipidFaBondType.UNDEFINED, False, 0)
        assert 1 == instanceOne.num_double_bonds
            

    @unittest.expectedFailure
    def test_wrong_bond_type(self):
        instanceZero = MolecularFattyAcid("FA1", 2, -1, 0, LipidFaBondType.UNDEFINED, False, 0)
            
            
    def test_name(self):
        instance = MolecularFattyAcid("FAX", 2, 0, 0, LipidFaBondType.UNDEFINED, False, 0)
        assert "FAX" == instance.name

        
    def test_position(self):
        instance = MolecularFattyAcid("FAX", 2, 0, 0, LipidFaBondType.UNDEFINED, False, 1);
        assert 1 == instance.position
        

    @unittest.expectedFailure
    def test_wrong_carbon(self):
        instanceZero = MolecularFattyAcid("FA1", 2, 0, 0, LipidFaBondType.UNDEFINED, False, -2)


    def test_carbon(self):
        instance = MolecularFattyAcid("FAX", 2, 0, 0, LipidFaBondType.UNDEFINED, False, 1)
        assert 2 == instance.num_carbon


    @unittest.expectedFailure
    def test_wrong_db(self):
        instance = MolecularFattyAcid("FAX", 1, 0, 0, LipidFaBondType.UNDEFINED, False, 1)
        

    def test_hydroxyl(self):
        instance = MolecularFattyAcid("FAX", 2, 0, 1, LipidFaBondType.UNDEFINED, False, 1)
        assert 1 == instance.num_hydroxyl


    @unittest.expectedFailure
    def test_wrong_hydroxyl(self):
        instance = MolecularFattyAcid("FAX", 2, 0, -1, LipidFaBondType.UNDEFINED, False, 1)
    
if __name__ == '__main__':
    unittest.main()