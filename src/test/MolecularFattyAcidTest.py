import unittest

from domain.LipidFaBondType import LipidFaBondType
from domain.MolecularFattyAcid import MolecularFattyAcid
from domain.LipidExceptions import *

class MolecularFattyAcidTest(unittest.TestCase):
    
    def test_instanceZero(self):
        instanceZero = MolecularFattyAcid("FA1", 0, 2, 0, 0, LipidFaBondType.UNDEFINED, False)
        assert 0 == instanceZero.num_double_bonds
        
        
    def test_instanceOne(self):
        instanceOne = MolecularFattyAcid("FA1", 0, 2, 0, 1, LipidFaBondType.UNDEFINED, False)
        assert 1 == instanceOne.num_double_bonds
            

    @unittest.expectedFailure
    def test_wrong_bond_type(self):
        instanceZero = MolecularFattyAcid("FA1", 0, 2, 0, -1, LipidFaBondType.UNDEFINED, False)
            
            
    def test_name(self):
        instance = MolecularFattyAcid("FAX", 0, 2, 0, 0, LipidFaBondType.UNDEFINED, False)
        assert "FAX" == instance.name

        
    def test_position(self):
        instance = MolecularFattyAcid("FAX", 1, 2, 0, 0, LipidFaBondType.UNDEFINED, False);
        assert 1 == instance.position
        

    @unittest.expectedFailure
    def test_wrong_carbon(self):
        instanceZero = MolecularFattyAcid("FA1", -2, 2, 0, 0, LipidFaBondType.UNDEFINED, False)


    def test_carbon(self):
        instance = MolecularFattyAcid("FAX", 1, 2, 0, 0, LipidFaBondType.UNDEFINED, False)
        assert 2 == instance.num_carbon


    @unittest.expectedFailure
    def test_wrong_db(self):
        instance = MolecularFattyAcid("FAX", 1, 1, 0, 0, LipidFaBondType.UNDEFINED, False)
        

    def test_hydroxyl(self):
        instance = MolecularFattyAcid("FAX", 1, 2, 1, 0, LipidFaBondType.UNDEFINED, False)
        assert 1 == instance.num_hydroxy


    @unittest.expectedFailure
    def test_wrong_hydroxyl(self):
        instance = MolecularFattyAcid("FAX", 1, 2, -1, 0, LipidFaBondType.UNDEFINED, False)
    
if __name__ == '__main__':
    unittest.main()