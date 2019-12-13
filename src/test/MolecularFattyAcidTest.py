from ..domain import *

def main():
    
    instanceZero = MolecularFattyAcid("FA1", 0, 2, 0, 0, LipidFaBondType.UNDEFINED, False)
    assert 0 == instanceZero.num_double_bonds
    
    instanceOne = MolecularFattyAcid("FA1", 0, 2, 0, 1, LipidFaBondType.UNDEFINED, False)
    assert 1 == instanceOne.num_double_bonds
        

    try:
        instanceZero = MolecularFattyAcid("FA1", 0, 2, 0, -1, LipidFaBondType.UNDEFINED, False)
    except ConstraintViolationException:
        pass
        

    instance = MolecularFattyAcid("FAX", 0, 2, 0, 0, LipidFaBondType.UNDEFINED, False)
    assert "FAX" == instance.name

    
    instance = MolecularFattyAcid("FAX", 1, 2, 0, 0, LipidFaBondType.UNDEFINED, False);
    assert 1 == instance.position
    

    try:
        instanceZero = MolecularFattyAcid("FA1", -2, 2, 0, 0, LipidFaBondType.UNDEFINED, False)
    except ConstraintViolationException:
        pass


    instance = MolecularFattyAcid("FAX", 1, 2, 0, 0, LipidFaBondType.UNDEFINED, False)
    assert 2 == instance.num_carbon


    try:
        instance = MolecularFattyAcid("FAX", 1, 1, 0, 0, LipidFaBondType.UNDEFINED, False)
    except ConstraintViolationException:
        pass
    

    instance = MolecularFattyAcid("FAX", 1, 2, 1, 0, LipidFaBondType.UNDEFINED, False)
    assert 1 == instance.num_hydroxy


    try:
        instance = MolecularFattyAcid("FAX", 1, 2, -1, 0, LipidFaBondType.UNDEFINED, False)
    except ConstraintViolationException:
        pass
    
    print("test passed without problems")
    
    
    
if __name__ == "__main__":
    main()
