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


from pygoslin.domain.LipidStructuralSubspecies import LipidStructuralSubspecies
from pygoslin.domain.LipidMolecularSubspecies import LipidMolecularSubspecies
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidSpeciesInfo import LipidSpeciesInfo
from pygoslin.domain.LipidExceptions import *
from pygoslin.domain.LipidLevel import LipidLevel

class LipidIsomericSubspecies(LipidStructuralSubspecies):


    def __init__(self, head_group, fa = []):
        super().__init__(head_group)
        num_carbon = 0
        num_hydroxyl = 0
        num_double_bonds = 0
        lipid_FA_bond_type = LipidFaBondType.UNDEFINED
        if len(fa) > 0:
            lipid_FA_bond_type = LipidFaBondType.ESTER
        
        for fas in fa:
            if fas.name in self.fa:
                raise ConstraintViolationException("FA names must be unique! FA with name %s was already added!" % fas.name)
            
            else:
                self.fa[fas.name] = fas
                self.fa_list.append(fas)
                num_carbon += fas.num_carbon
                num_hydroxyl += fas.num_hydroxyl
                num_double_bonds += fas.num_double_bonds
                
                if lipid_FA_bond_type == LipidFaBondType.ESTER and fas.lipid_FA_bond_type in (LipidFaBondType.ETHER_PLASMANYL, LipidFaBondType.ETHER_PLASMENYL):
                    lipid_FA_bond_type = fas.lipid_FA_bond_type
                    
                elif lipid_FA_bond_type != LipidFaBondType.ESTER and fas.lipid_FA_bond_type in (LipidFaBondType.ETHER_PLASMANYL, LipidFaBondType.ETHER_PLASMENYL):
                    raise ConstraintViolationException("Only one FA can define an ether bond to the head group! Tried to add %s over existing %s" % (fas.lipid_FA_bond_type, lipid_FA_bond_type))
                
        self.info = LipidSpeciesInfo()
        self.info.level = LipidLevel.ISOMERIC_SUBSPECIES
        self.info.num_carbon = num_carbon
        self.info.num_hydroxyl = num_hydroxyl
        self.info.num_double_bonds = num_double_bonds
        self.info.lipid_FA_bond_type = lipid_FA_bond_type
    
    

    def get_lipid_string(self, level = None):
        
        if level == None or level == LipidLevel.ISOMERIC_SUBSPECIES:
            return self.build_lipid_subspecies_name("/")
        
        elif level in (LipidLevel.STRUCTURAL_SUBSPECIES, LipidLevel.MOLECULAR_SUBSPECIES, LipidLevel.CATEGORY, LipidLevel.CLASS, LipidLevel.SPECIES):
            return super().get_lipid_string(level)
        
        else:
            raise Exception("LipidIsomericSubspecies does not know how to create a lipid string for level %s" % level)
