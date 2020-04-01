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


from enum import Enum

class LipidLevel(Enum):
    # Undefined / non-inferable lipid level
    UNDEFINED = 0
    
    # Mediators, Glycerolipids, Glycerophospholipids, Sphingolipids, Steroids, Prenols
    CATEGORY = 1
    
    # Glyerophospholipids -> Glycerophosphoinositols (PI)
    CLASS = 2
    
    # Phosphatidylinositol (16:0) or PI(16:0)
    SPECIES = 3
    
    # Phosphatidylinositol (8:0-8:0) or PI(8:0-8:0)
    MOLECULAR_SUBSPECIES = 4

    # Phosphatidylinositol (8:0/8:0) or PI(8:0/8:0)
    STRUCTURAL_SUBSPECIES = 5
    
    """
    1,2-dioctanoyl-sn-glycero-3-phospho-1D-myo-inositol
    PE(P-18:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z))
    Phosphatidylethanolamine (P-18:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z))
    """
    ISOMERIC_SUBSPECIES = 6
