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

class LipidFaBondType(Enum):
    UNDEFINED = 0
    ESTER = 1
    ETHER_PLASMANYL = 2
    ETHER_PLASMENYL = 3
    ETHER = 4
    LCB_REGULAR = 5
    LCB_EXCEPTION = 6
    AMIDE = 7
    ETHER_UNSPECIFIED = 8

    def prefix(self):
        if self == self.ETHER_PLASMANYL: return "O-"
        elif self == self.ETHER_PLASMENYL: return "P-"
        else: return ""
    

    def double_bond_correction(self):
        return 1 if self == self.ETHER_PLASMENYL else 0
    
    
