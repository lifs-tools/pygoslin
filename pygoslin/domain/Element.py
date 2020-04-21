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

class Element(Enum):
    C = 0
    C13 = 1
    H = 2
    H2 = 3
    N = 4
    N15 = 5
    O = 6
    O17 = 7
    O18 = 8
    P = 9
    P32 = 10
    S = 11
    S34 = 12
    S33 = 13
    
element_positions = {"C": Element.C,
                    "H": Element.H,
                    "N": Element.N,
                    "O": Element.O,
                    "P": Element.P,
                    "P'": Element.P32,
                    "S": Element.S,
                    "S'": Element.S34,
                    "S''": Element.S33,
                    "H'": Element.H2,
                    "C'": Element.C13,
                    "N'": Element.N15,
                    "O'": Element.O17,
                    "O''": Element.O18,
                    "2H": Element.H2,
                    "13C": Element.C13,
                    "15N": Element.N15,
                    "17O": Element.O17,
                    "18O": Element.O18,
                    "32P": Element.P32,
                    "34S": Element.S34,
                    "33S": Element.S33,
                    "H2": Element.H2,
                    "C13": Element.C13,
                    "N15": Element.N15,
                    "O17": Element.O17,
                    "O18": Element.O18,
                    "P32": Element.P32,
                    "S34": Element.S34,
                    "S33": Element.S33
                     }
