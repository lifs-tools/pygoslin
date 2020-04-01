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

class LipidCategory(Enum):

    UNDEFINED = 0
    GL = 1 # SLM:000117142 Glycerolipids
    GP = 2 # SLM:000001193 Glycerophospholipids
    SP = 3 # SLM:000000525 Sphingolipids
    ST = 4 # SLM:000500463 Steroids and derivatives
    FA = 5 # SLM:000390054 Fatty acyls and derivatives
    SL = 6 # Saccharolipids
    PK = 7 # Polyketides
    
category_string_to_category = {"GL": LipidCategory.GL,
                      "GP": LipidCategory.GP,
                      "SP": LipidCategory.SP,
                      "ST": LipidCategory.ST,
                      "FA": LipidCategory.FA,
                      "SL": LipidCategory.SL,
                      "PK": LipidCategory.PK,
                      "UNDEFINED": LipidCategory.UNDEFINED
                      }
    