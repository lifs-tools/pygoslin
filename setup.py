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


from setuptools import setup, find_packages
import setuptools

pyx_support = True
try:
    from Cython.Build import cythonize
except:
    pyx_support = False
    print("Warning: cython module is not installed, parsing performance will be lower since pure python code will be applied.")


import os

#os.environ["CFLAGS"] = "-O3 -Wall" 
#os.environ["CPPFLAGS"] = "-O3 -Wall -std=c++11" 

setup(
    name = 'pygoslin',
    version = '1.1.2',
    url = 'https://gitlib.isas.de/kopczynski/goslin',
    license = 'MIT',
    author = 'Dominik Kopczynski',
    author_email = 'dominik.kopczynski@isas.de',
    description = 'Python implementation for Goslin',
    long_description = open('README.md', encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    packages = setuptools.find_packages(),
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    #ext_modules= cythonize("pygoslin/parser/ParserCore.pyx", language_level = 3) if pyx_support else None,
    setup_requires = ["pytest-runner"],
    tests_require = ["pytest"],
    python_requires = '>=3.5',
    include_package_data = True,
    package_data = {
        '': ['data/goslin/*.g4', 'data/goslin/*.G4', 'data/goslin/lipid-list.csv'],
    }
)
