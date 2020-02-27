from setuptools import setup, find_packages
import setuptools

pyx_support = True
try:
    from Cython.Build import cythonize
except:
    pyx_support = False
    print("Warning: cython module is not installed, parsing performance will be lower since pure python code will be applied.")


import os

os.environ["CFLAGS"] = "-O3 -Wall -std=c++11" 

setup(
    name = 'pygoslin',
    version = '1.0',
    url = 'https://gitlib.isas.de/kopczynski/goslin',
    license = 'MIT',
    author = 'Dominik Kopczynski',
    author_email = 'dominik.kopczynski@isas.de',
    description = 'Python implementation for Goslin',
    long_description = open('README.md', encoding='utf-8').read(),
    packages = setuptools.find_packages(),
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    ext_modules= cythonize("pygoslin/parser/ParserCore.pyx", language="c++") if pyx_support else None,
    setup_requires = ["pytest-runner"],
    tests_require = ["pytest"],
    python_requires='>=3.5',
    include_package_data=True,
    package_data={
        '': ['data/goslin/*.g4', 'data/goslin/*.G4'], # If any package contains *.G4 files, include them
    }
)
