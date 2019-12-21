from setuptools import setup, find_packages
import setuptools
 
setup(
    name = 'pygoslin',
    version = '1.0',
    url = 'https://gitlib.isas.de/kopczynski/goslin',
    license = 'MIT',
    author = 'Dominik Kopczynski',
    author_email = 'dominik.kopczynski@isas.de',
    description = 'Python implementation for Goslin',
    long_description = open('README.md').read(),
    packages = setuptools.find_packages(),
#      packages = find_packages(exclude = ['tests']),
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    setup_requires = ["pytest-runner"],
    tests_require = ["pytest"],
    python_requires='>=3.5',
    include_package_data=True,
    package_data={
        '': ['data/goslin/*.g4', 'data/goslin/*.G4'], # If any package contains *.G4 files, include them
    }
)