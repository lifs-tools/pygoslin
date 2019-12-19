# Python3 implementation for Goslin

This is the Goslin reference implementation  for Python 3.

> **_NOTE:_**  This is an *early* development version, please use at your own risk and report issues to help improve it!

The pygoslin has been developed with regard the following general issues:

1. Ease the handling with lipid names for developers working on mass spectrometry-based lipidomics tools.
2. Offering a tool to unify all existing dialects of lipids.


## Installation

### Prerequisites
The pygoslin package uses pip to create an isolated and defined build environment. You will need Python >=3.5 to build the package.

```
python3-pip
make (optional)
```

To install the package globally in your python distribution, simply type:

```
make install
```
or

```
[sudo] python setup.py install
```

Be sure that you have root permissions.


### Testing pygoslin

To run the tests, please type:

```
make test
```
or

```
python3 -m unittest pygoslin.tests.MolecularFattyAcidTest
python3 -m unittest pygoslin.tests.ParserTest
```



### Using pygoslin

The two major function within pygoslin are the parsing and printing of lipid names. This example will demonstrate both functions. Open a Python shell and type in:

```
from pygoslin.parser.Parser import GoslinParser
goslin_parser = GoslinParser()
goslin_parser_event_handler = goslin_parser.event_handler

lipid_name = "Cer 18:1;2/12:0"
goslin_parser.parse(lipid_name)
if goslin_parser.word_in_grammar:
    lipid = goslin_parser_event_handler.lipid
    print(lipid.get_lipid_string())
```

