# Python3 implementation for Goslin

This is the Goslin reference implementation for Python 3.

> **_NOTE:_**  This is an *early* development version, please use at your own risk and report issues to help improve it!

The pygoslin has been developed with regard the following general issues:

1. Ease the handling with lipid names for developers working on mass spectrometry-based lipidomics tools.
2. Offering a tool to unify all existing dialects of lipid names.


## Installation

### Prerequisites
The pygoslin package uses pip to create an isolated and defined build environment. You will need Python >=3.5 to build the package.

```
python3-pip
make (optional)
```

To install the package globally in your Python distribution, simply type:

```
[sudo] make install
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

The two major functions within pygoslin are the parsing and printing of lipid names. You have several options, to access these functions. This example will demonstrate both functions the easiest way. Open a Python shell and type in:


```
from pygoslin.parser.Parser import LipidParser
lipid_parser = LipidParser()

lipid_name = "PE 16:1-12:0"
lipid = lipid_parser.parse(lipid_name)

if lipid != None: print(lipid.get_lipid_string())

```

Be aware, that this method is subsequencially using all available grammars until a lipid name can be described by the first grammar. Currently, three grammars are available, namely:
```
Goslin
GoslinFragment
LipidMaps
```

To use a specific grammar / parser, you can use the following code:


```
# using solely the Goslin parser
from pygoslin.parser.Parser import GoslinParser
goslin_parser = GoslinParser()

lipid_name = "Cer 18:1;2/12:0"
lipid = goslin_parser.parse(lipid_name)
if lipid != None:
    print(lipid.get_lipid_string())
    
    
    
# using solely the Goslin Fragment parser
from pygoslin.parser.Parser import GoslinFragmentParser
goslin_fragment_parser = GoslinFragmentParser()

lipid_name = "DAG 18:1-12:0 - -(H2O)"
lipid = goslin_fragment_parser.parse(lipid_name)
if lipid != None:
    print(lipid.get_lipid_string())
    print(lipid.get_lipid_fragment_string())
    
    
    
# using solely the LipidMaps parser
from pygoslin.parser.Parser import LipidMapsParser
lipid_maps_parser = LipidMapsParser()

lipid_name = "Cer(d18:1/12:0)"
lipid = lipid_maps_parser.parse(lipid_name)
if lipid != None:
    print(lipid.get_lipid_string())
```
To be as generic as possible, no treatment of validation of the fragment is conducted within the GoslinFragmentParser.



### Supported lipids
<table>
    <tr><th>Lipid category</th><th>Lipid class</th><th>Abbreviation</th></tr>

    <tr><td rowspan="30">Glycerophospholipids (GP)</td><td>Bismonoacylglycerophosphate</td><td>BMP</td></tr>
    <tr><td>CDP-diacylglycerol</td><td>CDP-DAG</td></tr>
    <tr><td>Cardiolipin</td><td>CL</td></tr>
    <tr><td>Dimethylphosphatidylethanolamine</td><td>DMPE</td></tr>
    <tr><td>Lysophosphatidic acid</td><td>LPA</td></tr>
    <tr><td>Lysophophatidylcholine</td><td>LPC</td></tr>
    <tr><td rowspan="2">Ether lysophosphatidic acid</td><td>LPC O-a</td></tr>
    <tr><td>LPC O-p</td></tr>
    <tr><td>Lysophosphatidylethanolamine</td><td>LPE</td></tr>
    <tr><td rowspan="2">Ether lysophosphatidylethanolamine</td><td>LPE O-a</td></tr>
    <tr><td>LPE O-p</td></tr>
    <tr><td>Lysophosphatidylglycerol</td><td>LPG</td></tr>
    <tr><td>Lysophosphatidylinositol</td><td>LPI</td></tr>
    <tr><td>Lysophosphatidylserine</td><td>LPS</td></tr>
    <tr><td>Monolysocardiolipin</td><td>MLCL</td></tr>
    <tr><td>Monomethylphosphatidylethanolamine</td><td>MMPE</td></tr>
    <tr><td>Phosphatidic acid</td><td>PA</td></tr>
    <tr><td>Phosphatidylcholine</td><td>PC</td></tr>
    <tr><td rowspan="2">Ether phosphatidylcholine</td><td>PC O-a</td></tr>
    <tr><td>PC O-p</td></tr>
    <tr><td>Phosphatidylethanolamine</td><td>PE</td></tr>
    <tr><td rowspan="2">Ether phosphatidylethanolamine</td><td>PE O-a</td></tr>
    <tr><td>PE O-p</td></tr>
    <tr><td>Phosphatidylethanol</td><td>PEt</td></tr>
    <tr><td>Phosphatidylglycerol</td><td>PG</td></tr>
    <tr><td>Phosphatidylinositol</td><td>PI</td></tr>
    <tr><td>Phosphatidylinositolphosphate</td><td>PIP</td></tr>
    <tr><td>Phosphatidylinositolbisphosphate</td><td>PIP2</td></tr>
    <tr><td>Phosphatidylinositoltrisphosphate</td><td>PIP3</td></tr>
    <tr><td>Phosphatidylserine</td><td>PS</td></tr>

</table>
