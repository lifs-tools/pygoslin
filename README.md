# Python3 implementation for Goslin
[![Build Status](https://github.com/lifs-tools/pygoslin/actions/workflows/python-package.yml/badge.svg)](https://github.com/lifs-tools/pygoslin/actions/workflows/python-package.yml)

This is the Goslin reference implementation for Python 3.

The pygoslin has been developed with regard the following general issues:

1. Ease the handling with lipid names for developers working on mass spectrometry-based lipidomics tools.
2. Offering a tool to unify all existing dialects of lipid names.

## Related Projects

- [This project](https://github.com/lifs-tools/pygoslin)
- [C++ implementation](https://github.com/lifs-tools/cppgoslin)
- [Goslin grammars and test files](http://github.com/lifs-tools/goslin)
- [Java implementation](https://github.com/lifs-tools/jgoslin)
- [R implementation](https://github.com/lifs-tools/rgoslin)
- [Webapplication and REST API](https://github.com/lifs-tools/goslin-webapp)

## Installation

### Installation via pip

The easiest way to install pygoslin is simply via pip:
```
pip3 install pygoslin
```

### Regular installation: prerequisites
The pygoslin package uses pip to create an isolated and defined build environment. You will need Python >=3.5 to build the package.

```
python3-pip
cython (module for Python 3)
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
python3 -m unittest pygoslin.tests.FattyAcidTest
python3 -m unittest pygoslin.tests.ParserTest
python3 -m unittest pygoslin.tests.SwissLipidsTest
python3 -m unittest pygoslin.tests.GoslinTest
python3 -m unittest pygoslin.tests.LipidMapsTest
python3 -m unittest pygoslin.tests.HmdbTest

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
SwissLipids
HMDB
```

To use a specific grammar / parser, you can use the following code:


```
# using solely the Goslin parser
from pygoslin.parser.Parser import GoslinParser
goslin_parser = GoslinParser()

lipid_name = "Cer 18:1;2/12:0"
try:
    lipid = goslin_parser.parse(lipid_name)
    print(lipid.get_lipid_string())
except Exception as e:
    print(e)
    



# using solely the Goslin Fragment parser
from pygoslin.parser.Parser import GoslinFragmentParser
goslin_fragment_parser = GoslinFragmentParser()

lipid_name = "Cer 18:1;2/12:0"
try:
    lipid = goslin_fragment_parser.parse(lipid_name)
    print(lipid.get_lipid_string())
except Exception as e:
    print(e)
    
    


# using solely the LipidMaps parser
from pygoslin.parser.Parser import LipidMapsParser
lipid_maps_parser = LipidMapsParser()

lipid_name = "Cer(d18:1/12:0)"
try:
    lipid = lipid_maps_parser.parse(lipid_name)
    print(lipid.get_lipid_string())
except Exception as e:
    print(e)




# using solely the SwissLipids parser
from pygoslin.parser.Parser import SwissLipidsParser
swiss_lipids_parser = SwissLipidsParser()

lipid_name = "Cer(d18:1/12:0)"
try:
    lipid = swiss_lipids_parser.parse(lipid_name)
    print(lipid.get_lipid_string())
except Exception as e:
    print(e)
    


# using solely the HMDB parser
from pygoslin.parser.Parser import HmdbParser
hmdb_parser = HmdbParser()

lipid_name = "Cer(d18:1/12:0)"
try:
    lipid = hmdb_parser.parse(lipid_name)
    print(lipid.get_lipid_string())
except Exception as e:
    print(e)
```
To be as generic as possible, no treatment of validation of the fragment is conducted within the GoslinFragmentParser.



### Supported lipids
<table>
<tr><th>Lipid category</th><th>Lipid class</th><th>Abbreviation</th></tr>

<tr><td rowspan="48">Glycerophospholipids (GP)</td><td>Bismonoacylglycerophosphate</td><td>BMP / LBPA</td></tr>
<tr><td>CDP-diacylglycerol</td><td>CDP-DAG</td></tr>
<tr><td>Dimethylphosphatidylethanolamine</td><td>DMPE</td></tr>
<tr><td>Monomethylphosphatidylethanolamine</td><td>MMPE</td></tr>
<tr><td>Phosphatidylinositol mannoside inositol phosphate</td><td>PIMIP</td></tr>
<tr><td>Lyso-CDP-diacylglycerol</td><td>LCDPDAG</td></tr>
<tr><td>Lysodimethylphosphatidylethanolamine</td><td>LDMPE</td></tr>
<tr><td>Lysomonomethylphosphatidylethanolamine</td><td>LMMPE</td></tr>
<tr><td>Lysophosphatidylinositol- mannosideinositolphosphate</td><td>LPIMIP</td></tr>
<tr><td>Lysophosphatidylinositol-glucosamine</td><td>LPIN</td></tr>
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
<tr><td>Cardiolipin</td><td>CL</td></tr>
<tr><td>Glycerophosphoglycerophosphoglycerols</td><td>DLCL</td></tr>
<tr><td>Monolysocardiolipin</td><td>MLCL</td></tr>
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
<tr><td>Phosphatidylinositolphosphate</td><td>PIP / PIP[3'] / PIP[4'] / PIP[5']</td></tr>
<tr><td>Phosphatidylinositolbisphosphate</td><td>PIP2 / PIP2[3',4'] / PIP2[3',5'] / PIP2[4',5']</td></tr>
<tr><td>Phosphatidylinositoltrisphosphate</td><td>PIP3 / PIP3[3',4',5']</td></tr>
<tr><td>Phosphatidylserine</td><td>PS</td></tr>
<tr><td>Phosphatidylinositol mannoside</td><td>PIM / PIM1 / PIM2 / PIM3<br>PIM4 / PIM5 / PIM6</td></tr>
<tr><td>Lysophosphatidylinositol mannoside</td><td>LPIM / LPIM1 / LPIM2 / LPIM3<br>LPIM4 / LPIM5 / LPIM6</td></tr>
<tr><td>Phosphatidylglycerol phosphate</td><td>PGP</td></tr>
<tr><td>Diacylglyceropyrophosphate</td><td>PPA</td></tr> 
<tr><td>Diacylglycosylglycerophospholipid</td><td>Glc-GP / 6-Ac-Glc-GP</td></tr>
<tr><td>Diacylglycerophosphonocholine</td><td>PnC</td></tr>
<tr><td>Diacylglycerophosphonoethanolamine</td><td>PnE</td></tr>
<tr><td>Diacylglycerophosphoethanolamine</td><td>PE-NMe / PE-NMe2</td></tr>
<tr><td>Diacylglycerophosphomonoradylglycerol</td><td>SLBPA</td></tr>
<tr><td>N-acylphosphatidylethanolamine</td><td>NAPE</td></tr>
<tr><td></td><td>CPA</td></tr>


<tr><td rowspan="19">Sphingolipids (SP)</td><td>Ceramide</td><td>Cer</td></tr>
<tr><td>Ceramide phosphate</td><td>CerP</td></tr>
<tr><td>Ethanolaminephosphoceramide</td><td>EPC</td></tr>
<tr><td>Ganglioside GB3</td><td>GB3</td></tr>
<tr><td>Ganglioside GB4</td><td>GB4</td></tr>
<tr><td>Ganglioside GD3</td><td>GD3</td></tr>
<tr><td>Ganglioside GM3</td><td>GM3</td></tr>
<tr><td>Ganglioside GM4</td><td>GM4</td></tr>
<tr><td>Dihexosylceramide</td><td>Hex2Cer</td></tr>
<tr><td>Hexosylceramide</td><td>HexCer</td></tr>
<tr><td>Inositolphosphoceramide</td><td>IPC</td></tr>
<tr><td>Long-chain base</td><td>LCB</td></tr>
<tr><td>Long-chain base phosphate</td><td>LCBP</td></tr>
<tr><td>Lysomonohexosylceramide</td><td>LHexCer</td></tr>
<tr><td>Lysosphingomyelin</td><td>LSM</td></tr>
<tr><td>Mannosyldiinositolphosphoceramide</td><td>M(IP)2C</td></tr>
<tr><td>Mannosylinositolphosphoceramide</td><td>MIPC</td></tr>
<tr><td>Sulfatide</td><td>SHexCer</td></tr>
<tr><td>Sphingomyelin</td><td>SM</td></tr>

<tr><td rowspan="2">Sterol lipids (ST)</td><td>Cholesterol</td><td>Ch</td></tr>
<tr><td>Cholesteryl ester</td><td>ChE</td></tr>

<tr><td rowspan="7">Glycerolipids (GL)</td><td>Diacylglycerol</td><td>DAG</td></tr>
<tr><td>Digalactosyldiacylglycerol</td><td>DGDG</td></tr>
<tr><td>Monoacylglycerol</td><td>MAG</td></tr>
<tr><td>Monogalactosyldiacylglycerol</td><td>MGDG</td></tr>
<tr><td>Sulfoquinovosyl monoacylglycerols</td><td>SQMG</td></tr>
<tr><td>Sulfoquinovosyl diacylglycerol</td><td>SQDG</td></tr>
<tr><td>Triacylglycerol</td><td>TAG</td></tr>


<tr><td rowspan="3">Saccharo lipids (SL)</td><td></td><td>DAT</td></tr>
<tr><td></td><td>AC2SGL</td></tr>
<tr><td></td><td>PAT16 / PAT18</td></tr>


<tr><td rowspan="61">Mediator (LM)</td><td rowspan="9">Docosanoids</td><td>10-HDoHE</td></tr>
<tr><td>11-HDoHE</td></tr>
<tr><td>16-HDoHE</td></tr>
<tr><td>8-HDoHE</td></tr>
<tr><td>Maresin 1</td></tr>
<tr><td>Resolvin D1</td></tr>
<tr><td>Resolvin D2</td></tr>
<tr><td>Resolvin D3</td></tr>
<tr><td>Resolvin D5</td></tr>

<tr><td rowspan="39">Docosanoids</td><td>11(12)-EET</td></tr>
<tr><td>11,12-DHET</td></tr>
<tr><td>11-HETE</td></tr>
<tr><td>12-HEPE</td></tr>
<tr><td>12-HETE</td></tr>
<tr><td>12-HHTrE</td></tr>
<tr><td>12-OxoETE</td></tr>
<tr><td>14(15)-EET</td></tr>
<tr><td>14(15)-EpETE</td></tr>
<tr><td>14,15-DHET</td></tr>
<tr><td>15d-PGJ2</td></tr>
<tr><td>15-HEPE</td></tr>
<tr><td>15-HETE}</td></tr>
<tr><td>16-HETE</td></tr>
<tr><td>18-HEPE</td></tr>
<tr><td>5(6)-EET</td></tr>
<tr><td>5,12-DiHETE</td></tr>
<tr><td>5,6,15-LXA4</td></tr>
<tr><td>5,6-DiHETE</td></tr>
<tr><td>5-HEPE</td></tr>
<tr><td>5-HETE</td></tr>
<tr><td>5-HpETE</td></tr>
<tr><td>5-OxoETE</td></tr>
<tr><td>8(9)-EET</td></tr>
<tr><td>8,9-DHET</td></tr>
<tr><td>8-HETE</td></tr>
<tr><td>9-HEPE</td></tr>
<tr><td>9-HETE</td></tr>
<tr><td>LTB4</td></tr>
<tr><td>LTC4</td></tr>
<tr><td>LTD4</td></tr>
<tr><td>PGB2</td></tr>
<tr><td>PGD2</td></tr>
<tr><td>PGE2</td></tr>
<tr><td>PGF2alpha</td></tr>
<tr><td>PGI2</td></tr>
<tr><td>TXB1</td></tr>
<tr><td>TXB2</td></tr>
<tr><td>TXB3</td></tr>

<tr><td rowspan="6">Octadecanoids</td><td>12(13)-EpOME</td></tr>
<tr><td>13-HODE</td></tr>
<tr><td>13-HOTrE</td></tr>
<tr><td>9(10)-EpOME</td></tr>
<tr><td>9-HODE</td></tr>
<tr><td>9-HOTrE</td></tr>

<tr><td rowspan="7">Fatty Acids and Conjugates</td><td>AA (Arachidonic acid)</td></tr>
<tr><td>ALA (α-Linolenic acid)</td></tr>
<tr><td>DHA (Docosahexaenoic acid)</td></tr>
<tr><td>EPA (Eicosapentaenoic acid)</td></tr>
<tr><td>Linoleic acid</td></tr>
<tr><td>Palmitic acid</td></tr>
<tr><td>Tetranor-12-HETE</td></tr>

</table>
