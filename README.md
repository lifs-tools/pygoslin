# Python3 implementation for Goslin

This is the mzTab-M reference implementation and validation API service for Python 3.

> **_NOTE:_**  This is an *early* development version, please use at your own risk and report issues to help improve it!

mzTab-M is intended as a reporting standard for quantitative results from metabolomics/lipodomics approaches. This format is further intended to provide local LIMS systems as well as MS metabolomics repositories a simple way to share and combine basic information.

mzTab-M has been developed with a view to support the following general tasks:

1. Facilitate the sharing of final experimental results, especially with researchers outside the field of metabolomics.
2. Export of results to external software, including programs such as Microsoft Excel® and Open Office Spreadsheet and statistical software / coding languages such as R.
3. Act as an output format of (web-) services that report MS-based results and thus can produce standardized result pages.
4. Be able to link to the external experimental evidence e.g. by referencing back to mzML files.




## Installation

### Prerequisites
The pymzTab-m package uses pipenv to create an isolated and defined build environment. You will need Python >=3.5 to build the package.

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

The two major function within pygoslin are the parsing and printing of lipid names. This example will demonstrate both functions. Open a python shell and type in:

```
from pygoslin.parser.Parser import Parser
from pygoslin.parser.GoslinParserEventHandler import GoslinParserEventHandler
goslin_parser_event_handler = GoslinParserEventHandler()
goslin_parser = Parser(goslin_parser_event_handler, "data/goslin/Goslin.g4")

lipid_name = "PE 16:1/12:0"
goslin_parser.parse(lipid_name)
goslin_parser.raise_events()
lipid = goslin_parser_event_handler.lipid

print(lipid.get_lipid_string())
```


