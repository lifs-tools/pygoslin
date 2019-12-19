main:
	echo ""

install:
	python3 setup.py install
	
test:
	python3 -m unittest pygoslin.tests.MolecularFattyAcidTest
	python3 -m unittest pygoslin.tests.ParserTest