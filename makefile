main:
	echo ""

install:
	python3 setup.py install
	
test:
	python3 -m unittest pygoslin.tests.FattyAcidTest
	python3 -m unittest pygoslin.tests.ParserTest
	python3 -m unittest pygoslin.tests.TestFormulas
	python3 -m unittest pygoslin.tests.MassesTest
	python3 -m unittest pygoslin.tests.SwissLipidsTest
	python3 -m unittest pygoslin.tests.GoslinTest
	python3 -m unittest pygoslin.tests.LipidMapsTest
	python3 -m unittest pygoslin.tests.HmdbTest
