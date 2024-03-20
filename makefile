main:
	echo ""

install:
	python3 setup.py install
	
test:
	python3 -m unittest pygoslin.tests.Parser_test
	python3 -m unittest pygoslin.tests.Formulas_test
	python3 -m unittest pygoslin.tests.Masses_test
	python3 -m unittest pygoslin.tests.FattyAcids_test
	python3 -m unittest pygoslin.tests.SwissLipids_test
	python3 -m unittest pygoslin.tests.Shorthand_test
	python3 -m unittest pygoslin.tests.Goslin_test
	python3 -m unittest pygoslin.tests.LipidMaps_test
	python3 -m unittest pygoslin.tests.Hmdb_test

distclean:
	rm -rf build dist pygoslin.egg-info .eggs
	
build:
	python3 -m build --sdist
	
uploadtest:
	python3 -m twine upload --repository testpypi dist/* -u __token__ -p PERSONAL_TOKEN
	
upload:
	python3 -m twine upload dist/* -u __token__ -p PERSONAL_TOKEN

# upload to pypi:
# distclean
# build
# uploadtest
# upload
