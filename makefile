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
	python3 -m build
	
uploadtest:
	twine upload --repository testpypi dist/* -u __token__ -p PERSONAL_TOKEN
	
upload:
	twine upload dist/* -u __token__ -p PERSONAL_TOKEN

# upload to pypi:
## system python3 must be enabled
# python3 -m venv venv
# source venv/bin/activate
# pip3 install twine build setuptools
# pip3 install --upgrade pip setuptools wheel build
# rm -rf build dist pygoslin.egg-info .eggs   (make distclean)
# python3 -m build  (make build)
# twine upload --repository testpypi dist/* -u __token__ -p PERSONAL_TOKEN  (make uploadtest)
# twine upload dist/* -u __token__ -p PERSONAL_TOKEN   (make upload)
