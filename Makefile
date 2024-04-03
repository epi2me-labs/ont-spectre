.PHONY: develop docs

PYTHON ?= python3

IN_VENV=. ./venv/bin/activate
PROJECT=spectre

venv/bin/activate:
	test -d venv || $(PYTHON) -m venv venv
	${IN_VENV} && pip install pip --upgrade
	${IN_VENV} && pip install -r requirements.txt

develop: venv/bin/activate
	${IN_VENV} && python setup.py develop

test: venv/bin/activate
	# we typically have more flake8 extensions here, but the project is in such a state...
	${IN_VENV} && pip install pytest pytest-cov
	${IN_VENV} && pytest tests --cov spectre \
		--cov-report html --cov-report term --cov-report term-missing

IN_BUILD=. ./pypi_build/bin/activate
pypi_build/bin/activate:
	test -d pypi_build || $(PYTHON) -m venv pypi_build --prompt "(pypi) "
	${IN_BUILD} && pip install pip --upgrade
	${IN_BUILD} && pip install --upgrade pip setuptools twine wheel readme_renderer[md] keyrings.alt

.PHONY: sdist
sdist: pypi_build/bin/activate
	${IN_BUILD} && python setup.py sdist

.PHONY: clean
clean:
	rm -rf __pycache__ dist build venv ${PROJECT}.egg-info tmp docs/_build
