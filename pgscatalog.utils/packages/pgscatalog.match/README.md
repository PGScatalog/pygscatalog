# pgscatalog.match

[![.github/workflows/match-pytest.yml](https://github.com/PGScatalog/pygscatalog/actions/workflows/match-pytest.yml/badge.svg?branch=main)](https://github.com/PGScatalog/pygscatalog/actions/workflows/match-pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/pygscatalog/badge/?version=latest)](https://pygscatalog.readthedocs.io/en/latest/autoapi/pgscatalog/match/index.html)

This Python package contains:

* CLI applications to match variants from PGS Catalog scoring files against `plink` variant information files and output files ready for `plink2 --score`
* core library classes and functions for matching variants and creating plink scoring files

If you want to write Python code to match genetic variants, the library may be helpful for you.

If you want to do automatic PGS calculation check out the [PGS Catalog Calculator](https://github.com/PGScatalog/pgsc_calc), which uses these tools internally.

## Installation 

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pgscatalog.match/README.html)

or [install via `pip`](https://pypi.org/project/pgscatalog.match/):

```
$ pipx install pgscatalog.match
```

## Documentation

Documentation is available at https://pygscatalog.readthedocs.io/.

## Developer instructions

You'll need [`nox`](https://nox.thea.codes/en/stable/index.html) and [`uv`](https://github.com/astral-sh/uv) installed. 

To get set up with a development environment run:

```
$ nox -s dev
$ source .venv/bin/activate
```

This will create a virtual environment in the current directory.

```
$ pgscatalog-match --help
```

`nox` can also be used to run tests and lint the package:

```
$ nox
```