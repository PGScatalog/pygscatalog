# `pgscatalog.validate`

This Python package contains:

* CLI applications to check/validate that the scoring files and harmonized scoring files match the PGS Catalog scoring file formats
* library classes and functions for working with scoring file validation

| Application           | Description            | Link                                                                               |
|-----------------------|------------------------|------------------------------------------------------------------------------------|
| `pgscatalog-validate` | Validate scoring files | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/validate.html) |

If you want to write Python code to work with PGS data, the library may be helpful for you.

## Installation 

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pgscatalog.core/README.html)

or [install via `pip`](https://pypi.org/project/pgscatalog.core/):

```
$ pipx install pgscatalog.validate
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
$ pgscatalog-validate --help
```

`nox` can also be used to run tests and lint the package:

```
$ nox
```