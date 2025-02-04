# `pgscatalog.validate`

This Python package contains:

* CLI applications to check/validate that the scoring files and harmonized scoring files match the PGS Catalog scoring file formats
* library classes and functions for working with scoring file validation

| Application           | Description            | Link                                                                               |
|-----------------------|------------------------|------------------------------------------------------------------------------------|
| `pgscatalog-validate` | Validate scoring files | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/validate.html) |

If you want to write Python code to work with PGS data, the library may be helpful for you.

Please note that this tool validates formatted scoring files which are expected to contain the PGS-Catalog standard header.
This header is not necessary for new scoring files submission from authors.

## Developer instructions

You'll need [`nox`](https://nox.thea.codes/en/stable/index.html) and [`uv`](https://github.com/astral-sh/uv) installed. 

To get set up with a development environment run:

```
# Download the source code
$ git clone https://github.com/PGScatalog/pygscatalog.git
$ cd pygscatalog/pgscatalog.validate
$ nox -s dev
$ source .venv/bin/activate
```
This creates a virtual environment in the same directory.

```
$ pgscatalog-validate --help
```

### Or install the package in the current environment

```
# Build the wheel package
$ nox -s build

# Install the built package in the current environment via pip
$ pip install dist/pgscatalog_validate-0.1-py3-none-any.whl

# Run
$ pgscatalog-validate --help
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
$ pgscatalog-download --help
```
