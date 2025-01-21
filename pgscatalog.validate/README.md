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

## Installation 

Install and run from source code via poetry and pip:

```
# Download the source code
git clone https://github.com/PGScatalog/pygscatalog.git
cd pygscatalog/pygscatalog.validate

# Install the dependencies
poetry install
```

### Run with a virtual environment

From the directory where `poetry install` was executed:

```
poetry shell
pgscatalog-validate --help
```

### Or install the package in the current environment

```
# Build the wheel package
poetry build

# Install the built package in the current environment via pip
pip install dist/pgscatalog_validate-0.1-py3-none-any.whl

# Run
pgscatalog-validate --help
```

## Documentation

Documentation is available at https://pygscatalog.readthedocs.io/.