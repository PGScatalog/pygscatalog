# `pgscatalog.calc`

[![Run CI on pgscatalog.calc](https://github.com/PGScatalog/pygscatalog/actions/workflows/calc-ci.yml/badge.svg)](https://github.com/PGScatalog/pygscatalog/actions/workflows/calc-ci.yml)
[![Documentation Status](https://readthedocs.org/projects/pygscatalog/badge/?version=latest)](https://pygscatalog.readthedocs.io/en/latest/autoapi/pgscatalog/calc/index.html)

> [!IMPORTANT]  
> * In `1.0.0-alpha.1` programs to calculate polygenic scores were added to this package
> * The `pgsc_calc load` and `pgsc_calc score` programs are **pre-release software** that support basic PGS calculation only (chromosomes 1-22, indexed VCF and BGEN files, no ancestry adjustment)
> * The simplest way to calculate polygenic scores with this software is to use v3 of the [PGS Catalog Calculator](https://github.com/PGScatalog/pgsc_calc)
> * Documentation for this package and more PGS calculation features will be added to future releases

This Python package contains:

* CLI applications to aggregate and adjust calculated polygenic scores (PGS) in the context of genetic ancestry similarity
* library classes and functions for working with calculated PGS and PCA data


| Application                  | Description                                              | Link                                                                                |
|------------------------------|----------------------------------------------------------|-------------------------------------------------------------------------------------|
| `pgsc_calc load`             | Query an indexed VCF/BGEN and create a zarr zip archive  | `pgsc_calc load --help`                                                             |
| `pgsc_calc score`            | Calculate polygenic scores from zarr zip archives        | `pgsc_calc score --help`                                                            |
| `pgscatalog-aggregate`       | Aggregate calculated PGS split across multiple files     | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/aggregate.html) |
| `pgscatalog-ancestry-adjust` | Adjust calculated PGS in the context of genetic ancestry | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/ancestry.html)  |

If you want to write Python code to work with PGS data, the library may be helpful for you.

## Installation 

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pgscatalog.calc/README.html)

or [install via `pip`](https://pypi.org/project/pgscatalog.calc/):

```
$ pipx install pgscatalog.calc
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
$ pgscatalog-aggregate --help
```

`nox` can also be used to run tests and lint the package:

```
$ nox
```
