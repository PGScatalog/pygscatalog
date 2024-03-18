# `pgscatalog-calcapp`

[![.github/workflows/calcapp-pytest.yml](https://github.com/PGScatalog/pygscatalog/actions/workflows/calcapp-pytest.yml/badge.svg?branch=main)](https://github.com/PGScatalog/pygscatalog/actions/workflows/calcapp-pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/pygscatalog/badge/?version=latest)](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/aggregate.html)

A CLI application that can:

* aggregate PGS calculated with `plink2 --score` that are split across multiple files (e.g. per chromosome)
* estimate [genetic ancestry similarity and adjust calculated PGS](https://pgsc-calc.readthedocs.io/en/latest/explanation/geneticancestry.html)

The [PGS Catalog Calculator](https://github.com/PGScatalog/pgsc_calc) uses this CLI application internally. 


## Installation 

```
$ pip install pgscatalog-calcapp
```
