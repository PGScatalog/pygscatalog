# pgscatalog-match

[![.github/workflows/match-pytest.yml](https://github.com/PGScatalog/pygscatalog/actions/workflows/match-pytest.yml/badge.svg?branch=main)](https://github.com/PGScatalog/pygscatalog/actions/workflows/match-pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/pygscatalog/badge/?version=latest)](https://pygscatalog.readthedocs.io/en/latest/autoapi/pgscatalog/match/index.html)

This Python package contains:

* CLI applications to match variants from PGS Catalog scoring files against `plink` variant information files and output files ready for `plink2 --score`
* core library classes and functions for matching variants and creating plink scoring files

If you want to write Python code to match genetic variants, the library may be helpful for you.

If you want to do automatic PGS calculation check out the [PGS Catalog Calculator](https://github.com/PGScatalog/pgsc_calc), which uses these tools internally.

## Installation 

```
$ pipx install pgscatalog-match
```

## Documentation

Documentation is available at https://pygscatalog.readthedocs.io/.

