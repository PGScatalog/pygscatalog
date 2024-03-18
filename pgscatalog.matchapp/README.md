# `pgscatalog-match`

[![.github/workflows/matchapp-pytest.yml](https://github.com/PGScatalog/pygscatalog/actions/workflows/matchapp-pytest.yml/badge.svg?branch=main)](https://github.com/PGScatalog/pygscatalog/actions/workflows/matchapp-pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/pygscatalog/badge/?version=latest)](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/match.html)

A CLI application that matches (intersects) variants in a normalised scoring file against variant information files. It outputs scoring files compatible with `plink2 --score`, a summary log, and a large match candidate information log.

Various match strategies are used to generate a list of match candidates. These candidates are filtered by user parameters to produce up to one "best match"  variant for each variant in the original scoring file.

Matching fails if not enough variants are present in the target variant information files. You can adjust `--min_overlap` to change the matching failure threshold, but it's not a good idea normally.

The [PGS Catalog Calculator](https://github.com/PGScatalog/pgsc_calc) uses this CLI application internally. 

## Installation 

```
$ pip install pgscatalog-match
```
