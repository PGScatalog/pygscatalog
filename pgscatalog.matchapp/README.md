# `pgscatalog-match`

A CLI application that matches (intersects) variants in a normalised scoring file against variant information files. It outputs scoring files compatible with `plink2 --score`, a summary log, and a large match candidate information log.

Various match strategies are used to generate a list of match candidates. These candidates are filtered by user parameters to produce up to one "best match"  variant for each variant in the original scoring file.

Matching fails if not enough variants are present in the target variant information files. You can adjust `--min_overlap` to change the matching failure threshold, but it's not a good idea normally.

The [PGS Catalog Calculator](https://github.com/PGScatalog/pgsc_calc) uses this CLI application internally. If you want to calculate polygenic scores using an automated workflow, the PGS Catalog Calculator is the best method to use.

## Installation 

```
$ pip install pgscatalog-match
```

## Usage

```
$ pgscatalog-match ...
```

## Documentation

See **link to docs**.

## Help

```
$ pgscatalog-match --help
```