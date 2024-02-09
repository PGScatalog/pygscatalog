# `pgscatalog-combine`

[![.github/workflows/combineapp-pytest.yml](https://github.com/PGScatalog/pygscatalog/actions/workflows/combineapp-pytest.yml/badge.svg)](https://github.com/PGScatalog/pygscatalog/actions/workflows/combineapp-pytest.yml)

A CLI application that combines polygenic score (PGS) scoring files into a standardised output. 

The process involves:

* Normalising scoring files to a consistent data format
* Doing some quality control checks on variants
* Writing a long format / melted file 

Input scoring files must follow [PGS Catalog standards](https://www.pgscatalog.org/downloads/#dl_scoring_files). The output file is useful for doing data science tasks, like matching variants across a scoring file and target genome. 

The [PGS Catalog Calculator](https://github.com/PGScatalog/pgsc_calc) uses this CLI application internally.

## Installation

```
$ pip install pgscatalog-combine
```

## Usage

```
$ pgscatalog-combine -s PGS000001_hmPOS_GRCh38.txt.gz PGS0001229_hmPOS_GRCh38.txt.gz -t GRCh38 -o combined.txt
```

## Documentation

See **link to docs**.

## Help

```
$ pgscatalog-combine --help
```
