# `pgscatalog-download`

[![.github/workflows/downloadapp-pytest.yml](https://github.com/PGScatalog/pygscatalog/actions/workflows/downloadapp-pytest.yml/badge.svg)](https://github.com/PGScatalog/pygscatalog/actions/workflows/downloadapp-pytest.yml)

`pgscatalog-download` is a CLI application that makes it easy to download scoring files from the
PGS Catalog with a mixture of PGS, publication, or trait accessions.

The [PGS Catalog Calculator](https://github.com/PGScatalog/pgsc_calc) uses this CLI application internally.

## Installation

```
$ pip install pgscatalog-download
```

## Usage 

```
$ mkdir scorefiles
$ pgscatalog-download --pgs PGS000822 --pgp PGP000001 PGP000002 --efo EFO_0004214 -b GRCh38 -o scorefiles/
$ tree scorefiles
scorefiles
├── PGS000001_hmPOS_GRCh38.txt.gz
├── PGS000002_hmPOS_GRCh38.txt.gz
├── PGS000003_hmPOS_GRCh38.txt.gz
...
```

## Documentation

See **link to docs**.

## Help

```
$ pgscatalog-download --help
```