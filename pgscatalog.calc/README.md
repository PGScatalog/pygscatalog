# `pgscatalog-calc`

[![Run pytest on pgscatalog.calc](https://github.com/PGScatalog/pygscatalog/actions/workflows/calc-pytest.yml/badge.svg)](https://github.com/PGScatalog/pygscatalog/actions/workflows/calc-pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/pygscatalog/badge/?version=latest)](https://pygscatalog.readthedocs.io/en/latest/autoapi/pgscatalog/calc/index.html)

This Python package contains:

* CLI applications to aggregate and adjust calculated PGS in the context of genetic ancestry similarity
* library classes and functions for working with calculated PGS and PCA data

| Application                  | Description                                              | Link                                                                                |
|------------------------------|----------------------------------------------------------|-------------------------------------------------------------------------------------|
| `pgscatalog-aggregate`       | Aggregate calculated PGS split across multiple files     | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/aggregate.html) |
| `pgscatalog-ancestry-adjust` | Adjust calculated PGS in the context of genetic ancestry | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/ancestry.html)  |

If you want to write Python code to work with PGS data, the library may be helpful for you.

Please note this package doesn't contain functionality to calculate a PGS from target genomes and scoring files. Instead, it's focused on aggregation and normalisation. If you want to do automatic PGS calculation check out the [PGS Catalog Calculator](https://github.com/PGScatalog/pgsc_calc), which uses these tools internally.

## Installation 

```
$ pipx install pgscatalog-calc
```

## Documentation

Documentation is available at https://pygscatalog.readthedocs.io/.
