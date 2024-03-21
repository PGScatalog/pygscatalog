# `pgscatalog-core`

[![Run pytest on pgscatalog.core](https://github.com/PGScatalog/pygscatalog/actions/workflows/core-pytest.yml/badge.svg)](https://github.com/PGScatalog/pygscatalog/actions/workflows/core-pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/pygscatalog/badge/?version=latest)](https://pygscatalog.readthedocs.io/en/latest/autoapi/pgscatalog/core/index.html)

This Python package contains:

* CLI applications to download scoring files from the PGS Catalog and normalise them into a consistent structure
* core library classes and functions for working with PGS data

| Application           | Description                                                            | Link                                                                               |
|-----------------------|------------------------------------------------------------------------|------------------------------------------------------------------------------------|
| `pgscatalog-download` | Download scoring files from the PGS Catalog in specific genome builds  | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/download.html) |
| `pgscatalog-combine`  | Combine multiple scoring files into a consistent structure             | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/combine.html)  |
| `pgscatalog-relabel`  | Relabel values in a column based on values in a column in another file | [README](README.md)                                                |

If you want to write Python code to work with PGS data, the library may be helpful for you.

If you want to do automatic PGS calculation check out the [PGS Catalog Calculator](https://github.com/PGScatalog/pgsc_calc), which uses these tools internally.

## Installation 

```
$ pipx install pgscatalog-core
```

## Documentation

Documentation is available at https://pygscatalog.readthedocs.io/.

