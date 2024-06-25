# `pgscatalog-utils`

This convenience package bundles every PGS Catalog application, including:

| Application                  | Description                                                            | Link                                                                                |
|------------------------------|------------------------------------------------------------------------|-------------------------------------------------------------------------------------|
| `pgscatalog-download`        | Download scoring files from the PGS Catalog in specific genome builds  | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/download.html)  |
| `pgscatalog-combine`         | Combine multiple scoring files into a consistent structure             | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/combine.html)   |
| `pgscatalog-relabel`         | Relabel values in a column based on values in a column in another file | [README](pgscatalog.core/README.md)                                                 |
| `pgscatalog-match`           | Match structured scoring file to variants in target genomes            | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/match.html)     |
| `pgscatalog-matchmerge`      | Merge variant match results, useful on larger datasets                 | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/match.html)     |
| `pgscatalog-intersect`       | Match variants across two different variant information files (e.g. reference & target genomes) | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/intersect.html) |
| `pgscatalog-aggregate`       | Aggregate calculated PGS split across multiple files                   | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/aggregate.html) |
| `pgscatalog-ancestry-adjust` | Adjust calculated PGS in the context of genetic ancestry               | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/ancestry.html)  |

Please note `v1.0.0` contains breaking changes: CLI applications have been renamed and the package has been significantly refactored.

See https://github.com/PGScatalog/pygscatalog for more details.

## Installation 

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pgscatalog-utils/README.html)

or [install via `pip`](https://pypi.org/project/pgscatalog-utils/):

```
$ pipx install pgscatalog-utils
```

## Documentation 

Documentation is available at https://pygscatalog.readthedocs.io/.
