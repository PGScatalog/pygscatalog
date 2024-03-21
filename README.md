# pygscatalog

[![CI](https://github.com/PGScatalog/pygscatalog/actions/workflows/pytest.yaml/badge.svg)](https://github.com/PGScatalog/pygscatalog/actions/workflows/pytest.yaml)
[![codecov](https://codecov.io/github/PGScatalog/pygscatalog/graph/badge.svg?token=EEAU59C8IK)](https://codecov.io/github/PGScatalog/pygscatalog)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)

This repository contains Python applications and libraries for working with polygenic scores (PGS :dna:) and the 
[PGS Catalog](https://www.pgscatalog.org/), an open database of polygenic scores and the relevant metadata required for 
accurate application and evaluation. It is based on a previous codebase of utilities ([pgscatalog_utils](https://github.com/PGScatalog/pgscatalog_utils)) 
that has been converted to namespace packages for modularity and re-use.

## User applications 

These CLI applications are used internally by the [PGS Catalog Calculator (`pgsc_calc`)](https://github.com/PGScatalog/pgsc_calc) 
workflow for calculating PGS and performing common adjustments for genetic ancestry. 

If you want an automatic method of calculating PGS, including genetic ancestry similarity estimation and PGS normalisation, 
the workflow is the easiest method.


| Application                  | Description                                                            | Install                         | Link                                 |
|------------------------------|------------------------------------------------------------------------|---------------------------------|--------------------------------------|
| `pgscatalog-download`        | Download scoring files from the PGS Catalog in specific genome builds  | `pipx install pgscatalog-core`  | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/download.html)  |
| `pgscatalog-combine`         | Combine multiple scoring files into a consistent structure             | `pipx install pgscatalog-core`  | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/combine.html)  |
| `pgscatalog-relabel`         | Relabel values in a column based on values in a column in another file | `pipx install pgscatalog-core`  | [README](pgscatalog.core/README.md)  |
| `pgscatalog-match`           | Match structured scoring file to variants in target genomes            | `pipx install pgscatalog-match` | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/match.html) |
| `pgscatalog-matchmerge`      | Merge variant match results, useful on larger datasets                 | `pipx install pgscatalog-match` | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/match.html) |
| `pgscatalog-aggregate`       | Aggregate calculated PGS split across multiple files                   | `pipx install pgscatalog-calc`  | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/aggregate.html)  |
| `pgscatalog-ancestry-adjust` | Adjust calculated PGS in the context of genetic ancestry               | `pipx install pgscatalog-calc`  | [README](https://pygscatalog.readthedocs.io/en/latest/how-to/guides/ancestry.html)  |


## Developer libraries

If you write Python code to work with PGS, the underlying libraries for the apps are documented and available for re-use:

| Library            | Description                                              | Link                                                                                              |
|--------------------|----------------------------------------------------------|---------------------------------------------------------------------------------------------------|
| `pgscatalog-core`  | Core classes and functions to work with PGS data         | [API reference](https://pygscatalog.readthedocs.io/en/latest/autoapi/pgscatalog/core/index.html)  |
| `pgscatalog-match` | Variant matching across scoring files and target genomes | [API reference](https://pygscatalog.readthedocs.io/en/latest/autoapi/pgscatalog/match/index.html) |
| `pgscatalog-calc`  | Genetic ancestry similarity estimation and normalisation | [API reference](https://pygscatalog.readthedocs.io/en/latest/autoapi/pgscatalog/calc/index.html)  |

## Documentation

Full documentation for the applications and libraries is available at [https://pygscatalog.readthedocs.io/](https://pygscatalog.readthedocs.io/en/latest/).

## Credits & Licence

`pygscatalog`(_aka_ `pgscatalog_utils`) is developed as part of the [PGS Catalog](https://www.pgscatalog.org/) project, a
collaboration between the University of Cambridge’s Department of Public Health and Primary Care (Michael Inouye, 
Samuel Lambert) and the European Bioinformatics Institute (Helen Parkinson, Laura Harris).

This package contains code libraries and apps for working with PGS Catalog data and calculating PGS within the 
[PGS Catalog Calculator (`pgsc_calc`)](https://github.com/PGScatalog/pgsc_calc) workflow, and is based on an earlier 
codebase ([pgscatalog_utils](https://github.com/PGScatalog/pgscatalog_utils)) with contributions and input from members 
of the PGS Catalog team (Samuel Lambert, Benjamin Wingfield, Aoife McMahon Laurent Gil) and Inouye lab 
(Rodrigo Canovas, Scott Ritchie, Jingqin Wu).

A manuscript describing this package and `pgsc_calc` pipeline is *in preparation*. In the meantime if you use the
tool we ask you to cite the repo and the paper describing the PGS Catalog resource:

- >PGS Catalog Calculator _(in preparation)_. PGS Catalog
  Team. [https://github.com/PGScatalog/pgsc_calc](https://github.com/PGScatalog/pgsc_calc)
- >Lambert _et al._ (2021) The Polygenic Score Catalog as an open database for
reproducibility and systematic evaluation.  Nature Genetics. 53:420–425
doi:[10.1038/s41588-021-00783-5](https://doi.org/10.1038/s41588-021-00783-5).

All of our code is open source and permissively licensed with [Apache 2](LICENSE).

This work has received funding from EMBL-EBI core funds, the Baker Institute, the University of Cambridge, 
Health Data Research UK (HDRUK), and the European  Union’s Horizon 2020 research and innovation programme under grant 
agreement No 101016775 INTERVENE.
