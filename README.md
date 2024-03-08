# pygscatalog

This repository contains Python applications and libraries for working with polygenic scores (PGS :dna:) and the [PGS Catalog](https://www.pgscatalog.org/), an open database of polygenic scores and the relevant metadata required for accurate application and evaluation.

## User applications 

These CLI applications are used internally by the [PGS Catalog Calculator](https://github.com/PGScatalog/pgsc_calc) workflow. 

If you want an automatic method of calculating PGS, including genetic ancestry similarity estimation and PGS normalisation, the workflow is the easiest method.


| Application                  | Description                                                | Link                                       |
|------------------------------|------------------------------------------------------------|--------------------------------------------|
| `pgscatalog-download`        | Download scoring files from the PGS Catalog                | [README](pgscatalog.downloadapp/README.md) |
| `pgscatalog-combine`         | Combine multiple scoring files into a consistent structure | [README](pgscatalog.combineapp/README.md) |
| `pgscatalog-match`           | Match scoring file variants to target genomes              | [README](pgscatalog.matchapp/README.md) |
| `pgscatalog-matchmerge`      | Merge variant match results, useful on larger datasets     | [README](pgscatalog.matchapp/README.md) |
| `pgscatalog-aggregate`       | Aggregate calculated PGS split across multiple files       | [README](pgscatalog.calcapp/README.md) |
| `pgscatalog-ancestry-adjust` | Adjust calculated PGS in the context of genetic ancestry   | [README](pgscatalog.calcapp/README.md) |

## Developer libraries

If you write Python code to work with PGS, we publish some libraries that might be helpful:

| Library               | Description                                              | Link                                   |
|-----------------------|----------------------------------------------------------|----------------------------------------|
| `pgscatalog-corelib`  | Core classes and functions to work with PGS data         | [README](pgscatalog.corelib/README.md) |
| `pgscatalog-matchlib` | Variant matching across scoring files and target genomes | [README](pgscatalog.corelib/README.md) |
| `pgscatalog-calclib`  | Genetic ancestry similarity estimation and normalisation | [README](pgscatalog.corelib/README.md) |

## Documentation

Full documentation is provided **LINK TO DOCS**

## License

All of our code is open source and permissively licensed with [Apache 2](LICENSE).
