# pygscatalog

This repository contains Python applications and libraries for working with polygenic scores (PGS) and the [PGS Catalog](https://www.pgscatalog.org/), an open database of polygenic scores and the relevant metadata required for accurate application and evaluation. 
## User applications 

These CLI applications are used by the PGS Catalog Calculator workflow. 

| Application           | Description                                    | Link                                                 |
|-----------------------|------------------------------------------------|------------------------------------------------------|
| `pgscatalog-download` | Download scoring files from the PGS Catalog    | [README](pgscatalog.downloadapp/README.md) |
| `pgscatalog-combine`  | Combine scoring files into a consistent format |


## Developer libraries

If you write  code to work with PGS, we publish some libraries that might be helpful:


| Library                | Description                                              | Link                                              |
|------------------------|----------------------------------------------------------|---------------------------------------------------|
| `pgscatalog-corelib`   | Core classes and functions                               | [README](pgscatalog.corelib/README.md) |
| `pgscatalog-matchlib`  | Variant matching across scoring files and target genomes |
| `pgscatalog-calclib`   | Ancestry estimation and normalisation                    |


## Installation

### pip

If you want to use the packages in this repository, use pip:

### Local install for developers

If you want to make changes to a package or application, it's simplest to clone the repository and install packages in editable mode.

```
$ git clone https://github.com/PGScatalog/pygscatalog.git
$ cd pygscatalog/pgscatalog.downloadapp # replace with the package you want to edit
$ poetry add --editable ../pgscatalog.corelib # downloadapp requires corelib
$ poetry install  
```

## Documentation

Full documentation is provided.. 

## License

All of our code is open source and licensed with [Apache 2](LICENSE).
