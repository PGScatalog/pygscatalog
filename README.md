# pygscatalog (namespace package test)

Test set up of a library subpackage and a CLI subpackage. 

The CLI subpackage imports a public object from the library subpackage and prints it. 

Both the CLI package and library package are in the namespace `pgscatalog`.

```
$ cd downloadapp
$ poetry shell
$ poetry install
$ poetry run pgscatalog-download
I imported a test object <pgscatalog.corelib.scoringfile.ScoringFile object at 0x104e561d0>
```

**TODO: look at importing subpackages in root pyproject.toml**
