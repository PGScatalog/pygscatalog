# pygscatalog (namespace package test)

Test set up of a library subpackage and a CLI subpackage. 

The CLI subpackage imports a public object from the library subpackage and prints it. 

Both the CLI package and library package are in the namespace `pgscatalog`.

```
$ cd pgscatalog.downloadapp
$ poetry shell
$ poetry install
$ poetry run pgscatalog-download
I imported a test object <pgscatalog.corelib.scoringfile.ScoringFile object at 0x100ab6310>
I imported a test object <pgscatalog.calclib.testclass.TestClass object at 0x100a464d0> from a different package
```

[Helpful naming conventions](https://peps.python.org/pep-0423/)
