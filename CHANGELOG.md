## 2024-06-13

Minor release [pgscatalog-match-0.2.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.match-0.2.0)

* Improve performance when reading target genomes and scoring files https://github.com/PGScatalog/pygscatalog/issues/20
* Fix raising MatchRateError when one or more scores fails matching but at least one score passes https://github.com/PGScatalog/pygscatalog/issues/21
* Fix warnings about map_element dtypes during labelling 

Minor release [pgscatalog-calc-0.2.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.calc-0.2.0)

* Improve aggregation to have clearer logic and logging https://github.com/PGScatalog/pygscatalog/pull/23 
  * Fix aggregating scores with different numbers of columns
* Fix adjusting averages during ancestry adjustment https://github.com/PGScatalog/pygscatalog/pull/19

## 2024-06-12

Minor release [pgscatalog-core-0.2.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.core-0.2.0)

* Drop pyarrow support when reading target genomes or scoring files 
* Add natural sorting functions to sort scoring file paths by effect type or chromosome
* Remove pyarrow dependency 
  
## 2024-05-23

Patch release [pgscatalog-core-0.1.2](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.core-0.1.2) and [pgscatalog.calc-0.1.1](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.calc-0.1.1):

* Bug fixes for PGS Catalog Calculator integration

Patch release [pgscatalog.match-0.1.1](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.match-0.1.1):

* Bug fixes for PGS Catalog Calculator integration
* New pgscatalog-intersect CLI application

Bump [pgscatalog-utils](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog-utils-1.0.2) to use patched dependencies

## 2024-04-09

Patch release [pgscatalog-core-0.1.1](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.core-0.1.1):

* Add NCBI35 to GenomeBuild
* Disable concurrency in `pgscatalog-combine` to prevent memory usage explosion (will be fixed and re-enabled later)
* Fix handling scoring file cases when only one optional effect type column is present
* Fix reading pvar files with complex headers (like ones made from a VCF)

Bump [pgscatalog-utils](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog-utils-1.0.1) to use patched core

## 2024-03-22

Initial release of new python packages:

* [pgscatalog.core-0.1.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.core-0.1.0)
* [pgscatalog.match-0.1.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.match-0.1.0)
* [pgscatalog.calc-0.1.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.calc-0.1.0)

Also [pgscatalog-utils](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog-utils-1.0.0) migrated to using these tools on the backend (v1.0.0)
