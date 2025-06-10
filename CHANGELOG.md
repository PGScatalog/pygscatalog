## 2025-06-10

Major release [pgscatalog-core-1.0.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog-core-1.0.0)

* Add multiprocessing support to `pgscatalog-format` (`pgscatalog-combine` is now deprecated), each input scoring file now gets written to a separate output file
* Fix `pgscatalog-download` sometimes leaving temporary files in the output directory

Minor release [pgscatalog-match-0.4.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog-match-0.4.0)

* Add support for multiple input scoring files to `pgscatalog-match` and `pgscatalog-matchcombine`

Major release [pgscatalog-utils-2.0.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog-utils-2.0.0)

* Bump dependencies

## 2024-10-29

Patch release [pgscatalog-utils-1.4.2](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog-utils-1.4.2)

* Bump core dependency

Patch release [pgscatalog-core-0.3.3](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.core-0.3.3)

* Remove rsID format check from `CatalogScoreVariant` model

## 2024-10-08

Patch release [pgscatalog-utils-1.4.1](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog-utils-1.4.1)

* Fix missing imports
* Bump dependencies

Patch release [pgscatalog-match-0.3.3](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.match-0.3.3)

* Migrate classes deleted in `pgscatalog.core` to `pgscatalog.match`
* Add missing test data

Patch release [pgscatalog-utils-1.4.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog-utils-1.4.0)

* Bump dependencies
 
## 2024-10-07

Patch release [pgscatalog-match-0.3.2](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.match-0.3.2)

* Bumps `pgscatalog.core` dependency
  
Minor release [pgscatalog-calc-0.3.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.calc-0.3.0)

* Add `--verify_variants` CLI option to `pgscatalog-aggregate` to confirm that scoring file variants (input) match the scored variants perfectly during score calculation 

Patch release  [pgscatalog-core-0.3.1](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.core-0.3.1)

* Convert all relative imports to absolute to fix broken autoapi documentation

Minor release  [pgscatalog-core-0.3.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.core-0.3.0)

* Implements unified PGS Catalog data models with pydantic (used in `pgscatalog-combine`)
* Rate limits `pgscatalog-download` to prevent throttling

## 2024-08-05

Patch release  [pgscatalog-utils-1.3.1](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog-utils-1.3.1)

* Bump dependencies to latest

Patch release  [pgscatalog-match-0.3.1](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.match-0.3.1)

* Fix splitting but revert implementation to original function. v0.3.0 would trigger ValueErrors consistently on larger datasets
* Fix variant match information getting dropped from the log if an accession was filtered from output scoring files 

## 2024-08-02

Minor release  [pgscatalog-utils-1.3.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog-utils-1.3.0)

* Bump dependencies to latest

Minor release  [pgscatalog-match-0.3.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.match-0.3.0)

* **Important** Fix splitting duplicate variant IDs across different scoring files https://github.com/PGScatalog/pygscatalog/pull/38
* Fix `--keep_ambiguous` https://github.com/PGScatalog/pgsc_calc/issues/346
* Fix `--ignore_strand_flips`
  
## 2024-07-29

Minor release  [pgscatalog-utils-1.2.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog-utils-1.2.0)

* Bump dependencies to latest

Patch release [pgscatalog-match-0.2.3](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.match-0.2.3)

* Update `pgscatalog-intersect` PCA_ELIGIBLE filters https://github.com/PGScatalog/pygscatalog/pull/33
   
Patch release [pgscatalog-calc-0.2.2](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.calc-0.2.2)

* Sample IDs now use FIDs if available https://github.com/PGScatalog/pygscatalog/pull/31

Patch release [pgscatalog-core-0.2.2](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.core-0.2.2)

* Fixes support for missing other_allele column in custom scoring files https://github.com/PGScatalog/pygscatalog/pull/30 
   
## 2024-06-19

Patch release [pgscatalog-match-0.2.2](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.match-0.2.1)

* Fixes `pgscatalog-intersect` with bim target input https://github.com/PGScatalog/pygscatalog/pull/25
 
## 2024-06-14

Patch releases to downgrade minimum Python version to 3.10:

*  [pgscatalog-utils-1.1.1](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog-utils-1.1.1)
*  [pgscatalog-match-0.2.1](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.match-0.2.1)
*  [pgscatalog-calc-0.2.1](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.calc-0.2.1)
*  [pgscatalog-core-0.2.1](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog.core-0.2.1)

## 2024-06-13

Minor release [pgscatalog-utils-1.1.0](https://github.com/PGScatalog/pygscatalog/releases/tag/pgscatalog-utils-1.1.0)

* Bump dependencies to latest

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
