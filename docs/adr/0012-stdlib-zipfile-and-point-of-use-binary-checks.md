# Archive the cache with stdlib `zipfile`, and check external binaries where they are used

## Status

Accepted

## Context

LOAD writes the cache as a directory store and archives it to `genotypes.zarr.zip`, which is what a Nextflow task stages to MATCH and SCORE ([ADR 0003](0003-shards-are-zarr-leaf-groups.md)). Both directions shelled out to `7z`, and `cli/load.py` refused to import at all unless `7z` and `bgenix` were both on `PATH`.

The import-time check is a defect independent of any archiving decision: `cli/main.py` imports `load_cli` to build its subparsers, so a machine without `bgenix` could not run `pgsc_calc combine`, or even `--help`, despite bgen input being nowhere in sight.

The archiving choice deserved re-examination because of what `7z` was actually asked to do. The archive is created with `-mx0` — no compression, because zarr has already compressed the arrays — so nothing is being traded for CPU time. Extraction of a single leaf group was expressed as `-ir!` include patterns plus a loop appending each ancestor's `zarr.json`, which means the set of members to extract was already being computed in Python and handed to `7z` as argv. The stated reason for the dependency was that "zipfile in standard library and other python libraries are pure Python. 7z is written in C++ so shouldn't have problems scaling" — a hypothesis about throughput that does not survive `-mx0`, where the bytes are moved by `shutil.copyfileobj` and checksummed by zlib, both C. What remains genuinely faster in `7z` is per-entry bookkeeping, which matters in proportion to the *number of chunks* in a store rather than its size.

Against that, `7z` is the one thing standing between a bare Python environment and a working install. That matters now because the four stages are also being made runnable in sequence in a single process, for environments where Nextflow cannot run (see the Serial run entry in `CONTEXT.md`), and the audience for that entry point is exactly the audience that cannot install system packages.

## Decision

`zip_zarr` and `unzip_zarr` use `zipfile` from the standard library, with `ZIP_STORED`. The `7z` dependency and both import-time `shutil.which` checks are removed.

`bgenix` is checked at the top of `run_bgenix_subprocess`, the function that invokes it. No check is performed for input formats that never reach it.

**A runtime fallback between `7z` and `zipfile` was rejected.** Two archiving paths means the one nobody runs locally is the one that breaks, for a benefit that has not been measured.

**Zarr's native `ZipStore` was rejected.** LOAD needs a directory store precisely so that an existing cache can be updated via `--zarr_zip_file`, which `ZipStore` does not support.

## Consequences

- `pip install pgscatalog.calc` yields a working VCF pipeline with no system packages. `pgsc_calc` still ships containers, and its own Nextflow modules still shell out to `7z` directly (`score/main.nf` zips both the cache and `variant_match_logs.zip`), so this removes a dependency from the Python distribution, not from the pipeline's images.
- Selective extraction becomes a filter over `namelist()`, so the ancestor-`zarr.json` rule is expressed as a prefix test rather than as accumulated argv.
- **The crossover point is unmeasured.** Stores with very large chunk counts are where stdlib per-entry overhead could show up, and the signal would be LOAD spending time after the arrays are written. If it ever appears, the answer is a larger chunk size rather than a second archiver.
- A missing `bgenix` is now discovered when a bgen file is read rather than at startup. This is later, and it is the first point at which the requirement is real.
