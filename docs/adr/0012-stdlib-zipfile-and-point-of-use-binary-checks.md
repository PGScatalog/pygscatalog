# Archive the cache with stdlib `zipfile` and check external binaries at use

## Status

Accepted

## Context

LOAD writes the cache as a directory store and archives it to
`genotypes.zarr.zip`. Nextflow stages that archive to MATCH and SCORE
([ADR 0003](0003-shards-are-zarr-leaf-groups.md)). Both archive operations invoked
`7z`. `cli/load.py` refused to import unless both `7z` and `bgenix` were on `PATH`.

The import-time check is a separate defect. `cli/main.py` imports `load_cli` to
build its subparsers. Without `bgenix`, a machine could not run `pgsc_calc combine`
or even `--help`. Neither command needs bgen input.

The archive uses `-mx0`, which disables compression. Zarr already compressed the
arrays, so compression creates no CPU trade-off here. Selective extraction used
`-ir!` include patterns. A loop added each ancestor's `zarr.json`. Python therefore
already calculated the members to extract and passed them to `7z` as argv.

The original dependency rationale compared pure-Python ZIP libraries with `7z`,
which uses C++. It assumed that C++ would avoid scaling problems. That throughput assumption
does not apply with `-mx0`. `shutil.copyfileobj` moves the bytes and zlib calculates
the checksums, both in C. `7z` remains faster at per-entry bookkeeping. That
advantage depends on the number of chunks rather than the store's size.

The `7z` dependency prevents installation in a bare Python environment. This
matters because the four stages are also becoming available sequentially in one
process. That entry point serves environments that cannot run Nextflow or install
system packages. See the Serial run entry in `CONTEXT.md`.

## Decision

`zip_zarr` and `unzip_zarr` use the standard library's `zipfile` with `ZIP_STORED`.
The change removes the `7z` dependency and both import-time `shutil.which` checks.

`run_bgenix_subprocess` checks for `bgenix` at its start, before invoking it. Input
formats that never reach this function need no check.

**A runtime fallback between `7z` and `zipfile` was rejected.** The archive path
that nobody uses locally is the one that breaks. Its benefit is unmeasured.

**Zarr's native `ZipStore` was rejected.** LOAD needs a directory store to update
an existing cache through `--zarr_zip_file`. `ZipStore` does not support that
operation.

## Consequences

- `pip install pgscatalog.calc` provides a working VCF pipeline without system
  packages. `pgsc_calc` still ships containers. Its Nextflow modules still invoke
  `7z` directly. In `score/main.nf`, `7z` archives both the cache and
  `variant_match_logs.zip`. The change removes a Python distribution dependency,
  not an image dependency.
- Selective extraction filters `namelist()`. A prefix test expresses the
  ancestor-`zarr.json` rule, instead of accumulating argv.
- **The crossover point is unmeasured.** Standard-library overhead per entry could
  matter for stores with very large chunk counts. The symptom would be LOAD
  spending time after it writes the arrays. If that occurs, increase the chunk
  size instead of adding a second archiver.
- A bgen read now detects a missing `bgenix`, rather than detecting it at startup.
  Detection happens later, at the first operation that actually requires it.
