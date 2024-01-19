# `pgscatalog-combine`

[![.github/workflows/combineapp-pytest.yml](https://github.com/PGScatalog/pygscatalog/actions/workflows/combineapp-pytest.yml/badge.svg)](https://github.com/PGScatalog/pygscatalog/actions/workflows/combineapp-pytest.yml)

A CLI application that combines polygenic score (PGS) scoring files into a standardised output. 

The process involves:

* Normalising scoring files to a consistent data format
* Doing some quality control checks on variants
* Writing a long format / melted file 

Input scoring files must follow [PGS Catalog standards](https://www.pgscatalog.org/downloads/#dl_scoring_files). The output file is useful for doing data science tasks, like matching variants across a scoring file and target genome. 

The [PGS Catalog Calculator](https://github.com/PGScatalog/pgsc_calc) uses this CLI application internally.

## Installation

```
$ pip install pgscatalog-combine
```

## Example output

Outputs a tab delimited (and optionally compressed) text file. 

```
chr_name  chr_position  effect_allele  other_allele  effect_weight         effect_type  is_duplicated  accession  row_nr
11        69516650      T              C             0.16220387987485377   additive     False          PGS000001  0
11        69564393      A              C             0.023618865598634027  additive     False          PGS000001  1
...
11        69516650      T              C             0.2104229869048294      additive     False          PGS000002  0
11        69564393      A              C             0.01872361399810251     additive     False          PGS000002  1
```

And a JSON log:

```
[
    {
        "PGS000001": {
            "pgs_id": "PGS000001",
            "pgp_id": "PGP000001",
            "pgs_name": "PRS77_BC",
            "genome_build": "GRCh38",
            "variants_number": "77",
...
```

## Help

```
$ pgscatalog-combine --help

Combine multiple scoring files in PGS Catalog format (see
https://www.pgscatalog.org/downloads/ for details) to a 'long' table of columns
needed for variant matching and subsequent calculation.

Custom scorefiles in PGS Catalog format can be combined with PGS Catalog scoring files, and
optionally liftover genomic coordinates to GRCh37 or GRCh38. The script can accept a mix of
unharmonised and harmonised PGS Catalog data. By default all variants are output (including
positions with duplicated data [often caused by rsID/liftover collions across builds]) and
variants with missing positions.

options:
  -h, --help            show this help message and exit
  -s SCOREFILES [SCOREFILES ...], --scorefiles SCOREFILES [SCOREFILES ...]
                        <Required> Scorefile paths
  --liftover            <Optional> Convert scoring file variants to target genome build?
  -t {GRCh37,GRCh38}, --target_build {GRCh37,GRCh38}
                        <Required> Build of target genome
  -c CHAIN_DIR, --chain_dir CHAIN_DIR
                        Path to directory containing chain files
  -m MIN_LIFT, --min_lift MIN_LIFT
                        <Optional> If liftover, minimum proportion of variants lifted over
  --drop_missing        <Optional> Drop variants with missing information (chr/pos) and non-standard alleles (e.g. HLA=P/N) from the output file.
  -o OUTFILE, --outfile OUTFILE
                        <Required> Output path to combined long scorefile [ will compress output if filename ends with .gz ]
  -l LOGFILE, --logfile LOGFILE
                        <Required> Name for the log file (score metadata) for combined scores.[ will write to identical directory as combined scorefile]
  -v, --verbose         <Optional> Extra logging information

The long table is used to simplify intersecting variants in target genotyping datasets
and the scoring files with the match_variants program.
```
