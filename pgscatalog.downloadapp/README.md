# `pgscatalog-download`

[![.github/workflows/downloadapp-pytest.yml](https://github.com/PGScatalog/pygscatalog/actions/workflows/downloadapp-pytest.yml/badge.svg)](https://github.com/PGScatalog/pygscatalog/actions/workflows/downloadapp-pytest.yml)

This simple CLI application makes it easy to download scoring files in bulk from the [PGS Catalog](https://www.pgscatalog.org/) with a mixture of PGS, publication, or trait accessions.

```
$ pip install pgscatalog-download
$ mkdir scorefiles
$ pgscatalog-download --pgs PGS000822 --pgp PGP000001 PGP000002 --efo EFO_0004214 -b GRCh38 -o scorefiles/
$ tree scorefiles
scorefiles
├── PGS000001_hmPOS_GRCh38.txt.gz
├── PGS000002_hmPOS_GRCh38.txt.gz
├── PGS000003_hmPOS_GRCh38.txt.gz
...
```

## Help

```
$ pgscatalog-download --help

usage: pgscatalog-download [-h] [-i PGS [PGS ...]] [-t EFO [EFO ...]] [-e] [-p PGP [PGP ...]] [-b {GRCh37,GRCh38}] -o OUTDIR [-w] [-c USER_AGENT] [-v]

Download a set of scoring files from the PGS Catalog using PGS Scoring
IDs, traits, or publication accessions.

The PGS Catalog API is queried to get a list of scoring file URLs.
Scoring files are downloaded asynchronously via HTTPS to a specified
directory. Downloaded files are automatically validated against an md5
checksum.

PGS Catalog scoring files are staged with the name:

    {PGS_ID}.txt.gz

If a valid build is specified harmonized files are downloaded as:

    {PGS_ID}_hmPOS_{genome_build}.txt.gz

These harmonised scoring files contain genomic coordinates, remapped
from author-submitted information such as rsIDs.

options:
  -h, --help            show this help message and exit
  -i PGS [PGS ...], --pgs PGS [PGS ...]
                        PGS Catalog ID(s) (e.g. PGS000001)
  -t EFO [EFO ...], --efo EFO [EFO ...]
                        Traits described by an EFO term(s) (e.g. EFO_0004611)
  -e, --efo_direct      <Optional> Return only PGS tagged with exact EFO term (e.g. no PGS for child/descendant terms in the ontology)
  -p PGP [PGP ...], --pgp PGP [PGP ...]
                        PGP publication ID(s) (e.g. PGP000007)
  -b {GRCh37,GRCh38}, --build {GRCh37,GRCh38}
                        Download harmonized scores with positions in genome build: GRCh37 or GRCh38
  -o OUTDIR, --outdir OUTDIR
                        <Required> Output directory to store downloaded files
  -w, --overwrite       <Optional> Overwrite existing Scoring File if a new version is available for download on the FTP
  -c USER_AGENT, --user_agent USER_AGENT
                        <Optional> Provide custom user agent when querying PGS Catalog API
  -v, --verbose         <Optional> Extra logging information
```