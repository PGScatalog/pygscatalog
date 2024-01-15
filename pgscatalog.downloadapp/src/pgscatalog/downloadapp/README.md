# `pgscatalog-download`

This simple CLI application makes it easy to download scoring files in bulk from the [PGS Catalog](https://www.pgscatalog.org/) with a mixture of PGS, publication, or trait accessions.

```
$ pip install pgscatalog-download
$ pgscatalog-download -i PGS000822 -b GRCh38 -o scorefiles/
```

## Help

```
Download a set of scoring files from the PGS Catalog using PGS
Scoring IDs, traits, or publication IDs.

The PGS Catalog API is queried to get a list of scoring file
URLs. Scoring files are downloaded via FTP to a specified
directory. PGS Catalog scoring files are staged with the name:

        {PGS_ID}.txt.gz

If a valid build is specified harmonized files are downloaded as:

    {PGS_ID}_hmPOS_{genome_build}.txt.gz

These harmonised scoring files contain genomic coordinates,
remapped from author-submitted information such as rsids.

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
                        Download Harmonized Scores with Positions in Genome build: GRCh37 or GRCh38
  -o OUTDIR, --outdir OUTDIR
                        <Required> Output directory to store downloaded files
  -w, --overwrite       <Optional> Overwrite existing Scoring File if a new version is available for download on the FTP
  -c PGSC_CALC, --pgsc_calc PGSC_CALC
                        <Optional> Provide information about downloading scoring files via pgsc_calc
```