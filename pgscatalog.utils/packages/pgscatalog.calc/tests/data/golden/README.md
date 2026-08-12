# Golden scores

`scores.txt.gz` is the output of the fixture data scored by a single process, before
SCORE was split into MATCH / MATCH_COMBINE / SCORE / SCORE_COMBINE (docs/adr/0003).
Once `score` only ever sees one shard there is no single-process reference left, so
this file is it.

Every target genome shape in `tests/integration/conftest.py` holds the same
(sample, variant) pairs, so all of them must reproduce these numbers.

To regenerate, run `pgsc_calc run` on the `unsplit_target` fixture (see
`run_serial` in `tests/integration/conftest.py`) and copy the resulting
`scores.txt.gz` here. Only do this if the change to the numbers is intended and
understood.
