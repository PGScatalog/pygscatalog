import csv

from xopen import xopen

from pgscatalog.core.lib.models import ScoreVariant


def _read_normalised_rows(path):
    with xopen(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield ScoreVariant(**row)


class NormalisedScoringFile:
    """This class represents a ScoringFile that's been normalised to have a consistent format

    Its main purpose is to provide a convenient way to iterate over variants

    # TODO: replace with a pydantic model in pgscatalog.core
    """

    def __init__(self, path):
        try:
            with xopen(path):
                pass
        except TypeError:
            self.is_path = False
            self.path = str(path)
        else:
            self.is_path = True
            self.path = path
        finally:
            # either a ScoringFile or a path to a combined file
            self._scoringfile = path

    def __iter__(self):
        yield from self.variants

    @property
    def variants(self):
        if self.is_path:
            # get a fresh generator from the file
            self._variants = _read_normalised_rows(self._scoringfile)
        else:
            # get a fresh generator from the normalise() method
            self._variants = self._scoringfile.normalise()

        return self._variants

    def __repr__(self):
        if self.is_path:
            x = f"{repr(str(self._scoringfile))}"
        else:
            x = f"{repr(self._scoringfile)}"

        return f"{type(self).__name__}({x})"
