import polars as pl


class MatchResult:
    """This class represents variant matching results"""

    def __init__(self, scorefile, dataset, matchresult=None, df=None):
        if matchresult is None and df is None:
            raise ValueError("Please provide either matchresult or df")

        self.df = df
        self.scorefile = scorefile
        self.dataset = dataset
        self._lazydf = matchresult

    def collect(self, outfile=None):
        """Compute match results and optionally save to file"""
        if self._lazydf is None:
            raise ValueError("Can't collect, missing matchresult")

        if outfile is not None:
            pl.concat(pl.collect_all(self._lazydf)).write_ipc(outfile)
            self.df = pl.scan_ipc(outfile)
        else:
            # call .lazy() to prevent eager queries downstream
            self.df = pl.concat(pl.collect_all(self._lazydf)).lazy()

        return self.df

    @classmethod
    def from_ipc(cls, matchresults_ipc_path, scorefile_ipc_path, dataset):
        return cls(
            df=pl.scan_ipc(matchresults_ipc_path),
            scorefile=pl.scan_ipc(scorefile_ipc_path),
            dataset=dataset,
        )

    def __add__(self, other):
        if not isinstance(other, MatchResult):
            return NotImplemented

        if self.df or other.df is None:
            raise ValueError("Can't add a missing df")

        # TODO: add different scoring files? append to a list?
        # MatchResult(df=pl.concat([self.df, other.df]))
