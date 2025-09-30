from __future__ import annotations

import os
import pathlib

type Pathish = str | pathlib.Path | os.PathLike[str]
type PathishList = list[Pathish]
