from __future__ import annotations

import enum


class GenomeFileType(str, enum.Enum):
    VCF = "vcf"
    BGEN = "bgen"
