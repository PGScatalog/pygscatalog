import enum


class PopulationType(enum.Enum):
    TARGET = "target"
    REFERENCE = "reference"


class PrincipalComponents:
    def __init__(
        self, pcs_path, psam_path, npcs_popcomp, npmcs_normalisation, pop_type
    ):
        self.pcs_path = pcs_path
        self.psam_path = psam_path
        self.max_pcs = max([10, npcs_popcomp, npmcs_normalisation])
        self.pop_type = pop_type
