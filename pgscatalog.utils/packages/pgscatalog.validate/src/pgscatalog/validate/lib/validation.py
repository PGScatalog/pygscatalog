from os import PathLike
from pathlib import Path
from typing import Dict

from pgscatalog.core.lib.models import VariantType, CatalogScoreHeader
from pgscatalog.validate.lib.models import *
from pgscatalog.validate.lib.parser import ScoreFileParser, SpreadsheetFileParser, TSVFileParser


class ScoringFileValidation:
    """Main class for performing the validation of PGS Scoring files.

    The validation is performed when the object is built.

    >>> validation = ScoringFileValidation("path_to_scoring_file", header=False)
    >>> validation.is_valid()
    True

    The validation errors are then available by reading the .errors attribute.

    >>> for validation_error in validation.errors:
    >>>     print(f"Row: {validation_error.row}, Error: {str(validation_error.error)}")
    """

    column_names: ColumnNames
    accession: str
    errors: list[ScoringFileValidationError]
    warnings: list[str]
    parser: ScoreFileParser

    def __init__(self, file_name: str | PathLike[str], header: bool = False, hm: bool = False, strict: bool = False):

        path = Path(file_name)
        if path.suffix == '.xlsx':
            self.parser = SpreadsheetFileParser(path)
        else:
            self.parser = TSVFileParser(path)

        self.accession = self.parser.file_path.stem
        self.errors = []
        self.warnings = []

        if header:
            try:
                CatalogScoreHeader.from_path(path)
            except ValidationError as e:
                self.errors.append(ScoringFileValidationError(0, e))

        # Setting static variables (meant to be changed in the future to something more elegant)
        ValidationModel.strict = strict
        ValidationModel.warnings = []

        try:
            column_names = self.parser.get_column_names()
            # Validate the column names, and interrupt if not valid
            try:
                self.column_names = ColumnNames(columns=column_names, hm=hm)
            except ValidationError as e:
                self.errors.append(ScoringFileValidationError(0, e))
                return

            # Validate the variants
            for i, line in enumerate(self.parser.get_variants(), 1):
                line = self.__remove_empty_values(line)
                try:
                    if 'variant_type' in line and line['variant_type'] in VariantType:
                        ComplexVariant(**line, **{"accession": self.accession, "row_nr": i})
                    else:
                        ValidationVariant(**line, **{"accession": self.accession, "row_nr": i})
                except ValidationError as e:
                    self.errors.append(ScoringFileValidationError(i, e))

            self.warnings = ValidationModel.warnings

        except FileNotFoundError:
            print(f"Error: File not found: {self.parser.file_path}")

    def is_valid(self) -> bool:
        """Checks if the validated file is valid. If not, the errors available through the .errors attribute."""
        return not self.errors

    @staticmethod
    def __remove_empty_values(line: dict) -> Dict:
        return {key: value for key, value in line.items() if value is not None and value != ""}