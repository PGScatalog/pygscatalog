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
    parser: ScoreFileParser

    def __init__(self, file_name: str, header: bool = False, hm: bool = False):

        if file_name.endswith('.xlsx'):
            self.parser = SpreadsheetFileParser(file_name)
        else:
            self.parser = TSVFileParser(file_name)

        self.accession = self.parser.file_path.stem
        self.errors = []

        if header:
            try:
                CatalogScoreHeader.from_path(file_name)
            except ValidationError as e:
                self.errors.append(ScoringFileValidationError(0, e))

        try:
            column_names = self.parser.get_column_names()
            # Validate the column names, and interrupt if not valid
            try:
                self.column_names = ColumnNames(columns=column_names, hm=hm)
            except ValidationError as e:
                self.errors.append(ScoringFileValidationError(0, e))
                return

            # Validate the variants
            index = 0
            for line in self.parser.get_variants():
                line = self.__remove_empty_values(line)
                index += 1
                try:
                    if 'variant_type' in line and line['variant_type'] in VariantType:
                        ComplexVariant(**line, **{"accession": self.accession, "row_nr": index})
                    else:
                        ValidationVariant(**line, **{"accession": self.accession, "row_nr": index})
                except ValidationError as e:
                    self.errors.append(ScoringFileValidationError(index, e))

        except FileNotFoundError:
            print(f"Error: File not found: {self.parser.file_path}")

    def is_valid(self) -> bool:
        """Checks if the validated file is valid. If not, the errors available through the .errors attribute."""
        return not self.errors

    @staticmethod
    def __remove_empty_values(line: dict) -> Dict:
        return {key: value for key, value in line.items() if value is not None and value != ""}