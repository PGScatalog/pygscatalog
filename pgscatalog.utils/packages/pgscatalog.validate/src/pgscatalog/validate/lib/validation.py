from os import PathLike
from pathlib import Path
from typing import Dict, Iterator

from pgscatalog.core.lib.models import VariantType, CatalogScoreHeader
from pgscatalog.validate.lib.models import *
from pgscatalog.validate.lib.parser import ScoreFileParser, SpreadsheetFileParser, TSVFileParser


class ScoringFileValidation:
    """Main class for performing the validation of PGS Scoring files.

    The validation is performed when the object is built.

    >>> validation = ScoringFileValidation("path_to_scoring_file", header=False)

    The validation errors are then available by reading the .errors attribute.

    >>> for validation_error in validation.get_errors():
    >>>     print(f"Row: {validation_error.row}:")
    >>>     for error in validation_error.errors:
    >>>         print(f"\t{error.attr}: {error.msg}")
    """

    column_names: ColumnNames
    accession: str
    parser: ScoreFileParser

    def __init__(self, file_name: str | PathLike[str], header: bool = False,
                 hm: bool = False, strict: bool = False) -> None:

        path = Path(file_name)
        self.path = path
        if path.suffix == '.xlsx':
            self.parser = SpreadsheetFileParser(path)
        else:
            self.parser = TSVFileParser(path)

        self.accession = self.parser.file_path.stem
        self.header = header
        self.hm = hm
        self.strict = strict

    def get_errors(self) -> Iterator[ScoringFileValidationError]:
        """Generator that yields the validation errors found in the file.
        The validation is performed when the generator is called, and the errors are yielded one by one.
        This allows to stop the validation after a certain number of errors."""

        if self.header:
            try:
                CatalogScoreHeader.from_path(self.path)
            except ValidationError as e:
                yield ScoringFileValidationError(0, e)

        # Setting "strict" static variable, so it can be read by the models during validation.
        ValidationModel.strict = self.strict

        try:
            column_names = self.parser.get_column_names()
            # Validate the column names, and interrupt if not valid
            try:
                self.column_names = ColumnNames(columns=column_names, hm=self.hm)
            except ValidationError as e:
                yield ScoringFileValidationError(0, e)

            # Validate the variants
            for i, line in enumerate(self.parser.get_variants(), 1):
                line = self.__remove_empty_values(line)

                if 'variant_type' in line and line['variant_type'] in VariantType:
                    yield from self.__validate_model(ComplexVariant, {**line, **{"accession": self.accession, "row_nr": i}})
                else:
                    yield from self.__validate_model(ValidationVariant, {**line, **{"accession": self.accession, "row_nr": i}})

        except FileNotFoundError:
            print(f"Error: File not found: {self.parser.file_path}")

    @staticmethod
    def __remove_empty_values(line: dict) -> Dict:
        return {key: value for key, value in line.items() if value is not None and value != ""}

    def __validate_model(self, model_type: type[ValidationModel], data: dict) -> Iterator[ScoringFileValidationError]:

        # The context is used to store the warnings that are generated during the validation,
        # as pydantic does not have a built-in way to handle warnings.
        # The context is passed to the model_validate method, and the warnings are added to it during the validation process.
        # After the validation, the warnings are checked and yielded as ScoringFileValidationError with level WARNING if there are any.
        context = {}

        try:
            model_type.model_validate(data, context=context)
        except ValidationError as e:
            yield ScoringFileValidationError(data.get("row_nr", 0), e)

        if "warnings" in context and context["warnings"]:
            yield ScoringFileValidationError(data.get("row_nr", 0), context["warnings"], ErrorLevel.WARNING)
