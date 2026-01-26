import argparse
import logging
import sys
import textwrap
import warnings
from pathlib import Path
from typing import TextIO

from pgscatalog.validate.lib.models import ScoringFileValidationError
from pgscatalog.validate.lib.validation import ScoringFileValidation

data_sum = {'valid': [], 'invalid': [], 'other': []}

val_types = ('formatted', 'hm_pos', 'raw')

logging.basicConfig(level=logging.INFO, format='(%(levelname)s): %(message)s')


def run() -> None:
    args = _parse_args()
    _check_args(args)

    validator_type: str = args.type
    filename: Path = args.filename
    dirname: Path = args.dir
    log_dir: Path = args.log_dir
    score_dir: Path = args.score_dir
    check_filename: bool = args.check_filename
    strict: bool = args.strict

    # Check PGS Catalog file name nomenclature
    if not check_filename and validator_type != 'raw':
        # Raw files can have any name
        print("WARNING: the parameter '--check_filename' is not present in the submitted command line, " +
              "therefore the validation of the scoring file name(s) won't be performed.")

    ## Select validator class ##
    match validator_type:
        case 'raw':
            header = False
        case _:
            header = True

    ## Run validator ##
    # One file
    if filename:
        _run_validator(filename, log_dir, score_dir, check_filename, header, strict)
    # Content of the directory
    elif dirname:
        count_files = 0
        valid_files = []
        invalid_files = []
        other_files = []
        # Browse directory: for each file run validator
        for filepath in sorted(dirname.glob("*.*")):
            count_files += 1
            try:
                summary = _run_validator(filepath, log_dir, score_dir, check_filename, header, strict)
            except Exception as e:
                warnings.warn(str(e))
                other_files.append(filepath)
                continue
            if summary['errors_count'] > 0:
                invalid_files.append(filepath)
            else:
                valid_files.append(filepath)

        # Print summary  + results
        print("\nSummary:")
        if valid_files:
            print(f"- Valid: {len(valid_files)}/{count_files}")
        if invalid_files:
            print(f"- Invalid: {len(invalid_files)}/{count_files}")
        if other_files:
            print(f"- Other issues: {len(other_files)}/{count_files}")

        if invalid_files:
            print("\nInvalid files:")
            print("\n".join(map(str, invalid_files)))


def _print_errors(validation_errors, file: TextIO = None):
    has_errors = False
    for validation_error in validation_errors:
        has_errors = True
        row = validation_error.row
        errors = validation_error.error.errors()
        for error in errors:
            print('Line {}{}: {}'.format(row, ' ' + error['loc'][0] if 'loc' in error and error['loc'] else '',
                                         error['msg']), file=file)
    return has_errors


def _report_errors_to_log_file(validation_errors: list[ScoringFileValidationError], log_file: str):
    with open(log_file, "w") as f:
        has_errors = _print_errors(validation_errors, f)
        if has_errors:
            print('File is invalid', file=f)
        else:
            print('File is valid', file=f)


def _report_errors_to_stdout(validation_errors: list[ScoringFileValidationError]):
    _print_errors(validation_errors, None)


def _report_warnings(warning_messages: list[str], file=sys.stdout) -> None:
    """Print the given warnings to stdout or the given file if provided."""
    print('- Warnings:', file=file)
    for warning in warning_messages:
        print(warning, file=file)


def _check_args(args: argparse.Namespace) -> None:
    ## Check parameters ##
    # Scoring files directory (only to compare with the harmonized files)
    if args.type == 'hm_pos' and not args.score_dir:
        print("WARNING: the parameter '--score_dir' is not present in the submitted command line,"
              " therefore the comparison of the number of data rows between the formatted scoring file(s)"
              " and the harmonized scoring file(s) won't be performed.")


def _run_validator(filepath: Path, log_dir: Path, score_dir: Path, check_filename: bool,
                   header: bool, strict: bool = False) -> dict:
    """Run the file validator"""
    # TODO: add check_filename and score_dir support
    file = filepath.name
    filename = filepath.stem
    print(f"# Filename: {file}")
    if log_dir:
        log_file = f'{log_dir}/{filename}_log.txt'
    else:
        log_file = None

    validation = ScoringFileValidation(filepath, header=header, strict=strict)

    summary = {'errors_count': len(validation.errors)}

    # Check log
    if log_file:
        _report_errors_to_log_file(validation.errors, log_file)
        summary['log_file'] = log_file
        if validation.warnings:
            with open(log_file, 'a') as f:
                _report_warnings(validation.warnings, file=f)
    else:
        # Report to stdout
        _report_errors_to_stdout(validation.errors)
        if validation.warnings:
            _report_warnings(validation.warnings)

    return summary


def _description_text() -> str:
    return textwrap.dedent('''\
    Validate a set of scoring files to match the PGS Catalog scoring file formats.
    It can validate:
    - The formatted scoring file format (https://www.pgscatalog.org/downloads/#dl_ftp_scoring)
    - The harmonized (Position) scoring file format (https://www.pgscatalog.org/downloads/#dl_ftp_scoring_hm_pos)
   ''')


def _epilog_text() -> str:
    return textwrap.dedent(f'''\
    You need to specify the type of file format to validate, using the parameter '-t' ({' or '.join(val_types)}).
   ''')


def _valid_file(file: str) -> Path:
    path = Path(file)
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"{file} is not a valid file")
    return path


def _valid_dir(directory: str) -> Path:
    path = Path(directory)
    if not path.is_dir():
        raise argparse.ArgumentTypeError(f"{directory} is not a directory")
    return path


def _parse_args(args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=_description_text(), epilog=_epilog_text(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-t", "--type",
                        help=f"Type of validator: {' or '.join(val_types)}", metavar='VALIDATOR_TYPE',
                        choices=val_types, default='raw')
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-f", "--filename",
                             type=_valid_file,
                             help="""The path to the polygenic scoring file to be validated 
                             (no need to use the [--dir] option)""",
                             metavar='SCORING_FILE_NAME')
    input_group.add_argument("-d", "--dir",
                             type=_valid_dir,
                             help="""The name of the directory containing the files that need to processed 
                             (no need to use the [-f] option""")
    parser.add_argument("-D", "--score_dir",
                        type=_valid_dir,
                        help="""<Optional> The name of the directory containing the formatted scoring files
                         to compare with harmonized scoring files (only if "--type hm_pos" is used)""")
    parser.add_argument("-l", "--log_dir",
                        type=_valid_dir,
                        help='<Optional> The name of the log directory where the log file(s) will be stored',
                        required=False)
    parser.add_argument("-c", "--check_filename",
                        help="""<Optional> Check that the file name match the PGS Catalog nomenclature. 
                        Only if using with --type [formatted|hm_pos]""",
                        required=False,
                        action='store_true')
    parser.add_argument("-s", "--strict",
                        help="""<Optional> Treat warnings (such as leading/trailing spaces) as errors""",
                        required=False,
                        action='store_true')
    return parser.parse_args(args)


if __name__ == '__main__':
    run()
