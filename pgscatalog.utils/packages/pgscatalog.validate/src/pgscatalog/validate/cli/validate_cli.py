import argparse
import glob
import logging
import os
import re
import textwrap
from typing import TextIO

from pgscatalog.validate.lib.models import ScoringFileValidationError
from pgscatalog.validate.lib.validation import ScoringFileValidation

data_sum = {'valid': [], 'invalid': [], 'other': []}

val_types = ('formatted', 'hm_pos', 'raw')

logging.basicConfig(level=logging.INFO, format='(%(levelname)s): %(message)s')


def run() -> None:
    global data_sum, score_dir
    args = _parse_args()
    _check_args(args)

    validator_type = args.type
    files_dir = args.dir
    log_dir = args.log_dir

    # Check PGS Catalog file name nomenclature
    check_filename = False
    if args.check_filename:
        check_filename = True
    elif validator_type != 'raw':
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
    if args.filename:
        _run_validator(args.filename, log_dir, score_dir, check_filename, header)
    # Content of the directory
    elif files_dir:
        count_files = 0
        # Browse directory: for each file run validator
        for filepath in sorted(glob.glob(files_dir+"/*.*")):
            _run_validator(filepath, log_dir, score_dir, check_filename, header)
            count_files += 1

        # Print summary  + results
        print("\nSummary:")
        if data_sum['valid']:
            print(f"- Valid: {len(data_sum['valid'])}/{count_files}")
        if data_sum['invalid']:
            print(f"- Invalid: {len(data_sum['invalid'])}/{count_files}")
        if data_sum['other']:
            print(f"- Other issues: {len(data_sum['other'])}/{count_files}")

        if data_sum['invalid']:
            print("Invalid files:")
            print("\n".join(data_sum['invalid']))


def _read_last_line(file: str) -> str:
    """
    Return the last line of the file
    """
    file_handle = open(file, "r")
    lines = file_handle.readlines()
    file_handle.close()
    return lines[-1]


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


def _file_validation_state(filename: str, log_file: str) -> None:
    global data_sum
    if os.path.exists(log_file):
        log_result = _read_last_line(log_file)
        if re.search("File is valid", log_result):
            print("> valid\n")
            data_sum['valid'].append(filename)
        elif re.search("File is invalid", log_result):
            print("#### invalid! ####\n")
            data_sum['invalid'].append(filename)
        else:
            print("!! validation process had an issue. Please look at the logs.\n")
            data_sum['other'].append(filename)
    else:
        print("!! validation process had an issue: the log file can't be found")
        data_sum['other'].append(filename)


def _check_args(args: argparse.Namespace) -> None:
    global score_dir
    ## Check parameters ##
    # Type of validator
    if args.type and args.type not in val_types:
        print(f"Error: Validator type (option -t) '{args.t}' is not in the list of recognized types: {val_types}.")
        exit(1)
    # Logs dir
    if args.log_dir and not os.path.isdir(args.log_dir):
        print(f"Error: Log dir '{args.log_dir}' can't be found!")
        exit(1)
    # File and directory parameters (only one of the '-f' and '--dir' can be used)
    if args.filename and args.dir:
        print("Error: you can't use both options [-f] - single scoring file and [--dir] - directory of scoring files."
              " Please use only 1 of these 2 options!")
        exit(1)
    elif not args.filename and not args.dir:
        print("Error: you need to provide a scoring file [-f] or a directory of scoring files [--dir]!")
        exit(1)
    elif args.filename and not os.path.isfile(args.filename):
        print(f"Error: Scoring file '{args.filename}' can't be found!")
        exit(1)
    elif args.dir and not os.path.isdir(args.dir):
        print(f"Error: the scoring file directory '{args.dir}' can't be found!")
        exit(1)
    # Scoring files directory (only to compare with the harmonized files)
    score_dir = None
    if args.score_dir:
        score_dir = args.score_dir
        if not os.path.isdir(score_dir):
            print(f"Error: Scoring file directory '{score_dir}' can't be found!")
            exit(1)
    elif args.type not in ['formatted', 'raw']:
        print("WARNING: the parameter '--score_dir' is not present in the submitted command line,"
              " therefore the comparison of the number of data rows between the formatted scoring file(s)"
              " and the harmonized scoring file(s) won't be performed.")


def _run_validator(filepath: str, log_dir: str, score_dir: str, check_filename: bool, header: bool) -> None:
    """Run the file validator"""
    # TODO: add check_filename and score_dir support
    file = os.path.basename(filepath)
    filename = file.split('.')[0]
    print(f"# Filename: {file}")
    if log_dir:
        log_file = f'{log_dir}/{filename}_log.txt'
    else:
        log_file = None

    validation = ScoringFileValidation(filepath, header=header)

    # Check log
    if log_file:
        _report_errors_to_log_file(validation.errors, log_file)
        _file_validation_state(file, log_file)
    else:
        # Report to stdout
        _report_errors_to_stdout(validation.errors)


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


def _parse_args(args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=_description_text(), epilog=_epilog_text(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-t", "--type",
                        help=f"Type of validator: {' or '.join(val_types)}", metavar='VALIDATOR_TYPE',
                        default='raw')
    parser.add_argument("-f", "--filename",
                        help="""The path to the polygenic scoring file to be validated 
                        (no need to use the [--dir] option)""",
                        metavar='SCORING_FILE_NAME')
    parser.add_argument("-d", "--dir",
                        help="""The name of the directory containing the files that need to processed 
                        (no need to use the [-f] option""")
    parser.add_argument("-D", "--score_dir",
                        help="""<Optional> The name of the directory containing the formatted scoring files
                         to compare with harmonized scoring files (only if "--type hm_pos" is used)""")
    parser.add_argument("-l", "--log_dir",
                        help='<Optional> The name of the log directory where the log file(s) will be stored',
                        required=False)
    parser.add_argument("-c", "--check_filename",
                        help="""<Optional> Check that the file name match the PGS Catalog nomenclature. 
                        Only if using with --type [formatted|hm_pos]""",
                        required=False,
                        action='store_true')
    return parser.parse_args(args)


if __name__ == '__main__':
    run()
