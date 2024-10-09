import os, glob, re
import argparse
import logging
import textwrap

data_sum = {'valid': [], 'invalid': [], 'other': []}

val_types = ('formatted', 'hm_pos')

logging.basicConfig(level=logging.INFO, format='(%(levelname)s): %(message)s')


def validate_scorefile() -> None:
    global data_sum, score_dir
    args = _parse_args()
    _check_args(args)

    # Check PGS Catalog file name nomenclature
    check_filename = False
    if args.check_filename:
        check_filename = True
    else:
        print("WARNING: the parameter '--check_filename' is not present in the submitted command line, therefore the validation of the scoring file name(s) won't be performed.")

    validator_type = args.t
    files_dir = args.dir
    log_dir = args.log_dir

    ## Select validator class ##
    if validator_type == 'formatted':
        import pgscatalog_utils.validate.formatted.validator as validator_package
    elif validator_type == 'hm_pos':
        import pgscatalog_utils.validate.harmonized_position.validator as validator_package

    ## Run validator ##
    # One file
    if args.f:
        _run_validator(args.f,log_dir,score_dir,validator_package,check_filename,validator_type)
    # Content of the directory
    elif files_dir:
        count_files = 0
        # Browse directory: for each file run validator
        for filepath in sorted(glob.glob(files_dir+"/*.*")):
            _run_validator(filepath,log_dir,score_dir,validator_package,check_filename,validator_type)
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
    '''
    Return the last line of the file
    '''
    fileHandle = open ( file,"r" )
    lineList = fileHandle.readlines()
    fileHandle.close()
    return lineList[-1]


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
        else:#
            print("!! validation process had an issue. Please look at the logs.\n")
            data_sum['other'].append(filename)
    else:
        print("!! validation process had an issue: the log file can't be found")
        data_sum['other'].append(filename)


def _check_args(args: argparse.Namespace) -> None:
    global score_dir

    ## Check parameters ##
    # Type of validator
    if args.t not in val_types:
        print(f"Error: Validator type (option -t) '{args.t}' is not in the list of recognized types: {val_types}.")
        exit(1)
    # Logs dir
    if not os.path.isdir(args.log_dir):
        print(f"Error: Log dir '{args.log_dir}' can't be found!")
        exit(1)
    # File and directory parameters (only one of the '-f' and '--dir' can be used)
    if args.f and args.dir:
        print("Error: you can't use both options [-f] - single scoring file and [--dir] - directory of scoring files. Please use only 1 of these 2 options!")
        exit(1)
    elif not args.f and not args.dir:
        print("Error: you need to provide a scoring file [-f] or a directory of scoring files [--dir]!")
        exit(1)
    elif args.f and not os.path.isfile(args.f):
        print(f"Error: Scoring file '{args.f}' can't be found!")
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
    elif args.t != 'formatted':
        print("WARNING: the parameter '--score_dir' is not present in the submitted command line, therefore the comparison of the number of data rows between the formatted scoring file(s) and the harmonized scoring file(s) won't be performed.")


def _run_validator(filepath: str, log_dir: str, score_dir: str, validator_package: object, check_filename: bool, validator_type: str) -> None:
    ''' Run the file validator '''
    file = os.path.basename(filepath)
    filename = file.split('.')[0]
    print(f"# Filename: {file}")
    log_file = f'{log_dir}/{filename}_log.txt'

    # Run validator
    validator = validator_package.init_validator(filepath,log_file,score_dir)
    if check_filename:
        validator.run_validator()
    else:
        validator.run_validator_skip_check_filename()

    # Check log
    _file_validation_state(file,log_file)


def _description_text() -> str:
    return textwrap.dedent('''\
    Validate a set of scoring files to match the PGS Catalog scoring file formats.
    It can validate:
    - The formatted scoring file format (https://www.pgscatalog.org/downloads/#dl_ftp_scoring)
    - The harmonized (Position) scoring file format (https://www.pgscatalog.org/downloads/#dl_ftp_scoring_hm_pos)
   ''')


def _epilog_text() -> str:
    return textwrap.dedent(f'''\
    You need to specify the type of file format to validate, using the paramter '-t' ({' or '.join(val_types)}).
   ''')


def _parse_args(args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=_description_text(), epilog=_epilog_text(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-t", help=f"Type of validator: {' or '.join(val_types)}", metavar='VALIDATOR_TYPE')
    parser.add_argument("-f", help='The path to the polygenic scoring file to be validated (no need to use the [--dir] option)', metavar='SCORING_FILE_NAME')
    parser.add_argument('--dir', help='The name of the directory containing the files that need to processed (no need to use the [-f] option')
    parser.add_argument('--score_dir', help='<Optional> The name of the directory containing the formatted scoring files to compare with harmonized scoring files')
    parser.add_argument('--log_dir', help='The name of the log directory where the log file(s) will be stored', required=True)
    parser.add_argument('--check_filename', help='<Optional> Check that the file name match the PGS Catalog nomenclature', required=False, action='store_true')
    return parser.parse_args(args)


if __name__ == '__main__':
    validate_scorefile()
