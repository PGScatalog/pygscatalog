import os, sys, gc
import gzip
import csv
import pathlib
import logging
import re
from typing import List
import pandas as pd
import pandas_schema
import warnings
from .schemas import *

'''
PGS Catalog file validator
- using pandas_schema https://github.com/TMiguelT/PandasSchema
'''


csv.field_size_limit(sys.maxsize)

class ValidatorBase:

    valid_extensions = VALID_FILE_EXTENSIONS
    schema_validators = GENERIC_VALIDATORS
    valid_cols = []
    valid_type = ''
    sep = '\t'

    def __init__(self, file, score_dir=None, logfile="VALIDATE.log", error_limit=0):
        self.file = file
        self.score_dir = score_dir
        self.schema = None
        self.header = []
        self.genomebuild = None
        self.comment_lines_count = 1 # Counting the header line
        self.cols_to_validate = []
        self.cols_to_read = []
        self.bad_rows = []
        self.row_errors = []
        self.errors_seen = {}
        self.logfile = logfile
        self.error_limit = int(error_limit)
        self.is_valid = True

        # Logging variables
        self.logger = logging.getLogger(__name__)
        self.handler = logging.FileHandler(self.logfile, 'w+')
        self.handler.setLevel(logging.INFO)
        self.logger.addHandler(self.handler)
        self.logger.propagate = False

        self.global_errors = 0
        self.variants_number = 0


    def validate_schema(self, schema: dict, dataframe_to_validate: pd.core.frame.DataFrame):
        '''
        Run the pandas_schema validation using the provided Schema and DataFrame
        '''
        self.schema = pandas_schema.Schema([schema[h] for h in self.cols_to_validate])
        with warnings.catch_warnings():
            # Ignore python warningd raised in the pandas_schema code
            warnings.simplefilter('ignore', UserWarning)
            errors = self.schema.validate(dataframe_to_validate)
            self.store_errors(errors)


    def setup_field_validation(self):
        '''
        Fetch the header and build the list of column to check/validate
        '''
        self.header = self.get_header()
        self.cols_to_validate = [h for h in self.header if h in self.valid_cols]
        self.cols_to_read = [h for h in self.header if h in self.valid_cols]


    def get_header(self):
        '''
        Fetch the header (i.e. column names) information from the harmonized scoring file and store the list in a variable
        '''
        first_row = pd.read_csv(self.file, sep=self.sep, comment='#', nrows=1, index_col=False)
        # Check if the column headers have leading and/or trailing spaces
        # The leading/trailing spaces should raise an error during the header validation
        has_trailing_spaces = self.check_leading_trailing_spaces(first_row.columns.values)
        if has_trailing_spaces:
            self.global_errors += 1
        return first_row.columns.values


    def get_genomebuild(self):
        ''' Retrieve the Genome Build from the comments '''
        if self.valid_type == 'hm_pos':
            self.genomebuild = self.get_comments_info('#HmPOS_build')
        else:
            self.genomebuild = self.get_comments_info('#Hm_genome_build')


    def get_pgs_id(self):
        ''' Retrieve the PGS ID from the comments '''
        self.pgs_id = self.get_comments_info('#pgs_id')


    def validate_content(self):
        ''' Validate the file content and verify that the number of variant lines corresponds to the number of variants in the headers '''
        variant_lines_count = 0
        meta_lines_count = 0
        
        with gzip.open( self.file, 'rb') as f:
            line_number = 0
            file_meta = []
            for line in f:
                line_number += 1
                line = line.decode('utf-8').rstrip()
                # Check Metadata
                if line.startswith('#'):
                    self.extract_specific_metadata(line)
                    # Check that we have all the meta information
                    for meta in self.meta_format:
                        if line.startswith(meta):
                            file_meta.append(meta)
                            meta_lines_count += 1
                            break
                
                # Check data
                else:
                    variant_lines_count += 1
                    if re.search(r'\w+', line): # Line not empty
                        cols_content = line.split(self.sep)
                        has_trailing_spaces = self.check_leading_trailing_spaces(cols_content,line_number)
                        if has_trailing_spaces:
                            self.global_errors += 1
                        
                        if line.startswith('rsID') or line.startswith('chr_name'):
                            continue
                        
                        self.validate_line_content(cols_content,variant_lines_count)
                    else:
                        self.logger.error(f'- Line {line_number} is empty')
                        self.global_errors += 1
        
        # Compare the number of metadata lines: read vs expected
        if meta_lines_count != len(self.meta_format):
            self.logger.error(f'- The number of metadata lines [i.e. starting with the "#" character] in the file ({meta_lines_count}) and the expected number of metadata lines ({len(self.meta_format)}) are different')
            diff_list = list(set(self.meta_format).difference(file_meta))
            self.logger.error(f"  > Missing metadata line(s): {', '.join(diff_list)}")
            self.global_errors += 1


    def validate_data(self) -> bool:
        ''' Validate the file: data format and data content '''
        self.logger.info("Validating data...")
        if not self.open_file_and_check_for_squareness():
            self.logger.error("Please fix the table. Some rows have different numbers of columns to the header")
            self.logger.info("Rows with different numbers of columns to the header are not validated")

        # Validate data content and check the consitence between the declared variants number and the actual number of variants in the file
        self.validate_content()
        for chunk in self.df_iterator(self.file):
            dataframe_to_validate = chunk[self.cols_to_read]
            dataframe_to_validate.columns = self.cols_to_validate # sets the headers to standard format if neeeded

            # Schema validation
            self.validate_schema(self.schema_validators,dataframe_to_validate)

            self.process_errors()
            if len(self.bad_rows) >= self.error_limit:
                break

        if not self.bad_rows and not self.global_errors and self.is_valid:
            self.logger.info("File is valid")
        else:
            self.logger.info("File is invalid - {} bad rows, limit set to {}".format(len(self.bad_rows), self.error_limit))
            self.set_file_is_invalid()
        return self.is_valid


    def is_file_valid(self) -> bool:
        ''' Method returning the boolean value: True if the file is valid, False if the file is invalid. '''
        return self.is_valid

    def set_file_is_invalid(self):
        ''' Set the flag "is_valid" to False. '''
        self.is_valid = False


    def process_errors(self):
        ''' Populate the logger error and the list of bad rows with the errors found. '''
        for error in self.row_errors:
            if len(self.bad_rows) < self.error_limit or self.error_limit < 1:
                self.logger.error(error)
                if error.row not in self.bad_rows:
                    self.bad_rows.append(error.row)
        self.row_errors = []


    def store_errors(self, errors: List[pandas_schema.validation_warning.ValidationWarning]):
        ''' Capture the errors found into a temporary structure before being processed. '''
        for error in errors:
            seen = 0
            row_number = error.row
            file_line_number = row_number + self.comment_lines_count + 1 # rows are 0 indexes
            error.row = str(row_number) + " (line "+str(file_line_number)+")"
            col = error.column
            # Avoid duplication as the errors can be detected several times
            if row_number in self.errors_seen.keys():
                if col in self.errors_seen[row_number].keys():
                    seen = 1
                else:
                    self.errors_seen[row_number][col] = 1
            else:
                self.errors_seen[row_number] = { col : 1 }
            if seen == 0:
                self.row_errors.append(error)


    def validate_file_extension(self):
        ''' Check/validate the file name extension. '''
        self.logger.info("Validating file extension...")
        check_exts = [self.check_ext(ext) for ext in self.valid_extensions]
        if not any(check_exts):
            self.valid_ext = False
            self.set_file_is_invalid()
            self.logger.info("Invalid file extension: {}".format(self.file))
            self.logger.error("File extension should be in {}".format(self.valid_extensions))
        else:
            self.valid_ext = True
        return self.valid_ext


    def compare_number_of_rows(self):
        ''' Compare the number of data rows between the harmonized and the formatted scoring files. '''
        # Harmonization file - length
        hm_rows_count = 0
        for chunk in self.df_iterator(self.file):
            hm_rows_count += len(chunk.index)
        gc.collect()

        # Formatted scoring file - length
        scoring_rows_count = 0
        scoring_file = f'{self.score_dir}/{self.pgs_id}.txt.gz'
        if os.path.isfile(scoring_file):
            for score_chunk in self.df_iterator(scoring_file):
                scoring_rows_count += len(score_chunk.index)
            gc.collect()

        comparison_status = True
        if scoring_rows_count == 0:
            self.logger.error(f"Can't find the Scoring file '{scoring_file}' to compare the number of rows with the harmonization file!")
            comparison_status = False
        elif hm_rows_count != scoring_rows_count:
            self.logger.error(f'The number of data rows between the Scoring file ({scoring_rows_count}) and the Harmonization POS file ({hm_rows_count}) are different')
            comparison_status = False
        return comparison_status


    def compare_with_filename(self):
        ''' Check that the filename matches the information present in the file metadata (PGS ID, genome build). '''
        self.logger.info("Comparing filename with metadata...")
        comparison_status = True
        if hasattr(self,'file_genomebuild') and hasattr(self,'file_pgs_id'):
            # Extract some metadata
            self.get_genomebuild()
            self.get_pgs_id()
            # Compare metadata with filename information
            if self.file_genomebuild != self.genomebuild:
                self.logger.error("Build: the genome build in the HmPOS_build header ({}) is different from the one on the filename ({})".format(self.genomebuild,self.file_genomebuild))
                comparison_status = False
            if self.file_pgs_id != self.pgs_id:
                self.logger.error("ID: the PGS ID of the header ({}) is different from the one on the filename ({})".format(self.pgs_id,self.file_pgs_id))
                comparison_status = False
            # Compare number of rows with Scoring file
            if self.score_dir:
                row_comparison_status = self.compare_number_of_rows()
                if row_comparison_status == False:
                    comparison_status = row_comparison_status
            else:
                self.logger.info("Comparison of the number of rows between Harmonized and Scoring file skipped!")
            if not comparison_status:
                self.logger.info("Discrepancies between filename information and metadata: {}".format(self.file))                
                self.set_file_is_invalid()
        return comparison_status


    def df_iterator(self, data_file: str):
        ''' Setup a pandas dataframe iterator. '''
        df = pd.read_csv(data_file,
                         sep=self.sep,
                         dtype=str,
                         comment='#',
                         chunksize=1000000)
        return df


    def check_file_is_square(self, csv_file: str):
        ''' Check that each row has the name number of columns. '''
        square = True
        csv_file.seek(0)
        reader = csv.reader(csv_file, delimiter=self.sep)
        count = 1
        for row in reader:
            if len(row) != 0:
                if row[0].startswith('#'):
                    self.comment_lines_count += 1
                    continue
            if (len(row) != len(self.header)):
                self.logger.error("Length of row {c} is: {l} instead of {h}".format(c=count, l=str(len(row)), h=str(len(self.header))))
                self.logger.error("ROW: "+str(row))
                square = False
            count += 1
        del csv_file
        return square


    def open_file_and_check_for_squareness(self):
        ''' Method to read the file in order to check that each row has the name number of columns. '''
        if pathlib.Path(self.file).suffix in [".gz", ".gzip"]:
             with gzip.open(self.file, 'rt') as f:
                 return self.check_file_is_square(f)
        else:
            with open(self.file) as f:
                 return self.check_file_is_square(f)


    def check_leading_trailing_spaces(self, cols:str, line_number:str = None):
        '''
        Check if the columns have leading and/or trailing spaces.
        The leading/trailing spaces should raise an error during the validation.
        '''
        leading_trailing_spaces = []
        found_trailing_spaces = False
        for idx, col in enumerate(cols):
            if col.startswith(' ') or col.endswith(' '):
                leading_trailing_spaces.append(self.header[idx]+' => |'+str(col)+'|')
        if len(leading_trailing_spaces):
            if line_number:
                line_name = f'line {line_number} has'
            else:
                line_name = 'following headers have'
            self.logger.error("The "+line_name+" leading and/or trailing spaces: "+' ; '.join(leading_trailing_spaces))
            found_trailing_spaces = True
        return found_trailing_spaces


    def check_ext(self, ext:str) -> bool:
        if self.file.endswith(ext):
            return True
        return False


    def check_build_is_legit(self, build:str) -> bool:
        if build in BUILD_LIST:
            return True
        return False

    
    def get_comments_info(self, type:str) -> str:
        ''' Retrieve information from the comments '''
        with gzip.open(self.file, 'rb') as f_in:
            for f_line in f_in:
                line = f_line.decode()
                # Update header
                if line.startswith(type):
                    info = (line.split('='))[1]
                    return info.strip()

    def run_generic_validator(self,check_filename):
        self.logger.propagate = False

        # Check files exist
        if not self.file or not self.logfile:
            self.logger.info("Missing file and/or logfile")
            self.set_file_is_invalid()
        elif self.file and not os.path.exists(self.file):
            self.logger.info("Error: the file '"+self.file+"' can't be found")
            self.set_file_is_invalid()

        # Validate file extension
        self.validate_file_extension()

        # Validate file name nomenclature
        if self.is_file_valid() and check_filename:
            self.validate_filename()

        # Only for harmonized files
        if self.is_file_valid() and type(self).__name__ != 'ValidatorFormatted':
            self.compare_with_filename()

        # Validate column headers
        if self.is_file_valid():
            self.validate_headers()

        # Validate data content
        if self.is_file_valid():
            self.validate_data()

        # Close log handler
        self.logger.removeHandler(self.handler)
        self.handler.close()

    def run_validator(self):
        self.run_generic_validator(True)

    def run_validator_skip_check_filename(self):
        self.run_generic_validator(False)


    def validate_filename(self):
        ''' Validate the file name structure. '''
        print("To be implemented in inherited classes")
        pass


    def validate_headers(self):
        ''' Validate the list of column names. '''
        print("To be implemented in inherited classes")
        pass


    def validate_line_content(self, cols_content:str, var_line_number:int):
        ''' Validate each data row. '''
        print("To be implemented in inherited classes")
        pass


    def extract_specific_metadata(self, line:str):
        ''' Extra method to extract and validate specific data. '''
        print("To be implemented in inherited classes")
        pass

