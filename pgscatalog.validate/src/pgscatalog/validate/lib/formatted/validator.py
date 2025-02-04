import gzip
import re
from pandas_schema import Schema
from ..schemas import *
from ..validator_base import *

'''
PGS Catalog Harmonized file validator
- using pandas_schema https://github.com/TMiguelT/PandasSchema
'''

class ValidatorFormatted(ValidatorBase):

    def __init__(self, file, score_dir=None, logfile="VALIDATE.log", error_limit=0):
        super().__init__(file, score_dir, logfile, error_limit)
        self.score_dir=None
        self.meta_format = FORMATTED_META_GENERIC
        self.schema_validators = FORMATTED_VALIDATORS
        self.valid_cols = VALID_COLS_FORMATTED
        self.valid_type = VALID_TYPE_FORMATTED
        self.setup_field_validation()

    
    def extract_specific_metadata(self,line):
        ''' Extract some of the metadata. '''
        match_variants_number = re.search(r'#variants_number=(\d+)', line)
        if match_variants_number:
            self.variants_number = int(match_variants_number.group(1))


    def get_and_check_variants_number(self):
        ''' Verify that the number of variant lines corresponds to the number of variants in the headers '''
        variant_lines = 0
        
        with gzip.open( self.file, 'rb') as f:
            line_number = 0
            for line in f:
                line_number += 1
                line = line.decode('utf-8').rstrip()
                if line.startswith('#'):
                    match_variants_number = re.search(r'#variants_number=(\d+)', line)
                    if match_variants_number:
                        self.variants_number = int(match_variants_number.group(1))
                else:
                    variant_lines += 1
                    if re.search(r'\w+', line): # Line not empty
                        cols = line.split(self.sep)
                        has_trailing_spaces = self.check_leading_trailing_spaces(cols,line_number)
                        if has_trailing_spaces:
                            self.global_errors += 1
                    else:
                        self.logger.error(f'- Line {line_number} is empty')
                        self.global_errors += 1
                            
        if self.variants_number:
            variant_lines -= 1 # Remove the header line from the count
            if self.variants_number != variant_lines:
                self.logger.error(f'- The number of variants lines in the file ({variant_lines}) and the number of variants declared in the headers ({self.variants_number}) are different')
                self.global_errors += 1
        else:
            self.logger.error("- Can't retrieve the number of variants from the headers")
            self.global_errors += 1


    def detect_duplicated_rows(self,dataframe_chunk):
        ''' Detect duplicated rows in the scoring file. '''
        # Columns of interest to compare the different rows
        cols_sel = []
        for col in ['rsID','chr_name','chr_position','effect_allele','other_allele']:
            if col in self.cols_to_validate:
                cols_sel.append(col)

        duplicate_status = dataframe_chunk.duplicated(cols_sel)
        if any(duplicate_status):
            duplicated_rows = dataframe_chunk[duplicate_status]
            self.logger.error(f'Duplicated row(s) found: {len(duplicated_rows.index)}\n\t-> {duplicated_rows.to_string(header=False,index=False)}')
            self.global_errors += 1
            for index in duplicated_rows.index:
                self.bad_rows.append(index)


    def validate_data(self) -> bool:
        ''' Validate the file: data format and data content '''
        self.logger.info("Validating data...")
        if not self.open_file_and_check_for_squareness():
            self.logger.error("Please fix the table. Some rows have different numbers of columns to the header")
            self.logger.info("Rows with different numbers of columns to the header are not validated")
        # Check the consitence between the declared variants number and the actual number of variants in the file
        self.get_and_check_variants_number()

        for chunk in self.df_iterator(self.file):
            dataframe_to_validate = chunk[self.cols_to_read]
            dataframe_to_validate.columns = self.cols_to_validate # sets the headers to standard format if neeeded

            # Detect duplicated rows
            self.detect_duplicated_rows(dataframe_to_validate)

            # validate the snp column if present
            if SNP_DSET in self.header:
                sub_schema = FORMATTED_VALIDATORS_SNP
                if CHR_DSET and BP_DSET in self.header:
                    sub_schema = FORMATTED_VALIDATORS_SNP_EMPTY
                self.validate_schema(sub_schema,dataframe_to_validate)

            if CHR_DSET and BP_DSET in self.header:
                self.validate_schema(FORMATTED_VALIDATORS_POS, dataframe_to_validate)

            if OR_DSET in self.header:
                self.validate_schema(FORMATTED_VALIDATORS_OR,dataframe_to_validate)

            if HR_DSET in self.header:
                self.validate_schema(FORMATTED_VALIDATORS_HR,dataframe_to_validate)

            self.process_errors()
            if len(self.bad_rows) >= self.error_limit:
                break
        if not self.bad_rows and not self.global_errors:
            if self.is_file_valid():
                self.logger.info("File is valid")
            else:
                self.logger.info("File is invalid")
        else:
            self.logger.info("File is invalid - {} bad rows, limit set to {}".format(len(self.bad_rows), self.error_limit))
            self.set_file_is_invalid()
        return self.is_file_valid()


    def validate_filename(self) -> bool:
        ''' Validate the file name structure. '''
        self.logger.info("Validating file name...")
        filename = self.file.split('/')[-1].split('.')[0]
        is_valid_filename = True
        if not re.match(r'^PGS\d{6}$', filename):
            self.logger.info("Invalid filename: {}".format(self.file))
            self.logger.error("Filename: {} should follow the pattern 'PGSXXXXXX.txt.gz', where the 'X' are the 6 digits of the PGS identifier (e.g. PGS000001)".format(filename))
            is_valid_filename = False
            self.set_file_is_invalid()

        return is_valid_filename


    def validate_headers(self) -> bool:
        ''' Validate the list of column names. '''
        self.logger.info("Validating headers...")
        self.detect_genomebuild_with_rsid()
        required_is_subset = set(STD_COLS_VAR_FORMATTED).issubset(self.header)
        if not required_is_subset:
            # check if everything but snp:
            required_is_subset = set(CHR_COLS_VAR_FORMATTED).issubset(self.header)
            if not required_is_subset:
                required_is_subset = set(SNP_COLS_VAR_FORMATTED).issubset(self.header)
            if not required_is_subset:
                self.logger.error("Required headers: {} are not in the file header: {}".format(STD_COLS_VAR_FORMATTED, self.header))

        # Check if at least one of the effect columns is there
        has_effect_col = 0
        for col in STD_COLS_EFFECT_FORMATTED:
            if set([col]).issubset(self.header):
                has_effect_col = 1
                break
        if not has_effect_col:
            self.logger.error("Required headers: at least one of the columns '{}' must be in the file header: {}".format(STD_COLS_EFFECT_FORMATTED, self.header))
            required_is_subset = None

        if not required_is_subset:
            self.logger.info("Invalid headers...exiting before any further checks")
            self.set_file_is_invalid()

        return required_is_subset


    def detect_genomebuild_with_rsid(self):
        ''' The column "rsID" should always be in the scoring file when the genome build is not reported (i.e. "NR") '''
        self.get_genomebuild()
        if self.genomebuild == 'NR':
            if SNP_DSET not in self.header:
                self.logger.error(f"- The combination: Genome Build = '{self.genomebuild}' & the missing column '{SNP_DSET}' in the header is not allowed as we have to manually guess the genome build.")
                self.global_errors += 1


    def get_genomebuild(self):
        ''' Retrieve the Genome Build from the comments '''
        with gzip.open(self.file, 'rb') as f_in:
            for f_line in f_in:
                line = f_line.decode()
                # Update header
                if line.startswith('#genome_build'):
                    gb = (line.split('='))[1]
                    self.genomebuild = gb.strip()
                    return


##################################################################

def init_validator(file, logfile, score_dir=None) -> ValidatorFormatted:
    validator = ValidatorFormatted(file=file, score_dir=score_dir, logfile=logfile)
    return validator