# Output error messages

# Columns
MISSING_COLUMN_RSID_OR_COORD = "Missing column either rsID, or both chr_name and chr_position"
MISSING_COLUMN_EFFECT_OR_DOSAGE_WEIGHT = "Missing column either effect_weight, or {'dosage_0_weight', 'dosage_1_weight', 'dosage_2_weight'}"
UNEXPECTED_COLUMNS = "The following columns are not in the schema: {columns}"

# Values
NON_PRINTABLE_CHAR = "Value contains non printable character(s). Input: {input}"
LEADING_OR_TRAILING_SPACE = "Value contains leading or trailing space(s). Input: {input}"
INVALID_CHR_NAME = "Invalid chromosome name. Input: {input}"
INVALID_VARIANT_TYPE = "Invalid variant type. Input: {input}"
NOT_A_BOOLEAN = "Not a valid boolean value"

# Complex variants
BOTH_DIPLOTYPE_AND_HAPLOTYPE = "Cannot be both diplotype and haplotype"
MISSING_DIPLOTYPE_ALLELE = "Variant is marked as diplotype, although it doesn't show 2 effect alleles ('/' separator)."
EXCLUSION_GROUP_NOT_IN_VAR_DESC = "Exclusion group \"{input}\" is not found in variant_description."
EXCLUSION_GROUPS_NOT_IN_EFFECT_ALLELE = "Exclusion group(s) \"{input}\" are declared in variant_description but cannot be found in effect_allele."
DUPLICATE_EXCLUSION_GROUPS = "Exclusion group \"{input}\" is declared twice"
LOCUS_NAME_NOT_APOE = "locus_name of APOE allele variant is not APOE"
