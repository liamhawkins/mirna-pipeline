# Specify company that did the sequencing, TORONTO or BC is supported
COMPANY='BC'

# Choose which directory you want to store the processed files in
ANALYSIS_DIR='/home/student/Desktop/pypipeline/test'

# Specify directory that contains the raw read (.fastq) files
RAW_FILES_DIR='/home/student/Desktop/nmranalysis/raw_files'

# Choose the prefix of files created by R script
PREFIX='NMRBRN'

# Specify files containing TORONTO and BC adapters
TORONTO_ADAPTERS = '/home/student/Desktop/pypipeline/reference_sequences/toronto_adapters.fasta'
BC_ADAPTERS = '/home/student/Desktop/pypipeline/reference_sequences/bc_adapters.fasta'

# Specify files containing Negative, Mature, and Hairpin reference sequences
NEG_FILE='/home/student/Desktop/pipeline/reference_sequences/neg_ref.fasta'
MATURE_FILE='/home/student/Desktop/pipeline/reference_sequences/mature_miRBase22_DNA.fa'
HP_FILE='/home/student/Desktop/pipeline/reference_sequences/hairpin_miRBase22_DNA.fa'

# Specify file containing file names and their corresponding experimental conditions
CONDITION_FILE='/home/student/Desktop/pipeline/reference_sequences/nmr_condition.csv'

# Specify name of control condition
CONTROL_NAME='Control'
