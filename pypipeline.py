#!/usr/bin/python3

import os
from datetime import datetime
from subprocess import Popen
import subprocess
from prompt_toolkit import print_formatted_text, HTML
try:
    from config import COMPANY, ANALYSIS_DIR, SEQ_DIR, RAW_FILES_DIR, NEG_FILE, MATURE_FILE, HP_FILE, CONDITION_FILE, PREFIX
except ImportError:
    print('ERROR importing config.py file. Check to make sure format is correct')
    exit(1)


class WrongFileTypeError(Exception):
    pass


LOG_FILE = os.path.join(ANALYSIS_DIR, datetime.now().strftime('%d-%m-%y_%H:%M') + '.log')
PIPELINE = 'MicroRNA-seq pypipeline'
FASTQC_DIR = os.path.join(ANALYSIS_DIR, 'fastqc/')
TRIMMED_DIR = os.path.join(ANALYSIS_DIR, 'trimmed/')
ADAPTER_DIR = os.path.join(ANALYSIS_DIR, 'adapters/')
BOWTIE_DIR = os.path.join('/home/student/Desktop/pypipeline', 'bowtie_index/')
NEG_IND_DIR = os.path.join(BOWTIE_DIR, 'neg_ref')
HP_IND_DIR = os.path.join(BOWTIE_DIR, 'hp_ref')
MATURE_IND_DIR = os.path.join(BOWTIE_DIR, 'mature_ref')
FLT_DIR = os.path.join(ANALYSIS_DIR, 'filtered/')
MATURE_DIR = os.path.join(ANALYSIS_DIR, 'matures/')
MATURE_ALIGNED_DIR = os.path.join(MATURE_DIR, 'aligned/')
MATURE_UNALIGNED_DIR = os.path.join(MATURE_DIR, 'unaligned/')
HP_DIR = os.path.join(ANALYSIS_DIR, 'hairpins/')
HP_ALIGNED_DIR = os.path.join(HP_DIR, 'aligned/')
READ_COUNT_DIR = os.path.join(ANALYSIS_DIR, 'read_counts/')
MATURE_READ_COUNT_DIR = os.path.join(READ_COUNT_DIR, 'mature/')
HAIRPIN_READ_COUNT_DIR = os.path.join(READ_COUNT_DIR, 'hairpin/')
R_DATA_FILE = PREFIX + '.Rdata'
CLUS_HEATMAP_FILE = PREFIX + '_clus_geat_map'
VOLCANO_FILE = PREFIX + '_volcano'
PCA_PDF = PREFIX + '_pca_pdf'
PCA_CSV = PREFIX + '_pca_csv'
# FIRST_FILE_BASENAME
# STRESS_NAME

GOOD = HTML('<green>GOOD</green>')
FILE_ALREADY_EXISTS = HTML('<yellow>FILE ALREADY EXISTS</yellow>')
NOT_BUILT = HTML('<yellow>NOT BUILT</yellow>')
BAD = HTML('<red>BAD</red>')
EXITING = HTML('<red>EXITING</red>')


def run_command(message, command):
    formatted_message = '[{}] '.format(PIPELINE) + message + '... '

    print(formatted_message, end='', flush=True)
    try:
        output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        with open(LOG_FILE, 'a') as f:
            f.write(formatted_message + 'GOOD\n')
            f.write(output.decode('utf-8'))
        print_formatted_text(GOOD)
        del output
    except subprocess.CalledProcessError as exc:
        print_formatted_text(BAD)
        print('ERROR:')
        print(exc.output.decode('utf-8'))
        with open(LOG_FILE, 'a') as f:
            f.write(formatted_message + 'BAD\n')
            f.write(exc.output.decode('utf-8'))
        exit(1)


def log_message(message, command_status=GOOD):
    formatted_message = '[{}] '.format(PIPELINE) + message + '... '
    print_formatted_text(HTML(formatted_message + command_status.value))

    if type(command_status) == str:
        log_text = command_status
    else:
        log_text = command_status.value

    with open(LOG_FILE, 'a') as f:
        f.write(formatted_message + log_text + '\n')


def validate_file(file_):
    if not os.path.isfile(file_):
        log_message('{} does not exist'.format(file_), command_status=EXITING)
        exit(1)


def validate_config():
    message = 'Creating log file {}'.format(LOG_FILE)
    command = 'touch {}'.format(LOG_FILE)
    run_command(message, command)

    print('[{}] Performing config validation... '.format(PIPELINE), end='', flush=True)
    validate_file(NEG_FILE)
    validate_file(MATURE_FILE)
    validate_file(HP_FILE)
    validate_file(CONDITION_FILE)

    if COMPANY.upper() not in ['BC', 'TORONTO']:
        log_message('COMPANY must be "BC" or "TORONTO"', command_status=EXITING)
        exit(1)

    if COMPANY.upper() == 'TORONTO':
        adapters = 'illumina_adapters.fasta'
    else:
        adapters = 'bc_adapters.fasta'

    if not os.path.isfile(os.path.join(SEQ_DIR, adapters)):
        log_message('{} is missing, exiting'.format(os.path.isfile(os.path.join(SEQ_DIR, adapters))), command_status=EXITING)
        exit(1)

    print_formatted_text(GOOD)

    print('\n')
    for filename in sorted(os.listdir(RAW_FILES_DIR)):
        if filename.endswith('.fastq') or filename.endswith('.fq'):
            print(filename)

    file_validation = input('\n[{}] Are these the files you want to process? (y/n): '.format(PIPELINE))
    if file_validation in ['n', 'N', 'no', 'NO', 'No']:
        print('[{}] Exiting'.format(PIPELINE))
        exit(0)


def check_program(program):
    print('[{}] Checking that {} is installed... '.format(PIPELINE, program), end='', flush=True)

    try:
        Popen([program], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).communicate()
        print_formatted_text(HTML('<green>GOOD</green>'))
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            print('BAD')
            print('[{}] The program {} was not found'.format(PIPELINE, program))
            exit()
        else:
            print('BAD')
            print('[{}] An unknown error occurred when looking for {}'.format(PIPELINE, program))
            raise


def fastqc_check(fastq_file):
    if not os.path.exists(FASTQC_DIR):
        os.makedirs(FASTQC_DIR)

    if not fastq_file.endswith(('.fastq', '.fq')):
        raise WrongFileTypeError('[{}] {} is not a .fastq file'.format(PIPELINE, fastq_file))

    message = 'Performing fastqc check on {}'.format(fastq_file)
    command = 'fastqc -q {} -o {}'.format(fastq_file, FASTQC_DIR)
    run_command(message, command)


def get_basename(path):
    return os.path.splitext(os.path.basename(path))[0].split('.')[0]


def trim_adapters(fastq_file, adapter_file, trim_6=False):
    if not os.path.exists(TRIMMED_DIR):
        os.makedirs(TRIMMED_DIR)

    output_file = os.path.join(TRIMMED_DIR, get_basename(fastq_file) + '.trimmed.fastq')
    temp_file = os.path.join(TRIMMED_DIR, 'temp.fastq')

    message = 'Trimming adapters from {}'.format(fastq_file)
    command = 'cutadapt -q 20 -m 10 -j 18 -b file:{0} {1} -o {2}'.format(adapter_file, fastq_file, temp_file)
    if os.path.exists(output_file):
        log_message(message, command_status=FILE_ALREADY_EXISTS)
    else:
        run_command(message, command)

        if trim_6:
            message = 'Trimming 6 nucleotides from {}'.format(temp_file)
            command = 'cutadapt -u 6 -j 18 {0} -o {1}'.format(temp_file, output_file)
            run_command(message, command)
            os.remove(temp_file)
        else:
            os.rename(temp_file, output_file)

    return output_file


def index_is_built(ind_prefix, index_name):
    try:
        os.listdir(os.path.dirname(ind_prefix))
    except FileNotFoundError:
        log_message('Checking if {} index is built'.format(index_name), command_status=NOT_BUILT)
        return False

    for filename in os.listdir(ind_prefix):
        if not filename.startswith(os.path.basename(ind_prefix)) or not filename.endswith('.ebwt'):
            log_message('Checking if {} index is built'.format(index_name), command_status=NOT_BUILT)
            return False

    log_message('Checking if {} index is built'.format(index_name))
    return True


def build_index(index_dir, index_name):
    if not index_is_built(index_dir, index_name):
        if not os.path.exists(index_dir):
            os.makedirs(index_dir)

        message = 'Building negative index'
        command = 'bowtie-build {} {}'.format(NEG_FILE, os.path.join(index_dir, os.path.basename(index_dir)))
        run_command(message, command)


def filter_out_neg(trimmed_file):
    os.makedirs(FLT_DIR, exist_ok=True)

    basename = get_basename(trimmed_file)
    negative_index = os.path.join(NEG_IND_DIR, os.path.basename(NEG_IND_DIR))
    output_file = os.path.join(FLT_DIR, basename + '.filtered.fastq')

    message = 'Filtering negative RNA species from {}'.format(trimmed_file)
    command = 'bowtie -p 18 -q {} {} --un {}'.format(negative_index, trimmed_file, output_file)
    if os.path.exists(output_file):
        log_message(message, command_status=FILE_ALREADY_EXISTS)
    else:
        run_command(message, command)

    return output_file


def align_mature(filtered_file):
    os.makedirs(MATURE_ALIGNED_DIR, exist_ok=True)
    os.makedirs(MATURE_UNALIGNED_DIR, exist_ok=True)

    basename = get_basename(filtered_file)
    mature_index = os.path.join(MATURE_IND_DIR, os.path.basename(MATURE_IND_DIR))
    aligned_sam = os.path.join(MATURE_ALIGNED_DIR, basename + '.aligned.sam')
    aligned_bam = os.path.join(MATURE_ALIGNED_DIR, basename + '.aligned.bam')
    unaligned_reads = os.path.join(MATURE_UNALIGNED_DIR, basename + '.unaligned.fastq')

    message = 'Aligning {} to mature index'.format(filtered_file)
    command = 'bowtie -p 18 -q -l 20 -n 0 -v 2 -a -S --best --strata {} {} --al -S {} --un {}'.format(mature_index, filtered_file, aligned_sam, unaligned_reads)
    if os.path.exists(aligned_sam) and os.path.exists(unaligned_reads):
        log_message(message, command_status=FILE_ALREADY_EXISTS)
    else:
        run_command(message, command)

    message = 'Converting SAM to BAM: {} to {}'.format(aligned_sam, aligned_bam)
    command = 'samtools view -S -b {} > {}'.format(aligned_sam, aligned_bam)
    if os.path.exists(aligned_bam):
        log_message(message, command_status=FILE_ALREADY_EXISTS)
    else:
        run_command(message, command)

    return aligned_bam, unaligned_reads


def align_hairpins(unaligned_reads):
    os.makedirs(HP_ALIGNED_DIR, exist_ok=True)

    basename = get_basename(unaligned_reads)
    hairpin_index = os.path.join(HP_IND_DIR, os.path.basename(HP_IND_DIR))
    aligned_sam = os.path.join(MATURE_ALIGNED_DIR, basename + '.aligned.sam')
    aligned_bam = os.path.join(MATURE_ALIGNED_DIR, basename + '.aligned.bam')

    message = 'Aligned {} to hairpin index'.format(unaligned_reads)
    command = 'bowtie -p 18 -q -l 20 -n 0 -v 2 -a -S --best --strata {} {} --al -S {}'.format(hairpin_index, unaligned_reads, aligned_sam)
    if os.path.exists(aligned_sam):
        log_message(message, command_status=FILE_ALREADY_EXISTS)
    else:
        run_command(message, command)

    message = 'Converting SAM to BAM: {} to {}'.format(aligned_sam, aligned_bam)
    command = 'samtools view -S -b {} > {}'.format(aligned_sam, aligned_bam)
    if os.path.exists(aligned_bam):
        log_message(message, command_status=FILE_ALREADY_EXISTS)
    else:
        run_command(message, command)

    return aligned_bam


def get_mature_read_counts(aligned_bam):
    os.makedirs(MATURE_READ_COUNT_DIR, exist_ok=True)

    basename = get_basename(aligned_bam)
    sorted_file = os.path.join(MATURE_READ_COUNT_DIR, basename + '.sorted')
    sorted_file_bam = os.path.join(MATURE_READ_COUNT_DIR, basename + '.sorted.bam')
    readcount = os.path.join(MATURE_READ_COUNT_DIR, basename + 'read_count.txt')

    message = 'Sorting {}'.format(aligned_bam)
    command = 'samtools sort -n {} {}'.format(aligned_bam, sorted_file)
    if os.path.exists(sorted_file_bam):
        log_message(message, command_status=FILE_ALREADY_EXISTS)
    else:
        run_command(message, command)

    message = 'Generating read count file from {}'.format(sorted_file_bam)
    command = "samtools view {sorted_file_bam} | awk '{print $3}' | sort | uniq -c | sort -nr > {readcount_file}".format(sorted_file_bam=sorted_file_bam, readcount_file=readcount)
    if os.path.exists(readcount):
        log_message(message, command_status=FILE_ALREADY_EXISTS)
    else:
        run_command(message, command)

    return readcount


def get_hp_read_counts(aligned_bam):
    os.makedirs(HAIRPIN_READ_COUNT_DIR, exist_ok=True)

    basename = get_basename(aligned_bam)
    sorted_file = os.path.join(HAIRPIN_READ_COUNT_DIR, basename + '.sorted')
    sorted_file_bam = os.path.join(HAIRPIN_READ_COUNT_DIR, basename + '.sorted.bam')
    readcount = os.path.join(HAIRPIN_READ_COUNT_DIR, basename + 'read_count.txt')

    message = 'Sorting {}'.format(aligned_bam)
    command = 'samtools sort -n {} {}'.format(aligned_bam, sorted_file)
    if os.path.exists(sorted_file_bam):
        log_message(message, command_status=FILE_ALREADY_EXISTS)
    else:
        run_command(message, command)

    message = 'Generating read count file from {}'.format(sorted_file_bam)
    command = "samtools view {sorted_file_bam} | awk '{print $3}' | sort | uniq -c | sort -nr > {readcount_file}".format(sorted_file_bam=sorted_file_bam, readcount_file=readcount)
    if os.path.exists(readcount):
        log_message(message, command_status=FILE_ALREADY_EXISTS)
    else:
        run_command(message, command)

    return readcount


def copy_read_counts():
    # TODO:
    pass


def run_r_analysis():
    # TODO:
    pass


if __name__ == '__main__':
    # Validate setup
    check_program('fastqc')
    check_program('fastq-mcf')
    check_program('cutadapt')
    check_program('bowtie-build')
    check_program('bowtie')
    check_program('samtools')
    check_program('Rscript')
    validate_config()

    if not os.path.exists(ANALYSIS_DIR):
        os.makedirs(ANALYSIS_DIR)

    # Set company specific variables
    if COMPANY == 'BC':
        trim_6 = True
        adapter_file = os.path.join(SEQ_DIR, 'bc_adapters.fasta')
    else:
        trim_6 = False
        adapter_file = os.path.join(SEQ_DIR, 'illumina_adapters.fasta')

    build_index(NEG_IND_DIR, 'negative')
    build_index(MATURE_IND_DIR, 'mature')
    build_index(HP_IND_DIR, 'hairpin')

    # Process files one at a time
    for filename in sorted(os.listdir(RAW_FILES_DIR)):
        fastq_file = os.path.join(RAW_FILES_DIR, filename)

        # fastqc_check(fastq_file=fastq_file)
        trimmed_file = trim_adapters(fastq_file=fastq_file, adapter_file=adapter_file, trim_6=trim_6)
        filtered_file = filter_out_neg(trimmed_file=trimmed_file)
        mature_aligned_bam, unaligned_reads = align_mature(filtered_file)
        hairpin_aligned_bam = align_hairpins(unaligned_reads)
        mature_readcount = get_mature_read_counts(mature_aligned_bam)
        hairpin_readcount = get_hp_read_counts(hairpin_aligned_bam)
        exit()
