#!/usr/bin/python3

# TODO: Validate all FASTQ files
import argparse
import csv
import os
import re
import shutil
import subprocess
from configparser import ConfigParser
from datetime import datetime

from prompt_toolkit import HTML, print_formatted_text
from prompt_toolkit.shortcuts import yes_no_dialog


class File:
    def __init__(self, raw_path, analysis_dir):
        """
        The `File` class creates and manages filepaths of a raw reads file and all
        intermediate processing files derived from it. A raw read fastq file is trimmed,
        filtered, aligned etc. during the pipeline and this class makes it easy to track
        the filepaths for each of these steps.

        :param raw_path: filepath of raw reads fastq file
        :type raw_path: str
        :param analysis_dir: pipeline root analysis directory
        :type analysis_dir: str
        """
        self.raw = raw_path
        self.analysis_dir = analysis_dir
        self.basename = self._get_basename(raw_path)
        self.trimmed = self._create_filepath('.trimmed.fastq')
        self.temp = self._create_filepath('.tmp.fastq')
        self.filtered = self._create_filepath('.filtered.fastq')
        self.mature_aligned_sam = self._create_filepath('_MATURE.aligned.sam')
        self.mature_aligned_bam = self._create_filepath('_MATURE.aligned.bam')
        self.unaligned = self._create_filepath('.unaligned.fastq')
        self.hairpin_aligned_sam = self._create_filepath('_HAIRPIN.aligned.sam')
        self.hairpin_aligned_bam = self._create_filepath('_HAIRPIN.aligned.bam')
        self.mature_basename = self._get_basename(self.mature_aligned_sam)
        self.hairpin_basename = self._get_basename(self.hairpin_aligned_sam)
        self.mature_sorted = self._create_filepath('.sorted.bam', file=self.mature_aligned_bam)
        self.hairpin_sorted = self._create_filepath('.sorted.bam', file=self.hairpin_aligned_bam)
        readcount_dir = os.path.join(self.analysis_dir, 'read_counts/')
        os.makedirs(readcount_dir, exist_ok=True)
        self.mature_readcount = self._create_filepath('.read_count.txt', file=self.mature_sorted, dir=readcount_dir)
        self.hairpin_readcount = self._create_filepath('.read_count.txt', file=self.hairpin_sorted, dir=readcount_dir)

    def __repr__(self):
        return '{}({}, {})'.format(self.__class__.__name__, self.raw, self.analysis_dir)

    @staticmethod
    def _get_basename(path):
        """
        Get the basename of a file ex. /some/path/basename.descriptor.txt -> basename

        Filenames in this pipeline have a "basename", an optional descriptor, and a file
        extension all separated by a period.

        :param path: filepath
        :type path: str
        :return: basename of filepath
        :rtype: str
        """
        return os.path.splitext(os.path.basename(path))[0].split('.')[0]

    def _create_filepath(self, suffix, file=None, dir=None):
        """
        Create a new filepath from `self.basename` + `suffix` (or basename of `file` + `suffix`)

        NOTE: This doesn't actually create a file, rather it just forms a string of the filepath

        :param suffix: string used as suffix of new filepath
        :type suffix: str
        :param file: optional file to get basename from
        :type file: str
        :param dir: optional directory to use in new filepath
        :type dir: str
        :return: filepath
        :rtype: str
        """
        if file is None:
            basename = self.basename
        else:
            basename = self._get_basename(file)

        if dir is None:
            dir = self.analysis_dir

        return os.path.join(dir, basename + suffix)

    def remove_intermediates(self):
        """
        Remove all intermediate files generated during the pipeline
        """
        for file_ in [self.trimmed,
                      self.filtered,
                      self.mature_aligned_sam,
                      self.mature_aligned_bam,
                      self.unaligned,
                      self.hairpin_aligned_sam,
                      self.hairpin_aligned_bam,
                      self.mature_sorted,
                      self.hairpin_sorted]:
            try:
                os.remove(file_)
            except FileNotFoundError:
                pass

    def read_counts_exist(self, mature_only=False, hairpin_only=False):
        """
        Check to see if read count file(s) have been created

        :param mature_only: If True, only check if mature readcount file exists
        :type mature_only: bool
        :param hairpin_only: If True, only check is hairpin readcount file exists
        :type hairpin_only: bool
        :return: boolean result of whether readcount file(s) exist
        :rtype: bool
        """
        if mature_only and hairpin_only:
            raise ValueError('mature_only and hairpin_only cannot both be True')
        elif mature_only:
            return os.path.isfile(self.mature_readcount)
        elif hairpin_only:
            return os.path.isfile(self.hairpin_readcount)
        else:
            # Defaults to check if both mature and hairpin readcount files exist
            return os.path.isfile(self.mature_readcount) and os.path.isfile(self.hairpin_readcount)

    def change_read_count_dir(self, dir_):
        """
        Change the directory of the readcount filepath string

        This method is useful when readcount files already exist in a specific directory

        :param dir_:  new directory to use in reacount filepaths
        :type dir_: str
        """
        self.mature_readcount = os.path.join(dir_, os.path.basename(self.mature_readcount))
        self.hairpin_readcount = os.path.join(dir_, os.path.basename(self.hairpin_readcount))


class PyPipeline:
    def __init__(self, config_file, no_prompts=False, no_fastqc=None, delete=None, no_analysis=False, read_count_dir=None):
        """
        A microRNA-seq processing and analysis pipeline

        This class contains all the methods to process and analyse raw reads fastq files and output microRNA read counts
        and publication quality figures using the RBioMir suite of R packages.

        :param config_file: filepath of config file withall necessary information about location of pipeline resources
        :type config_file: str
        :param no_prompts: silence all prompts
        :type no_prompts: bool
        :param no_fastqc: do not perform fastqc analysis on raw reads
        :type no_fastqc: bool
        :param delete:  delete intermediate processing files after completion of pipeline
        :type delete: bool
        :param no_analysis: do not perform R-based analysis of microRNA read counts
        :type no_analysis: bool
        :param read_count_dir: if supplied, processing of raw reads is skipped and analysis is performed on read counts
        :type read_count_dir: str
        """
        self.timestamp = lambda: datetime.now().strftime("%Y-%m-%d %H:%M")
        self.no_prompts = no_prompts
        self.no_fastqc = no_fastqc
        self.delete = delete
        self.no_analysis = no_analysis
        self.analysis_only = bool(read_count_dir)

        self.trim_summary = []
        self.filtering_bowtie_summary = []
        self.mature_bowtie_summary = []
        self.hairpin_bowtie_summary = []

        # Read in Config File
        config = ConfigParser()
        config.optionxform = str  # Have to replace optionxform function with str to preserve case sensitivity
        config.read(config_file)
        self.sample_conditions = {k: v for (k, v) in config['sample_conditions'].items()}
        config = config[config.sections()[0]]
        self.company = config['COMPANY']
        self.command_log = config.get('COMMAND_LOG', None)
        self.analysis_dir = config['ANALYSIS_DIR']
        self.raw_files_dir = config['RAW_FILES_DIR']
        self.toronto_adapters = config['TORONTO_ADAPTERS']
        self.bc_adapters = config['BC_ADAPTERS']
        self.negative_references = config['NEGATIVE_REFERENCE_FILE']
        self.mature_references = config['MATURE_REFERENCE_FILE']
        self.hairpin_references = config['HAIRPIN_REFERENCE_FILE']
        self.bowtie_dir = config['BOWTIE_DIR']
        self.kegg_id_file = config['KEGG_ID_FILE']
        self.go_bp_id_file = config['GO_BP_ID_FILE']
        self.go_mf_id_file = config['GO_MF_ID_FILE']
        self.go_cc_id_file = config['GO_CC_ID_FILE']
        if not self.no_analysis:
            self.rpipeline = config['R_SCRIPT']

        # Set up directories and filepaths
        self.log_file = os.path.join(self.analysis_dir, datetime.now().strftime('%d-%m-%y_%H:%M') + '.log')
        self.summary_file = os.path.join(self.analysis_dir, 'summary.txt')
        self.fastqc_dir = os.path.join(self.analysis_dir, 'fastqc/')
        self.negative_index_dir = os.path.join(self.bowtie_dir, 'neg_ref')
        self.hairpin_index_dir = os.path.join(self.bowtie_dir, 'hp_ref')
        self.mature_index_dir = os.path.join(self.bowtie_dir, 'mature_ref')
        self.figures = os.path.join(self.analysis_dir, 'figures/')
        self.mirna_targets_dir = os.path.join(self.figures, 'mirna_targets/')
        self.conditions_file = os.path.join(self.figures, 'conditions.csv')

        # Formatted strings
        self.GOOD = HTML('<green>GOOD</green>')
        self.FILE_ALREADY_EXISTS = HTML('<yellow>FILE ALREADY EXISTS</yellow>')
        self.NOT_BUILT = HTML('<yellow>NOT BUILT</yellow>')
        self.BAD = HTML('<red>BAD</red>')
        self.EXITING = HTML('<red>EXITING</red>')
        self.NONE = HTML('')
        self.F_PIPELINE = lambda: '<teal>{}</teal>'.format(self.timestamp())

        # Create log file
        os.makedirs(self.analysis_dir, exist_ok=True)
        self._create_log_file()

        # Create File object for each raw reads fastq file
        self.files = []
        for dirpath, _, filenames in os.walk(self.raw_files_dir):
            for f in sorted([f for f in filenames if f.endswith(('.fastq', '.fq'))]):
                abs_path = os.path.abspath(os.path.join(dirpath, f))
                self.files.append(File(abs_path, self.analysis_dir))
        if self.analysis_only:
            for file in self.files:
                file.change_read_count_dir(read_count_dir)

        # Set up config-dependent adapter variables
        self.adapters = None
        self.trim_6 = None
        self._validate_config()

    def _run_command(self, command, message, log_stdout=False, log_stderr=False):
        """
        Run CLI command, log and show message in terminal and whether command completed without error

        This method is used to run all CLI commands of the pipeline, giving an easy way to take the bash
        commands you would run and log them to the command log, log their results/output to the log file
        and print a message to the terminal

        :param command: bash command to run
        :type command: str
        :param message: message to print to terminal telling user what command/process is running
        :type message: str
        :param log_stdout: whether to log stdout produced by command
        :type log_stdout: bool
        :param log_stderr: whether to log stderr produced by commond
        :type log_stderr: bool
        """
        # Log command to command log
        if self.command_log:
            with open(self.command_log, 'a') as f:
                f.write(command + '\n')

        # Print message to screen and log message
        formatted_message = '[{}] {}...'.format(self.F_PIPELINE(), message)
        unformatted_message = '[{}] {}...'.format(self.timestamp(), message)
        print_formatted_text(HTML(formatted_message), end='', flush=True)
        with open(self.log_file, 'a') as f:
            f.write(unformatted_message + '\n')

        # Destinations added to command to direct stdout, stderr, or neither to the log
        if log_stdout and log_stderr:
            command += ' &>> {}'.format(self.log_file)
        elif log_stdout:
            command += ' >> {}'.format(self.log_file)
        elif log_stderr:
            command += ' 2>> {}'.format(self.log_file)

        # Call command, if error print to screen, log error, and exit with code 1
        try:
            subprocess.call(command, shell=True, stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as exc:
            print_formatted_text(self.BAD)
            print('ERROR:')
            print(exc.output.decode('utf-8'))
            with open(self.log_file, 'a') as f:
                f.write(exc.output.decode('utf-8'))
            exit(1)
        else:
            print_formatted_text(self.GOOD)

    def _log_message(self, message, command_status=None, **kwargs):
        """
        Log message and print to terminal

        :param message: Message to log and print to terminal
        :type message: str
        :param command_status: Optional command status (default=GOOD)
        :type command_status: HTML
        :param kwargs: kwargs passed to print_formatted_text function
        """
        # Default command status is GOOD
        if command_status is None:
            command_status = self.GOOD

        # Print message to screen and log message
        formatted_message = '[{}] {}...'.format(self.F_PIPELINE(), message)
        unformatted_message = '[{}] {}...'.format(self.timestamp(), message)
        print_formatted_text(HTML(formatted_message + command_status.value), **kwargs)
        with open(self.log_file, 'a') as f:
            f.write(unformatted_message + '\n')

    def _create_log_file(self):
        """
        Create log file using touch Unix command
        """
        message = 'Creating log file {}'.format(os.path.basename(self.log_file))
        command = 'touch {}'.format(self.log_file)
        self._run_command(command, message)

    def _create_summary_file(self):
        """
        Create summary file using touch Unix command
        """
        message = 'Creating summary file - {}'.format(os.path.basename(self.summary_file))
        command = 'touch {}'.format(self.summary_file)
        self._run_command(command, message)

    def _validate_file(self, file_):
        """
        Checks to see if file exists, if not exits with code 1

        :param file_: filepath to file
        :type file_: str
        """
        if not os.path.isfile(file_):
            self._log_message('{} does not exist'.format(file_), command_status=self.EXITING)
            exit(1)

    def _validate_sample_conditions(self):
        """
        Validate that all files are present in config file and are specified as `control` or `stress`,
        if not exits with code 1
        """
        for file in self.files:
            if file.basename not in self.sample_conditions.keys():
                self._log_message('Cannot find sample condition in config file: {}'.format(file.basename), command_status=self.EXITING)
                exit(1)

        if any([condition not in ['control', 'stress'] for condition in self.sample_conditions.values()]):
            self._log_message('Sample conditions can only be "control" or "stress" at this time', command_status=self.EXITING)
            exit(1)

    def _validate_config(self):
        """
        Validates config file, sets config-dependant adapter variables, prompts user with files that will be processed
        """
        self._log_message('Performing config validation', command_status=self.NONE, end='', flush=True)

        # Set config-dependant adapter variables, exits with code 1 if not BC or TORONTO
        if self.company.upper() == 'TORONTO':
            self.adapters = self.toronto_adapters
            self.trim_6 = False
        elif self.company.upper() == 'BC':
            self.adapters = self.bc_adapters
            self.trim_6 = True
        else:
            self._log_message('COMPANY must be "BC" or "TORONTO"', command_status=self.EXITING)
            exit(1)

        # Validates resource files specified in config
        self._validate_file(self.adapters)
        self._validate_file(self.negative_references)
        self._validate_file(self.mature_references)
        self._validate_file(self.hairpin_references)
        self._validate_file(self.kegg_id_file)
        self._validate_file(self.go_bp_id_file)
        self._validate_file(self.go_mf_id_file)
        self._validate_file(self.go_cc_id_file)
        if not self.no_analysis:
            self._validate_file(self.rpipeline)

        # Unless --no-prompts flag used, prompts user with list of found files
        if not self.no_prompts:
            files = '\n'.join([file for file in sorted(os.listdir(self.raw_files_dir)) if file.endswith('.fastq') or file.endswith('.fq')])
            continue_ = yes_no_dialog(title='File check', text='Are these the files you want to process?\n\n' + files)
            if not continue_:
                exit(0)

        self._validate_sample_conditions()

        print_formatted_text(self.GOOD)

    def _check_program(self, program):
        """
        Check to see if CLI program is installed

        If program is not found exists with code 1 and log error

        :param program: name of program (CLI command used to call program)
        :type program: str
        """
        self._log_message('Checking that {} is installed'.format(program), command_status=self.NONE, end='', flush=True)

        try:
            # Tries to communicate with program, exit with code 1 if not found
            subprocess.Popen([program], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).communicate()
            print_formatted_text(self.GOOD)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                self._log_message('The program {} was not found'.format(program), command_status=self.BAD)
                exit(1)
            else:
                self._log_message('An unknown error occurred when looking for {}'.format(program), command_status=self.BAD)
                raise

    def _index_is_built(self, dir_, name):
        """
        Check to see if bowtie index is built

        Given the directory of a specific bowtie index, checks to see that all files in directory have .ebwt extension

        :param dir_: directory of bowtie index
        :type dir_: str
        :param name: name of index to display in terminal and log
        :type name: str
        :return: True if index is built else False
        :rtype: bool
        """
        try:
            os.listdir(os.path.dirname(dir_))
        except FileNotFoundError:
            self._log_message('Checking if {} index is built'.format(name), command_status=self.NOT_BUILT)
            return False

        for filename in os.listdir(dir_):
            if not filename.startswith(os.path.basename(dir_)) or not filename.endswith('.ebwt'):
                self._log_message('Checking if {} index is built'.format(name), command_status=self.NOT_BUILT)
                return False

        self._log_message('Checking if {} index is built'.format(name))
        return True

    def _build_index(self, sequences, dir_, name):
        """
        Build bowtie index using bowtie-build CLI command

        :param sequences: filepath of fasta file containing sequences to build index from
        :type sequences: str
        :param dir_: directory to build bowtie index
        :type dir_: str
        :param name: name of index to use when printing and logging message
        :type name: str
        """
        if not self._index_is_built(dir_, name):
            os.makedirs(dir_, exist_ok=True)

            message = 'Building {} index'.format(name)
            command = 'bowtie-build {} {}'.format(sequences, os.path.join(dir_, os.path.basename(dir_)))
            self._run_command(command, message)

    def _fastqc_check(self, file):
        """
        Perform fastqc on file

        :param file: File object representing sample to be checked
        :type file: File
        """
        os.makedirs(self.fastqc_dir, exist_ok=True)
        message = 'Performing fastqc check on {}'.format(file.basename)
        command = 'fastqc -q {} -o {}'.format(file.raw, self.fastqc_dir)
        self._run_command(command, message)

    def _trim_adapters(self, file):
        """
        Perform adapter triming on file using cutadapt program

        If `self.trim_6` is True it trims adapter then next 6 nucleotides. This should be used when random 6 base
        barcode adapters are used.
        NOTE: only performed if trimmed file does not already exist

        :param file: File object representing sample to be trimmed
        :type file: File
        """
        message = '{}: Trimming adapters'.format(file.basename)
        command = 'cutadapt -q 20 -m 10 -j 18 -b file:{0} {1} -o {2}'.format(self.adapters, file.raw, file.temp)
        if os.path.exists(file.trimmed):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(command, message, log_stdout=True)
            self._get_trim_summary(self.log_file)

            if self.trim_6:
                message = '{}: Trimming 6 nucleotides'.format(file.basename)
                command = 'cutadapt -u 6 -j 18 {0} -o {1}'.format(file.temp, file.trimmed)
                self._run_command(command, message, log_stdout=True)
                os.remove(file.temp)
            else:
                os.rename(file.temp, file.trimmed)

    def _filter_out_neg(self, file):
        """
        Filter out negative reference sequences from trimmed fastq file

        For this pipeline, negative reference sequences are other small non-microRNA RNAs such as rRNA, tRNA, piRNA etc.
        NOTE: only performed if filtered file does not already exist

        :param file: File object representing sample to be filtered
        :type file: File
        """
        negative_index = os.path.join(self.negative_index_dir, os.path.basename(self.negative_index_dir))

        message = '{}: Filtering negative RNA species'.format(file.basename)
        command = 'bowtie -p 18 -q {} {} --un {}'.format(negative_index, file.trimmed, file.filtered)
        if os.path.exists(file.filtered):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(command, message, log_stderr=True)
            self._get_bowtie_summary(self.log_file, 'filtering')

    def _align_reads(self, file):
        """
        Align reads from filtered fastq file to mature and hairpin bowtie indexes and convert resulting SAM to BAM file

        Alignment to bowtie indexes is performed using bowtie, then samtools is used to convert SAM file to BAM file
        NOTE: only performed if aligned files do not already exist

        :param file: File object representing sample to be aligned
        :type file: File
        """
        mature_index = os.path.join(self.mature_index_dir, os.path.basename(self.mature_index_dir))
        hairpin_index = os.path.join(self.hairpin_index_dir, os.path.basename(self.hairpin_index_dir))

        message = '{}: Aligning to mature index'.format(file.basename)
        command = 'bowtie -p 18 -q -l 20 -n 0 -v 2 -a -S --best --strata {} {} --al -S {} --un {}'.format(
            mature_index, file.filtered, file.mature_aligned_sam, file.unaligned)
        if os.path.exists(file.mature_aligned_sam) and os.path.exists(file.unaligned):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(command, message, log_stderr=True)
            self._get_bowtie_summary(self.log_file, 'mature')

        message = '{}: Converting SAM to BAM'.format(file.basename)
        command = 'samtools view -S -b {} > {}'.format(file.mature_aligned_sam, file.mature_aligned_bam)
        if os.path.exists(file.mature_aligned_bam):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(command, message)

        message = '{}: Aligning to hairpin index'.format(file.basename)
        command = 'bowtie -p 18 -q -l 20 -n 0 -v 2 -a -S --best --strata {} {} --al -S {}'.format(hairpin_index, file.unaligned, file.hairpin_aligned_sam)
        if os.path.exists(file.hairpin_aligned_sam):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(command, message, log_stderr=True)
            self._get_bowtie_summary(self.log_file, 'hairpin')

        message = '{}: Converting SAM to BAM'.format(file.basename)
        command = 'samtools view -S -b {} > {}'.format(file.hairpin_aligned_sam, file.hairpin_aligned_bam)
        if os.path.exists(file.hairpin_aligned_bam):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(command, message)

    def _get_read_counts(self, file):
        """
        Create read count files for mature and hairpin sequences

        Uses samtools and other CLI commands to produce read count file for matures and file for hairpins. File format
        is a text file where the first column is the sequence identifier (e.g. microRNA name) and the second column
        is the number of reads.
        NOTE: only performed if read count files do not already exist

        :param file: File object representing sample to get read counts for
        :type file: File
        """
        message = '{}: Sorting BAM'.format(file.mature_basename)
        command = 'samtools sort -n {} -o {}'.format(file.mature_aligned_bam, file.mature_sorted)
        if os.path.exists(file.mature_sorted):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(command, message)

        message = '{}: Generating read count file'.format(file.mature_basename)
        command = "samtools view {sorted_file_bam} | awk '{{print $3}}' | sort | uniq -c | sort -nr > {readcount_file}".format(
            sorted_file_bam=file.mature_sorted, readcount_file=file.mature_readcount)
        if os.path.exists(file.mature_readcount):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(command, message)

        message = '{}: Sorting BAM'.format(file.hairpin_basename)
        command = 'samtools sort -n {} -o {}'.format(file.hairpin_aligned_bam, file.hairpin_sorted)
        if os.path.exists(file.hairpin_sorted):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(command, message)

        message = '{}: Generating read count file'.format(file.hairpin_basename)
        command = "samtools view {sorted_file_bam} | awk '{{print $3}}' | sort | uniq -c | sort -nr > {readcount_file}".format(
            sorted_file_bam=file.hairpin_sorted, readcount_file=file.hairpin_readcount)
        if os.path.exists(file.hairpin_readcount):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(command, message)

    @staticmethod
    def _run_successful(file):
        """
        Returns true if run was successful for a given sample

        :param file: File object representing sample to check if run was successful
        :type file: File
        """
        # TODO Implement more thoroughly than just checking if file is empty
        return os.stat(file.mature_readcount).st_size >= 0 and os.stat(file.hairpin_readcount).st_size >= 0

    @staticmethod
    def _tail(f, n):
        """
        Helper function to run tail Unix command

        :param f: filepath to run command on
        :type f: str
        :param n: number of lines
        :type n: int
        :return: Returns list of lines of output of tail command
        :rtype: list
        """
        proc = subprocess.Popen(['tail', '-n', str(n), f], stdout=subprocess.PIPE)
        return [line.decode("utf-8") for line in proc.stdout.readlines()]

    def _get_trim_summary(self, log_file):
        """
        Retrieve trimming summary from log file

        :param log_file: filepath to log file
        :type log_file: str
        """
        self.trim_summary = []
        with open(log_file, 'r') as f:
            lines_list = list(f)
            latest_summary_index = max([i for i, x in enumerate(lines_list) if x.startswith('=== Summary ===')])
            for line in lines_list[latest_summary_index+1:]:
                if line.startswith('==='):
                    break
                else:
                    self.trim_summary.append(line)

    def _get_bowtie_summary(self, log_file, bowtie_step):
        """
        Retrieve bowtie summary from log file

        :param log_file: filepath to log file
        :type log_file: str
        :param bowtie_step: 'filtering', 'mature', or 'hairpin'
        :type bowtie_step: str
        """
        if bowtie_step not in ['filtering', 'mature', 'hairpin']:
            raise ValueError('bowtie_step must be "filtering", "mature", or "hairpin"')

        if bowtie_step == 'filtering':
            self.filtering_bowtie_summary = self._tail(log_file, 4)
        elif bowtie_step == 'mature':
            self.mature_bowtie_summary = self._tail(log_file, 4)
        else:
            self.hairpin_bowtie_summary = self._tail(log_file, 4)

    def _write_summary(self, summary_file, file):
        """
        Write file containing summary of trimming, filtering, and aligning steps

        This file is useful for reporting number of reads that pass different stages of the pipeline

        :param summary_file: filepath to summary file
        :type summary_file: str
        :param file: File object to write summary about
        :type file: File
        """
        self._create_summary_file()
        with open(summary_file, 'a') as f:
            f.write('########## {} Processing Summary ##########\n'.format(file.basename))
            f.write('Adapter Trimming Results\n')
            f.write('------------------------\n')
            for line in self.trim_summary:
                f.write(line.replace('\n', '') + '\n')
            f.write('\nNegative Filtering Results\n')
            f.write('--------------------------\n')
            for line in self.filtering_bowtie_summary:
                f.write(line.replace('\n', '') + '\n')
            f.write('\nMature Aligning Results\n')
            f.write('-----------------------\n')
            for line in self.mature_bowtie_summary:
                f.write(line.replace('\n', '') + '\n')
            f.write('\nHairpin Aligning Results\n')
            f.write('------------------------\n')
            for line in self.hairpin_bowtie_summary:
                f.write(line.replace('\n', '') + '\n')
            f.write('\n')

    def run(self):
        """
        Run pipeline.

        This is the main function of the pipeline. After instance of PyPipeline object has been create, this method
        can be called directly. If `self.analysis_only` is True, processing of raw files is skipped and analysis is
        performed.
        """
        # If analysis is only being performed, check if Rscript is installed, run analysis, and return
        if self.analysis_only:
            self._check_program('Rscript')
            self._analyze()
            return

        # Validate all required programs are installed
        for program in ['fastqc', 'fastq-mcf', 'cutadapt', 'bowtie-build', 'bowtie', 'samtools', 'Rscript']:
            self._check_program(program)

        if not self.no_prompts and self.delete is None:
            self.delete = yes_no_dialog(title='Delete files', text='Do you want to delete intermediate files?')
        elif self.delete is None:
            self.delete = True

        if not self.no_prompts and self.no_fastqc is None:
            self.no_fastqc = not yes_no_dialog(title='FastQC', text='Do you want to perform FastQC on all files?')
        elif self.no_fastqc is None:
            self.no_fastqc = False

        # Build all three bowtie indexes
        self._build_index(self.negative_references, self.negative_index_dir, 'negative')
        self._build_index(self.mature_references, self.mature_index_dir, 'mature')
        self._build_index(self.hairpin_references, self.hairpin_index_dir, 'hairpin')

        if not self.no_fastqc:
            for file in self.files:
                self._fastqc_check(file)

        # Main file processing steps
        for file in self.files:
            if not file.read_counts_exist():
                self._trim_adapters(file)
                self._filter_out_neg(file)
                self._align_reads(file)
                self._get_read_counts(file)
                self._write_summary(self.summary_file, file)

                if self.delete and self._run_successful(file):
                    self._log_message('{}: Deleting intermediate files'.format(file.basename))
                    file.remove_intermediates()
                elif not self._run_successful(file):
                    self._log_message('{}: Run was not successful'.format(file.basename), command_status=self.EXITING)
                    exit(1)
                else:
                    pass
            else:
                self._log_message('{}: Read counts already created'.format(file.basename), command_status=self.FILE_ALREADY_EXISTS)

        if not self.no_analysis:
            self._analyze()

    def _create_conditions_file(self):
        """
        Create conditions file needed for analysis by RBioMir R packages

        File is csv file where first column is sample name and the second column is sample
        condition ('control' or 'stress')
        """
        tmp_list = sorted([[sample, condition] for (sample, condition) in self.sample_conditions.items()])
        csv_data = [['sample', 'condition']]
        csv_data.extend(tmp_list)

        with open(self.conditions_file, 'w') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerows(csv_data)

    def _copy_read_counts(self):
        """
        Copy read counts to what will be working directory for R script analysis
        """
        for file in self.files:
            shutil.copy(file.mature_readcount, self.figures)

    def _run_rscript(self):
        """
        Run R script

        This pipeline was designed to be analysed using an R script composed from RBioMir package functions. The path to
        the R script is specified in the config file and must take the following command line arguments in order:

        * `Path to working directory`
        * `Path to KEGG GMT file`
        * `Path to GO BP GMT file`
        * `Path to GO MF GMT file`
        * `Path to GO CC GMT file`
        """
        message = 'Running R analysis'
        command = 'Rscript {rpipeline} {wd} {kegg_ids} {go_bp_ids} {go_mf_ids} {go_cc_ids}'.format(rpipeline=self.rpipeline,
                                                                                                   wd=self.figures,
                                                                                                   kegg_ids=self.kegg_id_file,
                                                                                                   go_bp_ids=self.go_bp_id_file,
                                                                                                   go_mf_ids=self.go_mf_id_file,
                                                                                                   go_cc_ids=self.go_cc_id_file)
        self._run_command(command, message, log_stdout=True)

    @staticmethod
    def _move_files_by_regex(source, dest=None, pattern=None):
        """
        Move files matching regex pattern from source to destination

        :param source: source directory
        :type source: str
        :param dest: destination directory
        :param dest: str
        :param pattern: regex pattern
        :type pattern: str
        """
        for f in os.listdir(source):
            if re.search(pattern, f):
                if dest is None:
                    os.remove(os.path.join(source, f))
                else:
                    os.rename(os.path.join(source, f), os.path.join(dest, f))

    def _clean_up(self):
        """
        Clean up directory after R analysis completes

        This method deletes the conditions file created before R script is run, new copies of the read_count files,
        and moves the microRNA target files to the microRNA targets directory
        """
        self._log_message('Cleaning up directories')
        os.remove(self.conditions_file)
        self._move_files_by_regex(source=self.figures, dest=self.mirna_targets_dir, pattern=r'hsa.*\.csv')
        self._move_files_by_regex(source=self.figures, dest=None, pattern=r'.*read_count.txt')

    def _analyze(self):
        """
        Run analysis on read count files

        This method sets up the nessessary conditions for the R script to be run, runs the R analysis, then cleans up
        """
        os.makedirs(self.figures, exist_ok=True)
        os.makedirs(self.mirna_targets_dir, exist_ok=True)

        # Validate the sample conditions specified in config match specific files
        self._validate_sample_conditions()
        self._copy_read_counts()
        # "conditions.csv" files must be made for rbiomir
        self._create_conditions_file()
        self._run_rscript()
        self._clean_up()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-c', '--config', action='store', metavar='<config_file>', help='Path to config file')
    group.add_argument('-d', '--config-dir', action='store', metavar='<config_dir>', help='Directory containing config files')
    parser.add_argument('--no-prompts', action='store_true', default=False, help='Suppress user prompts')
    parser.add_argument('--no-fastqc', action='store_true', default=False, help='Do not perform fastqc on raw files')
    parser.add_argument('--delete', action='store_true', default=None, help='Delete intermediate processing files')
    analysis_group = parser.add_mutually_exclusive_group()
    analysis_group.add_argument('--no-analysis', action='store_true', default=False, help='Do not perform R analysis')
    analysis_group.add_argument('--analysis-only', action='store', dest='read_count_dir', metavar='<read_count_dir>', default=None, type=os.path.abspath, help='Run analysis only on read counts in supplied directory')
    args = parser.parse_args()

    # Path to config (or dir to multiple configs) can be passed as command line arguments
    if args.config:
        configs = [args.config]
    elif args.config_dir:
        configs = []
        for dirpath, _, filenames in os.walk(args.config_dir):
            for f in sorted([f for f in filenames if f.endswith('.ini')]):
                configs.append(os.path.abspath(os.path.join(dirpath, f)))
    else:
        configs = ['example_config.ini']

    # Set up pipelines
    pipelines = []
    for config in configs:
        pipelines.append(PyPipeline(config,
                                    no_prompts=args.no_prompts,
                                    no_fastqc=args.no_fastqc,
                                    delete=args.delete,
                                    no_analysis=args.no_analysis,
                                    read_count_dir=args.read_count_dir))

    for pipeline in pipelines:
        pipeline.run()
