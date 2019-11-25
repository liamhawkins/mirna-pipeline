#!/usr/bin/python3

# TODO: Validate all FASTQ files
import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from configparser import ConfigParser
from datetime import datetime

from prompt_toolkit import HTML, print_formatted_text
from prompt_toolkit.shortcuts import yes_no_dialog


class File:
    def __init__(self, raw_path, analysis_dir):
        self.raw = raw_path
        self.analysis_dir = analysis_dir
        self.basename = self._get_basename(raw_path)
        self.trimmed = self._create_file('.trimmed.fastq')
        self.temp = self._create_file('.tmp.fastq')
        self.filtered = self._create_file('.filtered.fastq')
        self.mature_aligned_sam = self._create_file('_MATURE.aligned.sam')
        self.mature_aligned_bam = self._create_file('_MATURE.aligned.bam')
        self.unaligned = self._create_file('.unaligned.fastq')
        self.hairpin_aligned_sam = self._create_file('_HAIRPIN.aligned.sam')
        self.hairpin_aligned_bam = self._create_file('_HAIRPIN.aligned.bam')
        self.mature_basename = self._get_basename(self.mature_aligned_sam)
        self.hairpin_basename = self._get_basename(self.hairpin_aligned_sam)
        self.mature_sorted = self._create_file('.sorted.bam', file=self.mature_aligned_bam)
        self.hairpin_sorted = self._create_file('.sorted.bam', file=self.hairpin_aligned_bam)
        readcount_dir = os.path.join(self.analysis_dir, 'read_counts/')
        os.makedirs(readcount_dir, exist_ok=True)
        self.mature_readcount = self._create_file('.read_count.txt', file=self.mature_sorted, dir=readcount_dir)
        self.hairpin_readcount = self._create_file('.read_count.txt', file=self.hairpin_sorted, dir=readcount_dir)

    def __repr__(self):
        return '{}({}, {})'.format(self.__class__.__name__, self.raw, self.analysis_dir)

    @staticmethod
    def _get_basename(path):
        return os.path.splitext(os.path.basename(path))[0].split('.')[0]

    def _create_file(self, postfix, file=None, dir=None):
        if file is None:
            basename = self.basename
        else:
            basename = self._get_basename(file)

        if dir is None:
            dir = self.analysis_dir

        return os.path.join(dir, basename + postfix)

    def remove_intermediates(self):
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

    def read_counts_exist(self):
        return os.path.isfile(self.mature_readcount) and os.path.isfile(self.hairpin_readcount)


class PyPipeline:
    def __init__(self, config_file, no_prompts=False, no_fastqc=None, delete=None, no_analysis=False):
        self.pipeline = lambda: datetime.now().strftime("%Y-%m-%d %H:%M")
        self.no_prompts = no_prompts
        self.no_fastqc = no_fastqc
        self.delete = delete
        self.no_analysis = no_analysis
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

        # Set up directories
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
        self.F_PIPELINE = lambda: '<teal>{}</teal>'.format(self.pipeline())

        # Create log file
        os.makedirs(self.analysis_dir, exist_ok=True)
        self._create_log_file()
        self._create_summary_file()

        self.files = []
        for dirpath, _, filenames in os.walk(self.raw_files_dir):
            for f in sorted([f for f in filenames if f.endswith(('.fastq', '.fq'))]):
                abs_path = os.path.abspath(os.path.join(dirpath, f))
                self.files.append(File(abs_path, self.analysis_dir))

        # Set up config-dependent variables
        self.adapters = None
        self.trim_6 = None
        self._validate_config()

    def _run_command(self, message, command, log_stdout=False, log_stderr=False):
        if self.command_log:
            with open(self.command_log, 'a') as f:
                f.write(command + '\n')

        formatted_message = '[{}] '.format(self.F_PIPELINE()) + message + '... '
        unformated_message = '[{}] '.format(self.pipeline()) + message + '... '
        print_formatted_text(HTML(formatted_message), end='', flush=True)
        with open(self.log_file, 'a') as f:
            f.write(unformated_message + '\n')

        try:
            if log_stdout and log_stderr:
                subprocess.call(command + ' &>> {}'.format(self.log_file), shell=True, stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL)
            elif log_stdout:
                subprocess.call(command + ' >> {}'.format(self.log_file), shell=True, stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL)
            elif log_stderr:
                subprocess.call(command + ' 2>> {}'.format(self.log_file), shell=True, stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL)
            else:
                subprocess.call(command, shell=True, stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError as exc:
            print_formatted_text(self.BAD)
            print('ERROR:')
            print(exc.output.decode('utf-8'))
            with open(self.log_file, 'a') as f:
                f.write(exc.output.decode('utf-8'))
            exit(1)

        print_formatted_text(self.GOOD)

    def _log_message(self, message, command_status=None, **kwargs):
        if command_status is None:
            command_status = self.GOOD

        formatted_message = '[{}] '.format(self.F_PIPELINE()) + message + '... '
        unformated_message = '[{}] '.format(self.pipeline()) + message + '... '
        print_formatted_text(HTML(formatted_message + command_status.value), **kwargs)

        with open(self.log_file, 'a') as f:
            f.write(unformated_message + '\n')

    def _create_log_file(self):
        message = 'Creating log file {}'.format(os.path.basename(self.log_file))
        command = 'touch {}'.format(self.log_file)
        self._run_command(message, command)

    def _create_summary_file(self):
        message = 'Creating summary file - {}'.format(os.path.basename(self.summary_file))
        command = 'touch {}'.format(self.summary_file)
        self._run_command(message, command)

    def _validate_file(self, file_):
        if not os.path.isfile(file_):
            self._log_message('{} does not exist'.format(file_), command_status=self.EXITING)
            exit(1)

    def _validate_sample_conditions(self):
        for file in self.files:
            if file.basename not in self.sample_conditions.keys():
                self._log_message('Cannot find sample condition in config file: {}'.format(file.basename), command_status=self.EXITING)
                exit(1)

        if any([condition not in ['control', 'stress'] for condition in self.sample_conditions.values()]):
            self._log_message('Sample conditions can only be "control" or "stress" at this time', command_status=self.EXITING)
            exit(1)

    def _validate_config(self):
        self._log_message('Performing config validation', command_status=self.NONE, end='', flush=True)

        if self.company.upper() == 'TORONTO':
            self.adapters = self.toronto_adapters
            self.trim_6 = False
        elif self.company.upper() == 'BC':
            self.adapters = self.bc_adapters
            self.trim_6 = True
        else:
            self._log_message('COMPANY must be "BC" or "TORONTO"', command_status=self.EXITING)
            exit(1)

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

        print_formatted_text(self.GOOD)

        files = '\n'.join([file for file in sorted(os.listdir(self.raw_files_dir)) if file.endswith('.fastq') or file.endswith('.fq')])
        if not self.no_prompts:
            continue_ = yes_no_dialog(title='File check', text='Are these the files you want to process?\n\n' + files)
            if not continue_:
                exit(0)

        self._validate_sample_conditions()

    def _check_program(self, program):
        self._log_message('Checking that {} is installed'.format(program), command_status=self.NONE, end='', flush=True)

        try:
            subprocess.Popen([program], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).communicate()
            print_formatted_text(self.GOOD)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                self._log_message('The program {} was not found'.format(program), command_status=self.BAD)
                exit()
            else:
                self._log_message('An unknown error occurred when looking for {}'.format(program), command_status=self.BAD)
                raise

    def _index_is_built(self, ind_prefix, index_name):
        try:
            os.listdir(os.path.dirname(ind_prefix))
        except FileNotFoundError:
            self._log_message('Checking if {} index is built'.format(index_name), command_status=self.NOT_BUILT)
            return False

        for filename in os.listdir(ind_prefix):
            if not filename.startswith(os.path.basename(ind_prefix)) or not filename.endswith('.ebwt'):
                self._log_message('Checking if {} index is built'.format(index_name), command_status=self.NOT_BUILT)
                return False

        self._log_message('Checking if {} index is built'.format(index_name))
        return True

    def _build_index(self, index_dir, index_name):
        if not self._index_is_built(index_dir, index_name):
            os.makedirs(index_dir, exist_ok=True)

            message = 'Building negative index'
            command = 'bowtie-build {} {}'.format(self.negative_references, os.path.join(index_dir, os.path.basename(index_dir)))
            self._run_command(message, command)

    def _fastqc_check(self, file):
        os.makedirs(self.fastqc_dir, exist_ok=True)
        message = 'Performing fastqc check on {}'.format(file.basename)
        command = 'fastqc -q {} -o {}'.format(file.raw, self.fastqc_dir)
        self._run_command(message, command)

    def _trim_adapters(self, file):
        message = '{}: Trimming adapters'.format(file.basename)
        command = 'cutadapt -q 20 -m 10 -j 18 -b file:{0} {1} -o {2}'.format(self.adapters, file.raw, file.temp)
        if os.path.exists(file.trimmed):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(message, command, log_stdout=True)
            self.get_trim_summary(self.log_file)

            if self.trim_6:
                message = '{}: Trimming 6 nucleotides'.format(file.basename)
                command = 'cutadapt -u 6 -j 18 {0} -o {1}'.format(file.temp, file.trimmed)
                self._run_command(message, command, log_stdout=True)
                os.remove(file.temp)
            else:
                os.rename(file.temp, file.trimmed)

    def _filter_out_neg(self, file):
        negative_index = os.path.join(self.negative_index_dir, os.path.basename(self.negative_index_dir))

        message = '{}: Filtering negative RNA species'.format(file.basename)
        command = 'bowtie -p 18 -q {} {} --un {}'.format(negative_index, file.trimmed, file.filtered)
        if os.path.exists(file.filtered):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(message, command, log_stderr=True)
            self.get_bowtie_summary(self.log_file, 'filtering')

    def _align_reads(self, file):
        mature_index = os.path.join(self.mature_index_dir, os.path.basename(self.mature_index_dir))
        hairpin_index = os.path.join(self.hairpin_index_dir, os.path.basename(self.hairpin_index_dir))

        message = '{}: Aligning to mature index'.format(file.basename)
        command = 'bowtie -p 18 -q -l 20 -n 0 -v 2 -a -S --best --strata {} {} --al -S {} --un {}'.format(
            mature_index, file.filtered, file.mature_aligned_sam, file.unaligned)
        if os.path.exists(file.mature_aligned_sam) and os.path.exists(file.unaligned):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(message, command, log_stderr=True)
            self.get_bowtie_summary(self.log_file, 'mature')

        message = '{}: Converting SAM to BAM'.format(file.basename)
        command = 'samtools view -S -b {} > {}'.format(file.mature_aligned_sam, file.mature_aligned_bam)
        if os.path.exists(file.mature_aligned_bam):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(message, command)

        message = '{}: Aligning to hairpin index'.format(file.basename)
        command = 'bowtie -p 18 -q -l 20 -n 0 -v 2 -a -S --best --strata {} {} --al -S {}'.format(hairpin_index, file.unaligned, file.hairpin_aligned_sam)
        if os.path.exists(file.hairpin_aligned_sam):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(message, command, log_stderr=True)
            self.get_bowtie_summary(self.log_file, 'hairpin')

        message = '{}: Converting SAM to BAM'.format(file.basename)
        command = 'samtools view -S -b {} > {}'.format(file.hairpin_aligned_sam, file.hairpin_aligned_bam)
        if os.path.exists(file.hairpin_aligned_bam):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(message, command)

    def _get_read_counts(self, file):
        message = '{}: Sorting BAM'.format(file.mature_basename)
        command = 'samtools sort -n {} -o {}'.format(file.mature_aligned_bam, file.mature_sorted)
        if os.path.exists(file.mature_sorted):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(message, command)

        message = '{}: Generating read count file'.format(file.mature_basename)
        command = "samtools view {sorted_file_bam} | awk '{{print $3}}' | sort | uniq -c | sort -nr > {readcount_file}".format(
            sorted_file_bam=file.mature_sorted, readcount_file=file.mature_readcount)
        if os.path.exists(file.mature_readcount):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(message, command)

        message = '{}: Sorting BAM'.format(file.hairpin_basename)
        command = 'samtools sort -n {} -o {}'.format(file.hairpin_aligned_bam, file.hairpin_sorted)
        if os.path.exists(file.hairpin_sorted):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(message, command)

        message = '{}: Generating read count file'.format(file.hairpin_basename)
        command = "samtools view {sorted_file_bam} | awk '{{print $3}}' | sort | uniq -c | sort -nr > {readcount_file}".format(
            sorted_file_bam=file.hairpin_sorted, readcount_file=file.hairpin_readcount)
        if os.path.exists(file.hairpin_readcount):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(message, command)

    @staticmethod
    def _run_successful(file):
        # TODO Implement more thoroughly than just checking if file is empty
        return os.stat(file.mature_readcount).st_size >= 0 and os.stat(file.hairpin_readcount).st_size >= 0

    @staticmethod
    def tail(f, n):
        proc = subprocess.Popen(['tail', '-n', str(n), f], stdout=subprocess.PIPE)
        return [line.decode("utf-8") for line in proc.stdout.readlines()]

    def get_trim_summary(self, log_file):
        self.trim_summary = []
        with open(log_file, 'r') as f:
            lines_list = list(f)
            latest_summary_index = max([i for i, x in enumerate(lines_list) if x.startswith('=== Summary ===')])
            for line in lines_list[latest_summary_index+1:]:
                if line.startswith('==='):
                    break
                else:
                    self.trim_summary.append(line)

    def get_bowtie_summary(self, log_file, bowtie_step):
        if bowtie_step not in ['filtering', 'mature', 'hairpin']:
            raise ValueError('bowtie_step must be "filtering", "mature", or "hairpin"')

        if bowtie_step == 'filtering':
            self.filtering_bowtie_summary = self.tail(log_file, 4)
        elif bowtie_step == 'mature':
            self.mature_bowtie_summary = self.tail(log_file, 4)
        else:
            self.hairpin_bowtie_summary = self.tail(log_file, 4)

    def write_summary(self, summary_file, file_basename):
        with open(summary_file, 'a') as f:
            f.write('########## {} Processing Summary ##########\n'.format(file_basename))
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
        # Validate all required programs are installed
        for program in ['fastqc', 'fastq-mcf', 'cutadapt', 'bowtie-build', 'bowtie', 'samtools', 'Rscript']:
            self._check_program(program)

        self._build_index(self.negative_index_dir, 'negative')
        self._build_index(self.mature_index_dir, 'mature')
        self._build_index(self.hairpin_index_dir, 'hairpin')

        if not self.no_prompts and self.delete is None:
            self.delete = yes_no_dialog(title='Delete files', text='Do you want to delete intermediate files?')
        elif self.delete is None:
            self.delete = True

        if not self.no_prompts and self.no_fastqc is None:
            self.no_fastqc = not yes_no_dialog(title='FastQC', text='Do you want to perform FastQC on all files?')
        elif self.no_fastqc is None:
            self.no_fastqc = False

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
                self.write_summary(self.summary_file, file.basename)

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
            self.analyze()

    def _create_conditions_file(self):
        tmp_list = sorted([[sample, condition] for (sample, condition) in self.sample_conditions.items()])
        csv_data = [['sample', 'condition']]
        csv_data.extend(tmp_list)

        with open(self.conditions_file, 'w') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerows(csv_data)

    def _copy_read_counts(self):
        for file in self.files:
            shutil.copy(file.mature_readcount, self.figures)

    def _run_rscript(self):
        message = 'Running R analyis'
        command = 'Rscript {rpipeline} {wd} {kegg_ids} {go_bp_ids} {go_mf_ids} {go_cc_ids}'.format(rpipeline=self.rpipeline,
                                                                                                   wd=self.figures,
                                                                                                   kegg_ids=self.kegg_id_file,
                                                                                                   go_bp_ids=self.go_bp_id_file,
                                                                                                   go_mf_ids=self.go_mf_id_file,
                                                                                                   go_cc_ids=self.go_cc_id_file)
        self._run_command(message, command, log_stdout=True)

    @staticmethod
    def _move_files_by_regex(source, dest=None, pattern=None):
        for f in os.listdir(source):
            if re.search(pattern, f):
                if dest is None:
                    os.remove(os.path.join(source, f))
                else:
                    os.rename(os.path.join(source, f), os.path.join(dest, f))

    def _clean_up(self):
        self._log_message('Cleaning up directories')
        os.remove(self.conditions_file)
        self._move_files_by_regex(source=self.figures, dest=self.mirna_targets_dir, pattern=r'hsa.*\.csv')
        self._move_files_by_regex(source=self.figures, dest=None, pattern=r'.*read_count.txt')

    def analyze(self):
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
    group.add_argument('-c', '--config', action='store', help='Path to config file')
    group.add_argument('-d', '--config-dir', action='store', help='Directory containing config files')
    parser.add_argument('--no-prompts', action='store_true', default=False, help='Suppress user prompts')
    parser.add_argument('--no-fastqc', action='store_true', default=False, help='Do not perform fastqc on raw files')
    parser.add_argument('--delete', action='store_true', default=None, help='Delete intermediate processing files')
    parser.add_argument('--no-analysis', action='store_true', default=False, help='Do not perform R analysis')
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
        pipelines.append(PyPipeline(config, no_prompts=args.no_prompts, no_fastqc=args.no_fastqc, delete=args.delete, no_analysis=args.no_analysis))

    for pipeline in pipelines:
        pipeline.run()
