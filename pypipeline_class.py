import os
import subprocess
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
        self.filtered = self._create_file('.filtered.fastq')
        self.mature_aligned_sam = self._create_file('._MATURE.aligned.sam')
        self.mature_aligned_bam = self._create_file('._MATURE.aligned.bam')
        self.unaligned = self._create_file('.unaligned.fastq')
        self.hairpin_aligned_sam = self._create_file('._HAIRPIN.aligned.sam')
        self.hairpin_aligned_bam = self._create_file('._HAIRPIN.aligned.bam')
        self.mature_sorted = self._create_file('.sorted.bam', file=self.mature_aligned_bam)
        self.hairpin_sorted = self._create_file('.sorted.bam', file=self.hairpin_aligned_bam)
        self.mature_readcount = self._create_file('.read_count.txt', file=self.mature_sorted)
        self.hairpin_readcount = self._create_file('.read_count.txt', file=self.mature_sorted)

    def __repr__(self):
        return '{}({}, {})'.format(self.__class__.__name__, self.raw, self.analysis_dir)

    @staticmethod
    def _get_basename(path):
        return os.path.splitext(os.path.basename(path))[0].split('.')[0]

    def _create_file(self, postfix, file=None):
        if file is None:
            basename = self.basename
        else:
            basename = self._get_basename(file)

        return os.path.join(self.analysis_dir, basename + postfix)

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
            os.remove(file_)


class PyPipeline:
    def __init__(self, config_file, debug=False):
        self.debug = debug
        self.delete = True

        # Read in Config File
        config = ConfigParser()
        config.read(config_file)
        config = config[config.sections()[0]]
        self.company = config['COMPANY']
        self.analysis_dir = config['ANALYSIS_DIR']
        self.raw_files_dir = config['RAW_FILES_DIR']
        self.toronto_adapters = config['TORONTO_ADAPTERS']
        self.bc_adapters = config['BC_ADAPTERS']
        self.negative_references = config['NEGATIVE_REFERENCE_FILE']
        self.mature_references = config['MATURE_REFERENCE_FILE']
        self.hairpin_references = config['HAIRPIN_REFERENCE_FILE']
        self.condition_file = config['CONDITION_FILE']
        self.control_name = config['CONTROL_NAME']

        # Set up directories
        self.pipeline = 'MicroRNA-seq PyPipeline'
        self.log_file = os.path.join(self.analysis_dir, datetime.now().strftime('%d-%m-%y_%H:%M') + '.log')
        self.fastqc_dir = os.path.join(self.analysis_dir, 'fastqc/')
        self.trimmed_dir = os.path.join(self.analysis_dir, 'trimmed/')
        self.adapter_dir = os.path.join(self.analysis_dir, 'adapters/')
        self.bowtie_dir = os.path.join('/home/student/Desktop/pypipeline', 'bowtie_index/')
        self.negative_index_dir = os.path.join(self.bowtie_dir, 'neg_ref')
        self.hairpin_index_dir = os.path.join(self.bowtie_dir, 'hp_ref')
        self.mature_index_dir = os.path.join(self.bowtie_dir, 'mature_ref')
        self.filtered_dir = os.path.join(self.analysis_dir, 'filtered/')
        self.mature_dir = os.path.join(self.analysis_dir, 'matures/')
        self.mature_aligned_dir = os.path.join(self.mature_dir, 'aligned/')
        self.mature_unaligned_dir = os.path.join(self.mature_dir, 'unaligned/')
        self.hairpin_dir = os.path.join(self.analysis_dir, 'hairpins/')
        self.hairpin_aligned_dir = os.path.join(self.hairpin_dir, 'aligned/')
        self.read_count_dir = os.path.join(self.analysis_dir, 'read_counts/')
        self.mature_read_count_dir = os.path.join(self.read_count_dir, 'mature/')
        self.hairpin_read_count_dir = os.path.join(self.read_count_dir, 'hairpin/')

        # Formatted strings
        self.GOOD = HTML('<green>GOOD</green>')
        self.FILE_ALREADY_EXISTS = HTML('<yellow>FILE ALREADY EXISTS</yellow>')
        self.NOT_BUILT = HTML('<yellow>NOT BUILT</yellow>')
        self.BAD = HTML('<red>BAD</red>')
        self.EXITING = HTML('<red>EXITING</red>')
        self.NONE = HTML('')
        self.F_PIPELINE = HTML('<teal>{}</teal>'.format(self.pipeline))

        # Create log file
        os.makedirs(self.analysis_dir, exist_ok=True)
        self._create_log_file()

        # Set up config-dependent variables
        self.adapters = None
        self.trim_6 = None
        self._validate_config()

    def _run_command(self, message, command, log_output=False):
        formatted_message = '[{}] '.format(self.pipeline) + message + '... '
        print(formatted_message, end='', flush=True)
        with open(self.log_file, 'a') as f:
            f.write(formatted_message)

        try:
            if log_output:
                subprocess.call(command + ' >> {}'.format(self.log_file), shell=True, stderr=subprocess.STDOUT,
                                stdout=subprocess.DEVNULL)
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

        formatted_message = '[{}] '.format(self.pipeline) + message + '... '
        print_formatted_text(HTML(formatted_message + command_status.value), **kwargs)

        with open(self.log_file, 'a') as f:
            f.write(formatted_message + '\n')

    def _create_log_file(self):
        message = 'Creating log file {}'.format(os.path.basename(self.log_file))
        command = 'touch {}'.format(self.log_file)
        self._run_command(message, command)

    def _validate_file(self, file_):
        if not os.path.isfile(file_):
            self._log_message('{} does not exist'.format(file_), command_status=self.EXITING)
            exit(1)

    def _validate_config(self):
        self._log_message('Performing config validation', command_status=self.NONE, end='', flush=True)
        self._validate_file(self.negative_references)
        self._validate_file(self.mature_references)
        self._validate_file(self.hairpin_references)
        self._validate_file(self.condition_file)

        if self.company.upper() == 'TORONTO':
            self.adapters = self.toronto_adapters
            self.trim_6 = False
        elif self.company.upper() == 'BC':
            self.adapters = self.bc_adapters
            self.trim_6 = True
        else:
            self._log_message('COMPANY must be "BC" or "TORONTO"', command_status=self.EXITING)
            exit(1)

        if not os.path.isfile(self.adapters):
            self._log_message('{} is missing, exiting'.format(self.adapters), command_status=self.EXITING)
            exit(1)

        print_formatted_text(self.GOOD)

        files = '\n'.join([file for file in sorted(os.listdir(self.raw_files_dir)) if file.endswith('.fastq') or file.endswith('.fq')])
        if not self.debug:
            continue_ = yes_no_dialog(title='File check', text='Are these the files you want to process?\n\n' + files)
        else:
            print(files)
            prompt = input('Are these the files you want to process?')
            if prompt.lower() in ['y', 'yes']:
                continue_ = True
            else:
                continue_ = False

        if not continue_:
            exit(0)

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

    def run(self):
        # Validate all required programs are installed
        for program in ['fastqc', 'fastq-mcf', 'cutadapt', 'bowtie-build', 'bowtie', 'samtools', 'Rscript']:
            self._check_program(program)


if __name__ == '__main__':
    pipeline = PyPipeline('config.ini', debug=True)
    pipeline.run()
