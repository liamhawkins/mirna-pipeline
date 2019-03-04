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
            try:
                os.remove(file_)
            except FileNotFoundError:
                pass

    def read_counts_exist(self):
        return os.path.isfile(self.mature_readcount) and os.path.isfile(self.hairpin_readcount)


class PyPipeline:
    def __init__(self, config_file, no_prompts=False, fastqc=None, delete=None):
        self.pipeline = 'MicroRNA-seq PyPipeline'
        self.no_prompts = no_prompts
        self.fastqc = fastqc
        self.delete = delete

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
        self.bowtie_dir = config['BOWTIE_DIR']
        self.condition_file = config['CONDITION_FILE']
        self.control_name = config['CONTROL_NAME']

        # Set up directories
        self.log_file = os.path.join(self.analysis_dir, datetime.now().strftime('%d-%m-%y_%H:%M') + '.log')
        self.fastqc_dir = os.path.join(self.analysis_dir, 'fastqc/')
        self.negative_index_dir = os.path.join(self.bowtie_dir, 'neg_ref')
        self.hairpin_index_dir = os.path.join(self.bowtie_dir, 'hp_ref')
        self.mature_index_dir = os.path.join(self.bowtie_dir, 'mature_ref')

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

        self.files = []
        for dirpath, _, filenames in os.walk(self.raw_files_dir):
            for f in sorted(filenames):
                abs_path = os.path.abspath(os.path.join(dirpath, f))
                self.files.append(File(abs_path, self.analysis_dir))

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
        if not self.no_prompts:
            continue_ = yes_no_dialog(title='File check', text='Are these the files you want to process?\n\n' + files)
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
            self._run_command(message, command, log_output=True)

            if self.trim_6:
                message = '{}: Trimming 6 nucleotides'.format(file.basename)
                command = 'cutadapt -u 6 -j 18 {0} -o {1}'.format(file.temp, file.trimmed)
                self._run_command(message, command, log_output=True)
                os.remove(file.temp)
            else:
                os.rename(file.temp, file.trimmed)

    def _filter_out_neg(self, file):
        negative_index = os.path.join(self.negative_index_dir, os.path.basename(self.negative_index_dir))

        message = '{}: Filtering negative RNA species'.format(file.basename)
        command = 'bowtie -p 18 --quiet -q {} {} --un {}'.format(negative_index, file.trimmed, file.filtered)
        if os.path.exists(file.filtered):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(message, command, log_output=True)

    def _align_reads(self, file):
        mature_index = os.path.join(self.mature_index_dir, os.path.basename(self.mature_index_dir))
        hairpin_index = os.path.join(self.hairpin_index_dir, os.path.basename(self.hairpin_index_dir))

        message = '{}: Aligning to mature index'.format(file.basename)
        command = 'bowtie -p 18 --quiet -q -l 20 -n 0 -v 2 -a -S --best --strata {} {} --al -S {} --un {}'.format(
            mature_index, file.filtered, file.mature_aligned_sam, file.unaligned)
        if os.path.exists(file.mature_aligned_sam) and os.path.exists(file.unaligned):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(message, command, log_output=True)

        message = '{}: Converting SAM to BAM'.format(file.basename)
        command = 'samtools view -S -b {} > {}'.format(file.mature_aligned_sam, file.mature_aligned_bam)
        if os.path.exists(file.mature_aligned_bam):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(message, command)

        message = '{}: Aligning to hairpin index'.format(file.basename)
        command = 'bowtie -p 18 --quiet -q -l 20 -n 0 -v 2 -a -S --best --strata {} {} --al -S {}'.format(hairpin_index, file.unaligned, file.hairpin_aligned_sam)
        if os.path.exists(file.hairpin_aligned_sam):
            self._log_message(message, command_status=self.FILE_ALREADY_EXISTS)
        else:
            self._run_command(message, command, log_output=True)

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

    def _run_successful(self):
        # TODO
        return True

    def run(self):
        # Validate all required programs are installed
        for program in ['fastqc', 'fastq-mcf', 'cutadapt', 'bowtie-build', 'bowtie', 'samtools', 'Rscript']:
            self._check_program(program)

        self._build_index(self.negative_index_dir, 'negative')
        self._build_index(self.mature_index_dir, 'mature')
        self._build_index(self.hairpin_index_dir, 'hairpin')

        if not self.no_prompts:
            self.delete = yes_no_dialog(title='Delete files', text='Do you want to delete intermediate files?')
        elif self.delete is None:
            self.delete = True

        if not self.no_prompts:
            self.fastqc = yes_no_dialog(title='FastQC', text='Do you want to perform FastQC on all files?')
        elif self.fastqc is None:
            self.fastqc = True

        if self.fastqc:
            for file in self.files:
                self._fastqc_check(file)

        # Main file processing steps
        for file in self.files:
            if not file.read_counts_exist():
                self._trim_adapters(file)
                self._filter_out_neg(file)
                self._align_reads(file)
                self._get_read_counts(file)

            if self.delete and self._run_successful(file):
                self._log_message('{}: Deleting intermediate files'.format(file.basename))
                file.remove_intermediates()
                self._log_message('Clearing log file')
                open(self.log_file, 'w').close()
            elif not self._run_successful(file):
                self._log_message('{}: Run was not successful'.format(file.basename), command_status=self.EXITING)
                exit(1)
            else:
                pass

        self._log_message('Deleting log file')
        os.remove(self.log_file)


if __name__ == '__main__':
    pipeline = PyPipeline('config.ini', no_prompts=True, fastqc=False)
    pipeline.run()
