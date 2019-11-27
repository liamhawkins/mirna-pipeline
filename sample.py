import os


class Sample:
    def __init__(self, raw_path, analysis_dir):
        """
        The `Sample` class creates and manages filepaths of a raw reads file and all
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


