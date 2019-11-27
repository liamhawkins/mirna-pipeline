import os
import shutil
import unittest

from sample import Sample


class TestSampleInit(unittest.TestCase):
    def setUp(self):
        self.raw_path = '/some/path/to/a/file1.fastq'.replace('/', os.sep)
        self.analysis_dir = '/path/to/analysis_dir'.replace('/', os.sep)
        self.sample = Sample(raw_path=self.raw_path, analysis_dir=self.analysis_dir)

    def test_repr(self):
        self.assertEqual(repr(self.sample), 'Sample(/some/path/to/a/file1.fastq, /path/to/analysis_dir)'.replace('/', os.sep))

    def test_sample_basename(self):
        self.assertEqual(self.sample.basename, 'file1')

    def test_mature_basename(self):
        self.assertEqual(self.sample.mature_basename, 'file1_MATURE')

    def test_hairpin_basename(self):
        self.assertEqual(self.sample.hairpin_basename, 'file1_HAIRPIN')

    def test_sample_trimmed_name(self):
        self.assertEqual(self.sample.trimmed, '/path/to/analysis_dir/file1.trimmed.fastq'.replace('/', os.sep))

    def test_sample_tmp_name(self):
        self.assertEqual(self.sample.temp, '/path/to/analysis_dir/file1.tmp.fastq'.replace('/', os.sep))

    def test_sample_filtered_name(self):
        self.assertEqual(self.sample.filtered, '/path/to/analysis_dir/file1.filtered.fastq'.replace('/', os.sep))

    def test_sample_mature_aligned_sam_name(self):
        self.assertEqual(self.sample.mature_aligned_sam, '/path/to/analysis_dir/file1_MATURE.aligned.sam'.replace('/', os.sep))

    def test_sample_mature_aligned_bam_name(self):
        self.assertEqual(self.sample.mature_aligned_bam, '/path/to/analysis_dir/file1_MATURE.aligned.bam'.replace('/', os.sep))

    def test_sample_hairpin_aligned_sam_name(self):
        self.assertEqual(self.sample.hairpin_aligned_sam, '/path/to/analysis_dir/file1_HAIRPIN.aligned.sam'.replace('/', os.sep))

    def test_sample_hairpin_aligned_bam_name(self):
        self.assertEqual(self.sample.hairpin_aligned_bam, '/path/to/analysis_dir/file1_HAIRPIN.aligned.bam'.replace('/', os.sep))

    def test_sample_mature_sorted_name(self):
        self.assertEqual(self.sample.mature_sorted, '/path/to/analysis_dir/file1_MATURE.sorted.bam'.replace('/', os.sep))

    def test_sample_hairpin_sorted_name(self):
        self.assertEqual(self.sample.hairpin_sorted, '/path/to/analysis_dir/file1_HAIRPIN.sorted.bam'.replace('/', os.sep))

    def test_sample_mature_readcount_name(self):
        self.assertEqual(self.sample.mature_readcount, '/path/to/analysis_dir/read_counts/file1_MATURE.read_count.txt'.replace('/', os.sep))

    def test_sample_hairpin_readcount_name(self):
        self.assertEqual(self.sample.hairpin_readcount, '/path/to/analysis_dir/read_counts/file1_HAIRPIN.read_count.txt'.replace('/', os.sep))


class TestSampleMethods(unittest.TestCase):
    raw_path = '/some/path/to/a/file1.fastq'.replace('/', os.sep)
    analysis_dir = '/path/to/analysis_dir'.replace('/', os.sep)
    sample = Sample(raw_path=raw_path, analysis_dir=analysis_dir)

    def test__get_basename(self):
        self.assertEqual(Sample._get_basename('/some/path/file1.fastq'.replace('/', os.sep)), 'file1')

    def test__create_filepath_no_file_no_dir(self):
        filepath = self.sample._create_filepath('.descriptor.extension')
        self.assertEqual(filepath, os.path.join(self.analysis_dir, 'file1.descriptor.extension'))

    def test__create_filepath_with_file_no_dir(self):
        file = '/new/path/to/new_file.fastq'
        filepath = self.sample._create_filepath('.descriptor.extension', file=file)
        self.assertEqual(filepath, os.path.join(self.analysis_dir, 'new_file.descriptor.extension'))

    def test__create_filepath_with_file_with_dir(self):
        file = '/new/path/to/new_file.fastq'
        dir = '/new/dir'
        filepath = self.sample._create_filepath('.descriptor.extension', file=file, dir=dir)
        self.assertEqual(filepath, os.path.join(dir, 'new_file.descriptor.extension'))

    def test_remove_intermediates(self):
        sample = Sample('tmp/file.fastq'.replace('/', os.sep), 'tmp/analysis_dir'.replace('/', os.sep))
        os.makedirs('tmp/analysis_dir'.replace('/', os.sep), exist_ok=True)

        trimmed_file = 'tmp/analysis_dir/file.trimmed.fastq'.replace('/', os.sep)
        with open(trimmed_file, 'w') as f:
            f.write('All these bases')

        self.assertTrue(os.path.isfile(trimmed_file))
        self.assertIsNone(sample.remove_intermediates())
        self.assertFalse(os.path.isfile(trimmed_file))

        shutil.rmtree('tmp'.replace('/', os.sep), ignore_errors=True)

    def test_read_counts_exist(self):
        sample = Sample('tmp/file.fastq'.replace('/', os.sep), 'tmp/analysis_dir'.replace('/', os.sep))
        os.makedirs('tmp/analysis_dir'.replace('/', os.sep), exist_ok=True)

        mature_read_count = 'tmp/analysis_dir/read_counts/file_MATURE.read_count.txt'.replace('/', os.sep)
        with open(mature_read_count, 'w') as f:
            f.write('mature read counts')
        hairpin_read_count = 'tmp/analysis_dir/read_counts/file_HAIRPIN.read_count.txt'.replace('/', os.sep)
        with open(hairpin_read_count, 'w') as f:
            f.write('hairpin read counts')

        self.assertTrue(sample.read_counts_exist())
        self.assertTrue(sample.read_counts_exist(mature_only=True))
        self.assertTrue(sample.read_counts_exist(hairpin_only=True))
        with self.assertRaises(ValueError):
            self.sample.read_counts_exist(mature_only=True, hairpin_only=True)

        shutil.rmtree('tmp'.replace('/', os.sep), ignore_errors=True)
        self.assertFalse(sample.read_counts_exist())
        self.assertFalse(sample.read_counts_exist(mature_only=True))
        self.assertFalse(sample.read_counts_exist(hairpin_only=True))

    def test_change_read_count_dir(self):
        self.sample.change_read_count_dir('/new/readcount_dir'.replace('/', os.sep))
        self.assertEqual(self.sample.mature_readcount, '/new/readcount_dir/file1_MATURE.read_count.txt'.replace('/', os.sep))
        self.assertEqual(self.sample.hairpin_readcount, '/new/readcount_dir/file1_HAIRPIN.read_count.txt'.replace('/', os.sep))


if __name__ == '__main__':
    unittest.main()
