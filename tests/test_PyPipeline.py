import os
import shutil
import unittest
from pypipeline import PyPipeline


class TestPyPipelineInit(unittest.TestCase):
    def setUp(self):
        self.pipeline = PyPipeline(config_file='resources/config.ini', no_prompts=True, testing=True)

    def test_init(self):
        self.assertTrue(self.pipeline)

    def test_sample_conditions(self):
        sample_conditions = {'SAMPLE1': 'control',
                             'SAMPLE2': 'stress'}
        self.assertEqual(self.pipeline.sample_conditions, sample_conditions)

    def test_samples(self):
        self.assertEqual(len(self.pipeline.samples), 2)

    def tearDown(self):
        shutil.rmtree('tmp/'.replace('/', os.sep), ignore_errors=True)


class TestRunCommand(unittest.TestCase):
    def setUp(self):
        self.pipeline = PyPipeline(config_file='resources/config.ini', no_prompts=True, testing=True)
        command = 'touch tmp/test.file'.replace('/', os.sep)
        message = 'Making test.file'
        self.pipeline._run_command(command, message)

    def test_command_log(self):
        with open(self.pipeline.command_log, 'r') as f:
            lines = f.read().splitlines()
            self.assertEqual('touch tmp/test.file'.replace('/', os.sep), lines[-1])

    def test_log(self):
        with open(self.pipeline.log_file, 'r') as f:
            lines = f.read().splitlines()
            last_line_message = lines[-1].split('] ')[1]
            self.assertEqual('Making test.file...', last_line_message)

    def test_file_made(self):
        self.assertTrue(os.path.isfile('tmp/test.file'.replace('/', os.sep)))

    def tearDown(self):
        shutil.rmtree('tmp/'.replace('/', os.sep), ignore_errors=True)


class TestMethods(unittest.TestCase):
    def setUp(self):
        self.pipeline = PyPipeline(config_file='resources/config.ini', no_prompts=True, testing=True)

    def test_no_command_status(self):
        self.pipeline._log_message('Test')
        with open(self.pipeline.log_file, 'r') as f:
            lines = f.read().splitlines()
            last_line_message = lines[-1].split('] ')[1]
            self.assertEqual('Test...', last_line_message)

    def test_command_status(self):
        self.pipeline._log_message('Test', command_status=self.pipeline.BAD)
        with open(self.pipeline.log_file, 'r') as f:
            lines = f.read().splitlines()
            last_line_message = lines[-1].split('] ')[1]
            self.assertEqual('Test...', last_line_message)

    def test_create_log_file(self):
        self.pipeline._create_log_file()
        self.assertTrue(os.path.isfile(self.pipeline.log_file))

    def test_create_summary_file(self):
        self.pipeline._create_summary_file()
        self.assertTrue(os.path.isfile(self.pipeline.summary_file))

    def test_validate_file_that_doesnt_exist(self):
        with self.assertRaises(SystemExit):
            self.pipeline._validate_file('file_that_doesnt_exist')

    def test_validate_file_that_does_exist(self):
        self.assertIsNone(self.pipeline._validate_file('resources/config.ini'))

    def test_existing_program(self):
        self.assertIsNone(self.pipeline._check_program('ls'))

    def test_nonexisting_program(self):
        with self.assertRaises(SystemExit):
            self.assertIsNone(self.pipeline._check_program('program_that_doesnt_exist'))

    def test_existing_index(self):
        path_to_index = os.path.join(self.pipeline.bowtie_dir, 'some_index')
        self.assertTrue(self.pipeline._index_is_built(path_to_index, 'some_index'))

    def test_bad_index(self):
        path_to_index = os.path.join(self.pipeline.bowtie_dir, 'bad_index')
        self.assertFalse(self.pipeline._index_is_built(path_to_index, 'bad_index'))

    def test_nonexisting_index(self):
        path_to_nonexisting_index = os.path.join(self.pipeline.bowtie_dir, 'some_non_existing_index')
        self.assertFalse(self.pipeline._index_is_built(path_to_nonexisting_index, 'some_non_existing_index'))

    def tearDown(self):
        shutil.rmtree('tmp/'.replace('/', os.sep), ignore_errors=True)


class TestConfigs(unittest.TestCase):
    def test_bad_sample_conditions(self):
        with self.assertRaises(SystemExit):
            PyPipeline(config_file='resources/config_bad_sample_conditions.ini', no_prompts=True, testing=True)

    def test_bad_sample_names(self):
        with self.assertRaises(SystemExit):
            PyPipeline(config_file='resources/config_bad_sample_names.ini', no_prompts=True, testing=True)

    def test_bc_config(self):
        PyPipeline(config_file='resources/config_bc.ini', no_prompts=True, testing=True)

    def test_not_bc_not_toronto_config(self):
        with self.assertRaises(SystemExit):
            PyPipeline(config_file='resources/config_not_bc_not_toronto.ini', no_prompts=True, testing=True)

    def tearDown(self):
        shutil.rmtree('tmp/'.replace('/', os.sep), ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
