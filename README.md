# mirna-pipeline

This tool serves as a preproccessing pipeline from raw miRNAseq fastq files to the
RBioMIR suite of analysis tools producing analysis of known miRNA in publication ready
figures.

![Screenshot](./screenshot.jpg?raw=true)

## Installation
Use git to download the package.
```bash
>>> git clone https://github.com/liamhawkins/mirna-pipeline.git
>>> cd mirna-pipeline
>>> pip install -r requirements.txt
```
This pipeline also requires the following programs to be installed on your system:
```
* fastqc
* fastq-mcf
* cutadapt
* bowtie-build
* bowtie
* samtools
* Rscript
```

## Usage
Create a config file (See `example_config.ini` for exact template) for each set of
analysis you wish to process.

The pipeline can then be run from the command line:
```bash
>>> pypipeline.py --config example_config.ini
```
#### Process and analyze multiple data sets
Multiple config files defining multiple analysis can
be run in sequence by supplying a directory containing config *.ini files:
```bash
>>> pypipeline.py --config-dir dir_containing_configs/
```
In this case it is useful to suppress user prompts with the `--no-prompts` flag:
```bash
>>> pypipeline.py --no-prompts --config-dir dir_containing_configs/
```
#### Performing analysis only
If read counts are already available, you can perform the R analysis only using the
`--analysis-only` flag:
```bash
>>> pypipeline.py --config example_config.ini --analysis-only dir_with_readcounts/
```
Readcount file names need to be in the following format:
`<sample_name_from_config>_MATURE.read_count.txt`
#### All command line arguments
A full list of command line options can be found using the help flag:
```bash
>>> pypipeline.py --help
usage: pypipeline.py [-h] [-c <config_file> | -d <config_dir>] [--no-prompts]
                     [--no-fastqc] [--delete]
                     [--no-analysis | --analysis-only <read_count_dir>]

optional arguments:
  -h, --help            show this help message and exit
  -c <config_file>, --config <config_file>
                        Path to config file
  -d <config_dir>, --config-dir <config_dir>
                        Directory containing config files
  --no-prompts          Suppress user prompts
  --no-fastqc           Do not perform fastqc on raw files
  --delete              Delete intermediate processing files
  --no-analysis         Do not perform R analysis
  --analysis-only <read_count_dir>
                        Run analysis only on read counts in supplied directory
```

## LICENSE
[Link](https://choosealicense.com/licenses/mit/)
