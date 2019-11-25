# mirna-pipeline

This tool serves as a preproccessing pipeline from raw miRNAseq fastq files to the
RBioMIR suite of analysis tools producing analysis of known miRNA in publication ready
figures.

![Screenshot](./screenshot.jpg?raw=true)

## Installation
Use git to download the package.
```bash
git clone https://github.com/liamhawkins/mirna-pipeline.git
cd mirna-pipeline
pip install -r requirements.txt
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
python pypipeline.py --config example_config.ini
```
#### Process and analyze multiple data sets
Multiple config files defining multiple analysis can
be run in sequence by supplying paths of each config file to the `configs` variable
in `pypipeline.py`. i.e:
```python
# Path to config can be passed as command line arguments
if args.config:
    configs = [args.config]
else:
    configs = ['config1.ini', 'config2.ini', 'config3.ini']
```
In this case it is useful to suppress user prompts with the `--no-prompts` flag:
```bash
python pypipeline.py --no-prompts
```

## LICENSE
[Link](https://choosealicense.com/licenses/mit/)
