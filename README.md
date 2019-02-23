#### FOCUS2: Agile and sensitive classification of metagenomics data using a reduced database
* [Installation](#installation)
    * [dependencies](#dependencies)
    * [pip3](#pip3)
    * [bioconda](#bioconda)
    * [git](#git)
* [Usage](#usage)
* [Output](#output)
* [Citing](#citing)

## Installation
### Dependencies
- [Python 3.6](http://www.python.org/download)
- [Setuptools 36.0.1](https://setuptools.readthedocs.io/en/latest/)
- [Jellyfish 2.2.6](https://github.com/gmarcais/Jellyfish/releases/tag/v2.2.6). if using macOS, use [bioconda](https://anaconda.org/bioconda/jellyfish)
- [Numpy 1.12.1](https://github.com/numpy/numpy)
- [SciPy 0.19.0](https://github.com/scipy/scipy)
- [Pysam](https://github.com/pysam-developers/pysam)
- unzip/curl

### pip3 (WIP
	pip3 install metagenomics-focus2

### Bioconda (WIP)
You can now easily install FOCUS2 using [conda](https://conda.io) via the
[Bioconda](https://bioconda.github.io/) channel. It is as easy as:

    # bioconda should handle all the dependencies
    conda create -n focus2 -c bioconda focus2
	source activate focus2

This will create a conda environment called `focus2` (as specified by the
`-n` argument), and install FOCUS2 along with all its dependencies. The second
line activates the newly created `focus2` conda environment.

### Git

	# clone focus2
	git clone git@github.com:metageni/FOCUS2.git

	# install focus2
	cd FOCUS2 && python setup.py install

## Usage

	focus2 [-h] [-v] -q QUERY -o OUTPUT_DIRECTORY [-p OUTPUT_PREFIX]
                 [-b ALTERNATE_DIRECTORY] [-t THREADS] [-l LOG]

    FOCUS2: Agile and sensitive classification of metagenomics data using a
    reduced database

    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
      -q QUERY, --query QUERY
                            Path to directory with FAST(A/Q) files
      -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                            Path to output files
      -p OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                            Output prefix (Default: output)
      -b ALTERNATE_DIRECTORY, --alternate_directory ALTERNATE_DIRECTORY
                            Alternate directory for your databases
      -t THREADS, --threads THREADS
                            Number Threads used during alignment (Default: 4)
      -l LOG, --log LOG     Path to log file (Default: STDOUT).

    example > focus2 -q INPUT_DIR -o OUTPUT


## Output (WIP)
FOCUS2 generates a tabular output per taxonomic level (`Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, `Species`, and `Strain`) and one with all levels which can be used as [STAMP](http://kiwi.cs.dal.ca/Software/STAMP)'s input for statistical analysis.

## Citing (WIP)
FOCUS2 was written by Genivaldo G. Z. Silva. Feel free to [contact me](mailto:genivaldo.gueiros@gmail.com)

If you use FOCUS2, please cite it:

    Silva, G. G. Z.,
