# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import logging
import os
import random

from pathlib import Path
from collections import defaultdict

from focus2_app import version


LOGGER_FORMAT = '[%(asctime)s - %(levelname)s] %(message)s'
logging.basicConfig(format=LOGGER_FORMAT, level=logging.INFO)
LOGGER = logging.getLogger(__name__)





def parse_args():
    """Parse args entered by the user.

    Returns:
        argparse.Namespace: parsed arguments

    """
    parser = argparse.ArgumentParser(
        description='FOCUS2: Agile and sensitive classification of metagenomics data using '
                    'a reduced database', epilog='example > focus2 -q INPUT_DIR -o OUTPUT')
    parser.add_argument('-v', '--version', action='version', version='FOCUS {}'.format(version))
    parser.add_argument('-q', '--query', help='Path to directory with FAST(A/Q) files', required=True)
    parser.add_argument('-o', '--output_directory', help='Path to output files', required=True)
    parser.add_argument('-p', '--output_prefix', help='Output prefix (Default: output)', default='output')
    parser.add_argument("-b", "--alternate_directory", help="Alternate directory for your databases", default="")
    parser.add_argument('-t', '--threads', help='Number Threads used during alignment (Default: 4)', default='4')

    return parser.parse_args()


def main():
    args = parse_args()

    # parameters and other variables
    query = Path(args.query)
    prefix = args.output_prefix
    output_directory = Path(args.output_directory)
    WORK_DIRECTORY = Path(args.alternate_directory) if args.alternate_directory else Path(__file__).parents[0]
    threads = args.threads
    aligner_path = which("bwa")

    LOGGER.info("FOCUS2: Agile and sensitive classification of metagenomics data using a reduced database")

    # check if output_directory is exists - if not, creates it
    if not output_directory.exists():
        Path(output_directory).mkdir(parents=True, mode=511)
        LOGGER.info("OUTPUT: {} does not exist - just created it :)".format(output_directory))

    # check if at least one of the queries is valid
    if not query.is_dir():
        LOGGER.critical("QUERY: {} is not a directory".format(query))

    # check if at least one of the queries is valid
    if is_wanted_file(os.listdir(query)) == []:
        LOGGER.critical("QUERY: {} does not have any FASTA/FNA/FASTQ file".format(query))

    # check if query is exists
    elif not query.exists():
        LOGGER.critical("QUERY: {} does not exist".format(query))

    # check if work directory exists
    elif WORK_DIRECTORY != WORK_DIRECTORY or not WORK_DIRECTORY.exists():
        LOGGER.critical("WORK_DIRECTORY: {} does not exist".format(WORK_DIRECTORY))


    else:
        LOGGER.info("1) RUNNING FOCUS2")

    LOGGER.info('Done'.format(output_directory))


if __name__ == "__main__":
    main()
