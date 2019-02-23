# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging
import os

from pathlib import Path
from collections import defaultdict
from shutil import which

from focus2_app import version
from wrappers import is_wanted_file, bwa_alignment, samtools_view, samtools_bam2fq, mock_focus_args
from do_alignment import parse_alignments

import focus_app.focus as focus

LOGGER_FORMAT = '[%(asctime)s - %(levelname)s] %(message)s'


def parse_args():
    """Parse args entered by the user.

    Returns:
        argparse.Namespace: Parsed arguments.

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
    parser.add_argument('-l', '--log', help='Path to log file (Default: STDOUT).', required=False)

    return parser.parse_args()


def main():
    args = parse_args()

    # parameters and other variables
    query = Path(args.query)
    prefix = args.output_prefix
    output_directory = Path(args.output_directory)
    reference = "db/SPN032672.fasta"
    work_directory = Path(args.alternate_directory) if args.alternate_directory else Path(__file__).parents[0]
    threads = args.threads

    if args.log:
        logging.basicConfig(format=LOGGER_FORMAT, level=logging.INFO, filename=args.log)
    else:
        logging.basicConfig(format=LOGGER_FORMAT, level=logging.INFO)

    logger = logging.getLogger(__name__)

    logger.info("FOCUS2: Agile and sensitive classification of metagenomics data using a reduced database.\n")

    # check if bwa was installed
    if not which("bwa"):
        logger.critical("Install BWA !!!")

    # check if output_directory is exists - if not, creates it
    if not output_directory.exists():
        Path(output_directory).mkdir(parents=True, mode=511)
        logger.info("OUTPUT: {} does not exist - just created it :)".format(output_directory))

    # check if at least one of the queries is valid
    if not query.is_dir():
        logger.critical("QUERY: {} is not a directory".format(query))

    # check if at least one of the queries is valid
    if is_wanted_file(os.listdir(query)) == []:
        logger.critical("QUERY: {} does not have any FASTA/FNA/FASTQ file".format(query))

    # check if query is exists
    elif not query.exists():
        logger.critical("QUERY: {} does not exist".format(query))

    # check if work directory exists
    elif work_directory != work_directory or not work_directory.exists():
        logger.critical("work_directory: {} does not exist".format(work_directory))

    else:
        logger.info("Running FOCUS to quickly predict taxa in your query file(s)")

        kmer_size = "6"
        mock_args = mock_focus_args(query, threads, output_directory, kmer_size, args.log)

        focus_result = focus.main(mock_args)
        logger.info("Done with FOCUS\n")

        logger.info("Running FOCUS2")

        abundance_counts = defaultdict(dict)

        for counter, target_file in enumerate(os.listdir(query)):
            counter += 1

            logger.info("{}) Working on: {}".format(counter, target_file))

            logger.info('   {}.1) Aligning reads in {}'.format(counter, target_file))

            # input and output names
            input_reads = ["{}/{}".format(query, target_file)]
            alignment_output = "{}/out_{}.sam".format(output_directory, target_file)

            # run alignment
            bwa_alignment(reference, input_reads, alignment_output, threads)

            # get best hits
            logger.info('   {}.2) Parsing Best hits in file: {}'.format(counter, alignment_output))
            clean_alignment_output = "{}/clean_out_{}.sam".format(output_directory, target_file)
            samtools_view(alignment_output, "-F 3844", clean_alignment_output)

            # get unmapped alignment and create FASTQ
            logger.info('   {}.3) Generating FASTQ with unmapped reads: {}'.format(counter, alignment_output))
            unmapped_alignment_output = "{}/unmapped_out_{}.sam".format(output_directory, target_file)
            unmapped_reads_output = "{}/unmapped_out_{}.fastq".format(output_directory, target_file)
            samtools_view(alignment_output, "-f 4", unmapped_alignment_output)
            samtools_bam2fq(unmapped_alignment_output, unmapped_reads_output)

            # get abundance counts
            logger.info('   {}.4) Counting abundance in output file: {}'.format(counter, clean_alignment_output))
            abundance_counts[target_file] = parse_alignments(clean_alignment_output)

    for i in abundance_counts:
        print(i, abundance_counts[i])

    logger.info('Done'.format(output_directory))


if __name__ == "__main__":
    main()
