# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import os


def which(program_name):
    """Python implementation of unix 'which' function.

    Args:
        program_name (str): Program name.

    Returns:
        None or str: Program path.

    """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program_name)
    if fpath:
        if is_exe(program_name):
            return program_name
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program_name)
            if is_exe(exe_file):
                return exe_file


def is_wanted_file(queries):
    """Remove files from query files that not have extension FASTA/FASTA/FNA.

    Args:
        queries (list): List with query names

    Returns:
        list: sorted list with only FASTA/FASTA/FNA files

    """
    queries = [query for query in queries if query.split(".")[-1].lower() in ["fna", "fasta", "fastq"]]
    queries.sort()

    return queries


def bwa_alignment(reference, input_reads, output):
    """Run bwa with either a single file or pair or files (FASTQs).

    Args:
        reference (str): Path to reference.
        input_reads (list): List with input reads path(s).
        output (str): Path to output file.

    Raises:
        Exception: More than two files given as input file.

    """
    # Single file
    if len(input_reads) == 1:
        os.system('bwa mem {} {} > {}'.format(reference, input_reads[0], output))

    # Pair files
    elif len(input_reads) == 2:
        os.system('bwa mem {} {} {} > {}'.format(reference, input_reads[0], input_reads[1], output))

    # not accepted
    else:
        raise Exception('Please input either one or a pair of files')


def bwa_index_db(reference_fasta):
    """Index reference and output it into same folder as reference

    Args:
        reference_fasta (str): Path to reference FASTA which will be formatted and outputted.

    """
    os.system('bwa index {}'.format(reference_fasta))


def samtools_view(alignment_file, parameter_flag_target, output_name):
    """Run samtools view.

    Args:
        alignment_file (str): Path to BAM/SAM
        flag_target (str): Parameter and flag target e.g -F 3844 (best hits) or -f 4 (unmapped flag)
        output_name (str): Path to output FASTQ file

    """
    os.system('samtools view -h {} {} > {}'.format(parameter_flag_target, alignment_file, output_name))


def samtools_bam2fq(alignment_file, output_name):
    """Run samtools bam2fq to get FASTQ from BAM/SAM.

    Args:
        alignment_file (str): Path to BAM/SAM
        output_name (str): Path to output FASTQ file
        .
    """
    os.system('samtools bam2fq {} > {}'.format(alignment_file, output_name))


def uncompress(input_file):
    """Uncompress ZIP file.

    Args:
        input_file (str): Input file to be uncompressed.

    Raises:
        Exception: unzip not installed.

    """
    # unzip is installed
    if which('unzip'):
        os.system('unzip {}'.format(input_file))

    else:
        raise Exception('Please install unzip')


def donwload(url):
    """Download database.

    Args:
        url (str): URL to be downloaded.

    Raises:
        Exception: neither wget nor curl is installed.

    """
    # wget is installed
    if which('wget'):
        os.system('wget {}'.format(url))

    # curl is installed
    elif which('curl'):
        os.system('curl -O {}'.format(url))
    else:
        raise Exception('Please install either wget or curl')
