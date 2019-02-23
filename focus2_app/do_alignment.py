# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict

from pysam import idxstats


def parse_alignments(alignments_path):
    """Using the pysam wrapper for the function idxstats, this function creates a dict with the counts per reference.

    Args:
        alignments_path (str): Path to SAM/BAM file.

    Returns:
        dict: counts per reference.

    """
    abundance_counts = defaultdict(int)

    for row in idxstats(alignments_path).split("\n"):
        row = row.split()
        if len(row) == 4 and row[0] != "*":
            reference_name = row[0]
            counts = int(row[2])
            abundance_counts[reference_name] = counts

    return abundance_counts
