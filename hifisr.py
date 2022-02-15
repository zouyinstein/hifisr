#!/usr/bin/env python3
"""
Author : zouyi <zouyi.nju@gmail.com>
Date   : 2022-02-15
Purpose: The Python user interface of the HiFi-SR pipeline
"""

import argparse
from typing import NamedTuple
import os
import subprocess

class Args(NamedTuple):
    """ Command-line arguments """
    mode: str   # single or merge
    sample: str
    thread: int
    input_type: str
    type_number: int
    samples_txt_prefix: str

# --------------------------------------------------

def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='The Python user interface of the HiFi-SR pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('mode',
                        metavar='mode',
                        type=str,
                        help='"single" to run hifisr for each sample; "merge" to merge results of multiple samples')

    parser.add_argument('-s',
                        '--sample',
                        help='In "single" mode, sample name, for example sample_1',
                        metavar='str',
                        type=str,
                        default='')

    parser.add_argument('-t',
                        '--thread',
                        help='In "single" mode, number of threads to run the pipeline',
                        metavar='int',
                        type=int,
                        default=8)

    parser.add_argument('-i',
                        '--input_type',
                        help='In "single" mode, fastq or fasta',
                        metavar='str',
                        type=str,
                        default='fastq')

    parser.add_argument('-n',
                        '--type_number',
                        help='In "single" mode, analyze reads up to four rearrangements in a single read',
                        metavar='int',
                        type=int,
                        default=5)

    parser.add_argument('-m',
                        '--samples_txt_prefix',
                        help='In "merge" mode, prefix of txt file containing multiple sample names, for example "merge_1"',
                        metavar='str',
                        type=str,
                        default='')

    args = parser.parse_args()

    return Args(args.mode, args.sample, args.thread, args.input_type, args.type_number, args.samples_txt_prefix)

# --------------------------------------------------

def main() -> None:
    """ Make a jazz noise here """

    args = get_args()
    mode_arg = args.mode
    sample_arg = args.sample
    thread_arg = args.thread
    input_type_arg = args.input_type
    type_number_arg = args.type_number
    samples_txt_prefix_arg = args.samples_txt_prefix

    if mode_arg == "single" and sample_arg == "":
        print()
        print("Missing sample name")
        print()
        exit()
    elif mode_arg == "single":
        print()
        print(f"Running HiFi-SR pipeline for {sample_arg}")
        print()
        ret = subprocess.call(f"bash ../../scripts/run_stage_all.sh {sample_arg} {thread_arg} {input_type_arg}", shell=True)
        print()
        print(f"Finished running HiFi-SR for {sample_arg}!")
        print()

    if mode_arg == "merge" and samples_txt_prefix_arg == "":
        print()
        print("Missing txt file prefix")
        print()
        exit()
    elif mode_arg == "merge":
        print()
        print(f"Merge HiFi-SR reports for multiple samples in {samples_txt_prefix_arg}")
        print()
        ret = subprocess.call(f"bash ../scripts/merge_samples.sh {samples_txt_prefix_arg}", shell=True)
        print()  
        print(f"Finished merging {samples_txt_prefix_arg}!")
        print()

# --------------------------------------------------

if __name__ == '__main__':
    main()