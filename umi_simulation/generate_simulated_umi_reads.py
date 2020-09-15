#!/usr/bin/env python3

import math
import re
import argparse
import random


UMI_PATTERN = "NNN"
BASES = {"N": ["A", "G", "C", "T"],
         "R": ["A", "G"],
         "Y": ["C", "T"]}


def get_umi(umi_length=22):
    pattern = (UMI_PATTERN * math.ceil(umi_length/len(UMI_PATTERN)))[:umi_length]

    umi = []
    for base_type in pattern:
        umi += random.choices(BASES[base_type])

    return "".join(umi)


def generate_umi_reads(umi_count, umi_lengths, add_quality=False):
    umi_reads = {}

    for umi_length in umi_lengths:
        umi_reads[umi_length] = []

        if add_quality:
            quality_score = f"\n+\n{'I' * umi_length}"
        else:
            quality_score = ""

        for _ in range(umi_count):
            umi_sequence = get_umi(umi_length)
            umi_reads[umi_length].append(f">UMI-{umi_sequence}\n{umi_sequence}"+quality_score)

    return umi_reads

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("umi_count", help="Number of UMI reads to generate", type=int)
    parser.add_argument("output_prefix", help="Output prefix for simulated UMI reads")
    parser.add_argument("--umi_lengths", help="UMI lengths to be simulated "
                                               "(multiple read lengths should be specified with comma)",
                        default="22")
    parser.add_argument("--output_type", help="Type of output - fastq or fasta file", default="fastq",
                        choices=["fastq", "fasta"])
    args = parser.parse_args()

    umi_lengths = []
    for umi_length in args.umi_lengths.split(","):
        if umi_length.isdigit():
            umi_lengths.append(int(umi_length))
        else:
            raise ValueError(f"Invalid read length specified - {umi_length}")

    umi_reads = generate_umi_reads(umi_count=args.umi_count,
                                   umi_lengths=umi_lengths,
                                   add_quality=True if args.output_type == "fastq" else False)

    for umi_length in umi_reads:
        with open(f"{args.output_prefix}_{umi_length}bp_umi.{args.output_type}", "w") as f:
            f.write("\n".join(umi_reads[umi_length]))
