#!/usr/bin/env python3

import math
import re
import argparse
import random
from Bio import SeqIO, Seq


CUTOFF = 2000  # Ignore reads shorter than this cutoff value (2kb)
UMI_PATTERN = "NRY"


def get_umi(umi_length=22):
    pattern = (UMI_PATTERN * math.ceil(umi_length/len(UMI_PATTERN)))[:umi_length]
    bases = {"N": ["A", "G", "C", "T"],
             "R": ["A", "G"],
             "Y": ["C", "T"]}
    umi = []
    for base_type in pattern:
        umi += random.choices(bases[base_type])

    return "".join(umi)


def get_adapter(umi_length=22):
    cell_barcode = "CACAAAGACACCGACAACTTTCTT"
    adapter_template = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG{cell_barcode}AGTGGTATC{umi}ACTGGCCGTCGTTTTACATGGC" \
                        "GTAGCGGGTTCGAGCGCACCGCAGGGTATCCGGCTATTTTTTTTTTTTTTT".rstrip("T")

    return Seq.Seq(adapter_template.format(cell_barcode=cell_barcode, umi=get_umi(umi_length)))


def generate_reads(records, read_lengths, window_overlap, add_quality=False, umi_length=22, add_adapter=False):
    reads = {"full": []}
    for read_length in read_lengths:
        reads[read_length] = []

    read_lengths.append(float("Inf"))

    for record in records:
        seq_length = len(record)

        if seq_length < CUTOFF:
            continue

        sequence_id, chromosome_location, *_ = record.description.split()
        chrom, seq_start, seq_end = re.search(r"(chr[^:]+):([0-9]+)-([0-9]+)", chromosome_location).groups()
        seq_start = int(seq_start)
        seq_end = int(seq_end)
        strand = re.search(r"strand=([+-])", record.description).groups()[0]

        for read_length in read_lengths:
            read_dict_key = read_length
            if read_length == float("Inf"):
                read_dict_key = "full"
                read_length = seq_length

            if read_length > seq_length:
                continue

            number_of_reads = math.ceil((seq_length - (read_length - window_overlap)) / window_overlap)

            # We will generate from the 3' end following the likely capture scenario of RNA molecules
            for i in range(number_of_reads):
                end = seq_length - (i * window_overlap)
                start = end - read_length

                if start < 0:  # If we go slightly overboard, correct the start to actual sequence start
                    start = 0
                    end = read_length

                read_sequence = record.seq[start:end]
                if strand == "+":
                    read_start = seq_start + start
                    read_end = read_start + read_length - 1
                    read_position = (read_start, read_end)
                    read_strand = "+strand"

                else:
                    read_start = seq_end - start
                    read_end = read_start - read_length + 1
                    read_position = (read_end, read_start)
                    read_strand = "-strand"

                if add_adapter:
                    if strand == "+":
                        read_sequence = read_sequence.upper() + get_adapter().reverse_complement()
                    else:
                        read_sequence = get_adapter(umi_length) + read_sequence.upper()

                if read_dict_key == "full":
                    read_strand = "full"

                read_info = f">{sequence_id}_{chrom}:{read_position[0]}-{read_position[1]}_" \
                            f"{read_length}bp_{read_strand}"

                if not add_quality:
                    reads[read_dict_key].append(f"{read_info}\n{read_sequence}")
                else:
                    reads[read_dict_key].append(f"{read_info}\n{read_sequence}\n+\n{'I' * len(read_sequence)}")

    return reads


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="FASTA file for simulating reads")
    parser.add_argument("output_prefix", help="Output prefix for simulated reads")
    parser.add_argument("--read_lengths", help="Read lengths to be simulated "
                                               "(multiple read lengths should be specified with comma)",
                        default="1000,2000,3000")
    parser.add_argument("--window_overlap", help="Length of overlap between read window", default=500, type=int)
    parser.add_argument("--output_type", help="Type of output - fastq or fasta file", default="fastq",
                        choices=["fastq", "fasta"])
    parser.add_argument("--umi_length", help="UMI length to be simulated with the set pattern", default=22, type=int)
    parser.add_argument("--add_adapter", help="Include read adapters for simulated reads", action='store_true')
    args = parser.parse_args()

    read_lengths = []
    for read_length in args.read_lengths.split(","):
        if read_length.isdigit():
            read_lengths.append(int(read_length))
        else:
            raise ValueError(f"Invalid read length specified - {read_length}")

    reads = generate_reads(list(SeqIO.parse(args.fasta, "fasta")),
                           read_lengths=read_lengths,
                           window_overlap=args.window_overlap,
                           add_quality=True if args.output_type == "fastq" else False,
                           umi_length=args.umi_length,
                           add_adapter=args.add_adapter)
    print(len(reads["full"]))

    for read_length in read_lengths:
        if read_length != float("Inf"):
            with open(f"{args.output_prefix}_{read_length}bp.{args.output_type}", "w") as f:
                f.write("\n".join(reads[read_length] + reads["full"]))
