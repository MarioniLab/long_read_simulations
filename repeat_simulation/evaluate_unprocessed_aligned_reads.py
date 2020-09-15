import argparse
import csv
import re

import pysam
from collections import defaultdict

import intervaltree


def analyse_sam(sam_records):
    read_list = []
    seq_intervaltrees = defaultdict(intervaltree.IntervalTree)  # For full length reads only
    seq_mapping_summary = defaultdict(lambda: defaultdict(int))

    for record in sam_records:
        if record.is_supplementary:
            continue

        read_name = record.query_name.split(",")[0]
        seq_family, reference_loc, read_length, read_type = re.findall(r"(.+)_(chr.+:[0-9-]+)_([0-9]+bp)_(.+)",
                                                                     read_name)[0]
        reference_chr, reference_pos = reference_loc.split(":")
        reference_chr = reference_chr if "_" not in reference_chr else reference_chr.split("_")[1]
        # We will need to deal with mismatch between alignment and reference for the scaffolds versioning
        if "." in reference_chr or "v" in reference_chr:
            sep = "." if "." in reference_chr else "v"
            reference_chr = reference_chr.rsplit(sep, 1)[0]
        reference_pos = sorted([int(i) for i in reference_pos.split("-")])
        read_length = int(read_length.rstrip("bp"))

        alignment_chr = ""
        alignment_pos = [0, 0]
        alignment_length = 0
        alignment_relative_pos = [0.0, 0.0]
        overlap = 0
        alignment_type = "unmapped"

        if not record.is_unmapped:
            alignment_chr = record.reference_name
            if "." in alignment_chr or "v" in alignment_chr:
                sep = "." if "." in alignment_chr else "v"
                alignment_chr = alignment_chr.rsplit(sep, 1)[0]
            alignment_pos = sorted([record.reference_start, record.reference_end])
            alignment_length = record.alen

            if reference_chr == alignment_chr:
                overlap = max(min([reference_pos[1], alignment_pos[1]]) - max([reference_pos[0], alignment_pos[0]]),
                              0)

            if alignment_length <= read_length * 1.1:
                if overlap > 0:
                    alignment_type = "mapped"
                    if read_type == "full":
                        alignment_relative_pos = [(alignment_pos[i] - reference_pos[i]) / read_length for i in range(2)]
                else:
                    alignment_type = "mismapped"

        read = {
            "name": read_name,
            "seq_name": "",
            "seq_family": seq_family,
            "read_type": read_type,
            "reference_chr": reference_chr,
            "reference_start": reference_pos[0],
            "reference_end": reference_pos[1],
            "reference_length": read_length,
            "alignment_seq_name": "",
            "alignment_chr": alignment_chr,
            "alignment_start": alignment_pos[0],
            "alignment_end": alignment_pos[1],
            "alignment_length": alignment_length,
            "alignment_relative_start": alignment_relative_pos[0],
            "alignment_relative_end": alignment_relative_pos[1],
            "alignment_type": alignment_type,
            "overlap": overlap,
            "overlap_percentage": overlap / read_length * 100,
        }

        if read_type == "full":
            read["seq_name"] = read_name.rsplit("_", 2)[0]
            seq_intervaltrees[reference_chr].addi(reference_pos[0], reference_pos[1], read)

        read_list.append(read)

    for read in read_list:
        seq_origins = list(seq_intervaltrees[read["reference_chr"]].overlap(read["reference_start"],
                                                                            read["reference_end"]))

        seq_info = None
        for seq_origin in seq_origins:
            if seq_origin.data["seq_family"] == read["seq_family"]:
                seq_info = seq_origin.data

        if seq_info == None:
            print(read, seq_origins)
            exit()

        seq_name = seq_info["seq_name"]

        seq_mapping_summary[seq_name]["seq_family"] = seq_info["seq_family"]
        seq_mapping_summary[seq_name][read["alignment_type"]] += 1
        seq_mapping_summary[seq_name]["total"] += 1
        read["seq_name"] = seq_name

        if read["alignment_type"] == "mapped":
            seq_reference_pos = [seq_info["reference_start"], seq_info["reference_end"]]
            seq_length = seq_info["reference_length"]
            alignment_pos = [read["alignment_start"], read["alignment_end"]]
            alignment_relative_pos = [(alignment_pos[i] - seq_reference_pos[i]) / seq_length for i in range(2)]
            read["alignment_start"] = alignment_relative_pos[0]
            read["alignment_end"] = alignment_relative_pos[1]

        if read["alignment_type"] == "mismapped":
            alignment_origins = list(seq_intervaltrees[read["reference_chr"]].overlap(read["alignment_start"],
                                                                                      read["alignment_end"]))

            for alignment_origin in alignment_origins:
                alignment_info = alignment_origin.data
                read["alignment_seq_name"] = alignment_info["seq_name"]
                seq_mapping_summary[alignment_info["seq_name"]]["false_positive_mapped"] += 1

    return read_list, seq_mapping_summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("sam", help="SAM file of simulated reads")
    parser.add_argument("output_prefix", help="Output prefix for info matrices")
    args = parser.parse_args()

    read_list, seq_mapping_summary = analyse_sam(pysam.AlignmentFile(args.sam))

    with open(args.output_prefix + "_reads_info.csv", "w") as read_info_out:
        dict_writer = csv.DictWriter(read_info_out, fieldnames=["name", "seq_name", "seq_family", "read_type",
                                                                "reference_chr", "reference_start", "reference_end",
                                                                "reference_length", "alignment_seq_name",
                                                                "alignment_chr", "alignment_start", "alignment_end",
                                                                "alignment_length", "alignment_relative_start",
                                                                "alignment_relative_end", "alignment_type", "overlap",
                                                                "overlap_percentage"])
        dict_writer.writeheader()
        for read in read_list:
            dict_writer.writerow(read)

    with open(args.output_prefix + "_seq_info.csv", "w") as seq_info_out:
        dict_writer = csv.DictWriter(seq_info_out, fieldnames=["seq_name", "seq_family", "mapped",
                                                               "mismapped", "unmapped", "total",
                                                               "false_positive_mapped"])
        dict_writer.writeheader()
        for seq_name in seq_mapping_summary:
            seq_info = seq_mapping_summary[seq_name]
            seq_info["seq_name"] = seq_name
            dict_writer.writerow(seq_info)

