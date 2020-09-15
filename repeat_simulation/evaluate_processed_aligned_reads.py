import argparse
import csv
import re

import pysam
from collections import defaultdict

import intervaltree


def analyse_sam(sam_records, read_names, full_reads=None):
    read_list = []
    seq_mapping_summary = defaultdict(lambda: defaultdict(int))

    read_names_list = []
    with open(read_names) as f:
        for line in f:
            read_names_list.append(line.strip())

    cons_read_id = 0
    past_cons_read_name = ""
    for record in sam_records:
        if record.is_supplementary:
            continue

        reference_chr = ""
        reference_pos = [0, 0]
        alignment_chr = ""
        alignment_pos = [0, 0]
        alignment_length = 0
        alignment_relative_pos = [0.0, 0.0]
        overlap = 0
        read_length = 0
        alignment_type = "unmapped"

        if past_cons_read_name == "":
            past_cons_read_name = record.query_name
        else:
            if past_cons_read_name != record.query_name:
                cons_read_id += 1
                past_cons_read_name = record.query_name
        read_name = read_names_list[cons_read_id]

        if read_name != "MIXED":
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

                if overlap > 0 and read_type == "full":
                    alignment_relative_pos = [(alignment_pos[i] - reference_pos[i]) / read_length for i in range(2)]

                if alignment_length <= read_length * 1.1:
                    if overlap > 0:
                        alignment_type = "mapped"
                    else:
                        alignment_type = "mismapped"
        else:
            seq_family = "mixed"
            read_type = "unresolved"

            if record.is_unmapped:
                alignment_type = "unmapped_unresolved"
            else:
                alignment_type = "mapped_unresolved"

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
            "overlap_percentage": overlap / read_length * 100 if read_length != 0 else 0,
        }

        read_list.append(read)

    if full_reads:
        seq_intervaltrees = defaultdict(intervaltree.IntervalTree)  # For full length reads only

        with open(full_reads) as f:
            for line in f:
                read_name = line.strip()
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

                full_read = {
                    "seq_name": read_name.rsplit("_", 2)[0],
                    "seq_family": seq_family,
                    "reference_chr": reference_chr,
                    "reference_start": reference_pos[0],
                    "reference_end": reference_pos[1],
                    "reference_length": read_length,
                }

                seq_intervaltrees[reference_chr].addi(reference_pos[0], reference_pos[1], full_read)

        for read in read_list:
            seq_info = None
            if read["name"] != "MIXED":
                seq_origins = list(seq_intervaltrees[read["reference_chr"]].overlap(read["reference_start"],
                                                                                    read["reference_end"]))

                for seq_origin in seq_origins:
                    if seq_origin.data["seq_family"] == read["seq_family"]:
                        seq_info = seq_origin.data

                if seq_info == None:
                    raise ValueError("Missing sequence origin for " + seq_origins)

                seq_name = seq_info["seq_name"]
                seq_family = seq_info["seq_family"]
            else:
                seq_name = "mixed"
                seq_family = "mixed"

            seq_mapping_summary[seq_name]["seq_family"] = seq_family
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

            if read["alignment_type"] == "mismapped" or read["alignment_type"] == "mapped_unresolved":
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
    parser.add_argument("read_names", help="Text file containing the true read names of the processed reads")
    parser.add_argument("output_prefix", help="Output prefix for info matrices")
    parser.add_argument("-f", "--full_reads", help="Text file containing the read names for full length reads")
    args = parser.parse_args()

    read_list, seq_mapping_summary = analyse_sam(pysam.AlignmentFile(args.sam), args.read_names, args.full_reads)

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

    if args.full_reads:
        with open(args.output_prefix + "_seq_info.csv", "w") as seq_info_out:
            dict_writer = csv.DictWriter(seq_info_out, fieldnames=["seq_name", "seq_family", "mapped",
                                                                   "mismapped", "unmapped", "mapped_unresolved",
                                                                   "unmapped_unresolved", "total",
                                                                   "false_positive_mapped"])
            dict_writer.writeheader()
            for seq_name in seq_mapping_summary:
                seq_info = seq_mapping_summary[seq_name]
                seq_info["seq_name"] = seq_name
                dict_writer.writerow(seq_info)

