import argparse
import csv
import re

import pysam


def analyse_sam(sam_records, read_names):
    read_list = []

    read_names_list = []
    with open(read_names) as f:
        for line in f:
            read_names_list.append(line.strip())

    cons_read_id = 0
    past_cons_read_name = ""
    for record in sam_records:
        if record.is_supplementary:
            continue

        alignment_length = 0
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

        cigar_stats = record.get_cigar_stats()

        read = {
            "name": read_name,
            "seq_name": "",
            "seq_family": seq_family,
            "read_type": read_type,
            "reference_length": read_length,
            "alignment_length": alignment_length,
            "alignment_type": alignment_type,
            "overlap": overlap,
            "overlap_percentage": overlap / read_length * 100 if read_length != 0 else 0,
            "mapping_quality": record.mapping_quality,
            "alignment_score": record.get_tag("AS") if record.has_tag("AS") else -1,
            "gap-compressed_seq_divergece": record.get_tag("de") if record.has_tag("de") else -1,
            "insertion_base_count": cigar_stats[0][1],
            "deletion_base_count": cigar_stats[0][1],
            "match_base_count": cigar_stats[0][7],
            "mismatch_base_count": cigar_stats[0][8],
            "mismatch_indel_base_count": cigar_stats[0][10],
            "insertion_event_count": cigar_stats[1][1],
            "deletion_event_count": cigar_stats[1][1],
            "match_event_count": cigar_stats[1][7],
            "mismatch_event_count": cigar_stats[1][8],
            "read_identity_match_mismatch": cigar_stats[0][7]/(cigar_stats[0][7] + cigar_stats[0][8]) * 100 if alignment_type.startswith("mapped") else 0,
            "read_identity_match_mismatch_gap": cigar_stats[0][7] / (cigar_stats[0][7] + cigar_stats[0][10]) * 100 if alignment_type.startswith("mapped") else 0
        }

        read_list.append(read)

    return read_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("sam", help="SAM file of simulated reads")
    parser.add_argument("read_names", help="Text file containing the true read names of the processed reads")
    parser.add_argument("output_prefix", help="Output prefix for info matrices")
    args = parser.parse_args()

    read_list = analyse_sam(pysam.AlignmentFile(args.sam), args.read_names)

    with open(args.output_prefix + "_identity_info.csv", "w") as read_info_out:
        dict_writer = csv.DictWriter(read_info_out, fieldnames=["name",
                                                                "seq_name",
                                                                "seq_family",
                                                                "read_type",
                                                                "reference_length",
                                                                "alignment_length",
                                                                "alignment_type",
                                                                "overlap",
                                                                "overlap_percentage",
                                                                "mapping_quality",
                                                                "alignment_score",
                                                                "gap-compressed_seq_divergece",
                                                                "insertion_base_count",
                                                                "deletion_base_count",
                                                                "match_base_count",
                                                                "mismatch_base_count",
                                                                "mismatch_indel_base_count",
                                                                "insertion_event_count",
                                                                "deletion_event_count",
                                                                "match_event_count",
                                                                "mismatch_event_count",
                                                                "read_identity_match_mismatch",
                                                                "read_identity_match_mismatch_gap"])
        dict_writer.writeheader()
        for read in read_list:
            dict_writer.writerow(read)

