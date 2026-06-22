import csv


def read_fasta(path):
    records = {}
    rid = None
    seq = []

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if rid:
                    records[rid] = "".join(seq)

                rid = line[1:].split(" ")[0]
                seq = []
            else:
                seq.append(line.strip())

    if rid:
        records[rid] = "".join(seq)

    return records


def write_perbase_logprobs_tsv(rows, output_path):
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")

        writer.writerow(
            [
                "record_tag",
                "allele",
                "pos_index_0based",
                "rel_to_center",
                "base",
                "logprob",
            ]
        )

        for row in rows:
            writer.writerow(
                [
                    row["record_tag"],
                    row["allele"],
                    row["pos_index_0based"],
                    row["rel_to_center"],
                    row["base"],
                    f"{row['logprob']:.6f}",
                ]
            )