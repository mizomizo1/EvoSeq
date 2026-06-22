from Bio import SeqIO


def read_fasta_as_dict(path):
    records = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        records[rec.id] = str(rec.seq).upper()
    return records


def check_fasta_pair(ref_records, mut_records):
    ref_ids = set(ref_records)
    mut_ids = set(mut_records)

    ref_only = ref_ids - mut_ids
    mut_only = mut_ids - ref_ids
    common = ref_ids & mut_ids

    if ref_only:
        raise ValueError(f"Reference-only IDs exist: { list(sorted(ref_only))[:5] }")
    if mut_only:
        raise ValueError(f"Mutant-only IDs exist: { list(sorted(mut_only))[:5] }")

    return common


def wrap_sequence(seq, width=80):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))