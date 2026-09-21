"""Minimal FASTA/FASTQ readers — just enough for pipelines and tests."""


def read_fasta(path):
    """Yield (record_id, sequence) for each entry in a FASTA file."""
    rid, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if rid is not None:
                    yield rid, "".join(seq)
                rid = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.upper())
    if rid is not None:
        yield rid, "".join(seq)


def write_fasta(records, path, wrap=60):
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n")
            for i in range(0, len(seq), wrap):
                fh.write(seq[i:i + wrap] + "\n")
