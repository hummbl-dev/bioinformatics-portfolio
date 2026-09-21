"""Generate a synthetic test dataset for the ngs-variant pipeline.

Creates:
  data/test/reference.fa      - 20 kb synthetic reference
  data/test/truth.vcf         - the planted variants (ground truth)
  data/test/sample_R{1,2}.fastq - paired reads with the variants present

A correct pipeline run should recover a high fraction of truth.vcf
calls. Run: python pipelines/ngs-variant/scripts/make_test_data.py
"""

import random
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "data" / "test"
REF_LEN = 20_000
READ_LEN = 150
DEPTH = 30
N_VARIANTS = 40


def main():
    rng = random.Random(13)
    OUT.mkdir(parents=True, exist_ok=True)

    ref = "".join(rng.choice("ACGT") for _ in range(REF_LEN))

    # plant heterozygous SNVs on ~half the reads' worth of signal
    variants = {}
    positions = sorted(rng.sample(range(1000, REF_LEN - 1000), N_VARIANTS))
    for pos in positions:
        alt = rng.choice([b for b in "ACGT" if b != ref[pos]])
        variants[pos] = (ref[pos], alt)

    sample = list(ref)
    for pos, (_, alt) in variants.items():
        sample[pos] = alt
    sample = "".join(sample)

    def reads(seq, tag):
        """Unpaired reads spanning the sequence at ~DEPTH/2 coverage."""
        for i in range(0, len(seq) - READ_LEN,
                       max(1, READ_LEN * 2 // DEPTH)):
            yield f"@{tag}_{i}", seq[i:i + READ_LEN]

    # diploid: R1 drawn from ref, R2 from mutated haplotype
    with open(OUT / "sample_R1.fastq", "w") as f:
        for name, seq in reads(ref, "r1"):
            f.write(f"{name}\n{seq}\n+\n{'I' * len(seq)}\n")
    with open(OUT / "sample_R2.fastq", "w") as f:
        for name, seq in reads(sample, "r2"):
            f.write(f"{name}\n{seq}\n+\n{'I' * len(seq)}\n")

    with open(OUT / "reference.fa", "w") as f:
        f.write(">synth\n")
        for i in range(0, REF_LEN, 60):
            f.write(ref[i:i + 60] + "\n")

    with open(OUT / "truth.vcf", "w") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\n")
        for pos, (r, a) in variants.items():
            f.write(f"synth\t{pos + 1}\t.\t{r}\t{a}\t100\n")

    print(f"wrote {OUT}: ref {REF_LEN}bp, {N_VARIANTS} planted SNVs, "
          f"{len(list(reads(ref, 'x')))} reads/haplotype")


if __name__ == "__main__":
    main()
