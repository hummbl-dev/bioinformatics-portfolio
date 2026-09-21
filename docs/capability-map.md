# Capability Map — artifacts ↔ program requirements

Public companion to the private programs matrix. Each artifact answers
the question an admissions committee asks: *has this person already done
the work our first-years do?*

## Unified spec (distilled from 20 programs' published requirements)

| Domain | What programs require | Proven by |
|---|---|---|
| D1 Algorithms | Sequence algorithms course (Stanford 214, UCSD 282, Duke 561, NYU 7453, MIT core) | `algorithms/` — implementations + benchmarks |
| D2 NGS pipelines | NGS analysis cores (NYU 7653, GA Tech 6150/7210, Georgetown 5990) | `pipelines/ngs-variant/` — QC→call→annotate w/ concordance |
| D3 Programming | Python/R/UNIX literacy everywhere | All repos + containers + CI |
| D4 Biostatistics | Stats-for-genomics cores (BU MA 681, UNC 720, UW 560, Emory 510, Penn 5330) | Flagship stats section + `ml/` baselines |
| D5 Biology fluency | Genetics/mol-bio cores (Harvard GENETICS 201, Duke prereqs, Baylor 6600) | Flagship interpretation + `docs/literature-notes/` |
| D6 ML for biology | ML cores/tracks (Stanford, MIT, Emory 534, Penn 5200, GA Tech CX 4803) | `ml/` project w/ baselines + error analysis |
| D7 Track domains | Proteomics/popgen/translational/transcriptomics electives | Flagship extension + manifest candidates |
| D8 Research capability | Rotations → quals (often grant format) → dissertation; UNC expects 1st-author pub; Emory wants a technical report | Flagship → bioRxiv + 3-aim proposal + OSS PRs |

## Program-specific notes

- **UNC BCB**: "1st-author publication" appears in degree progression —
  the bioRxiv preprint is the closest pre-admission analog.
- **UW Genome Sciences**: general exam = 6-page, 3-aim written proposal
  in F31 format. `docs/` will carry a matching proposal.
- **Emory BMI**: qualifying component is a *publication-quality
  technical report* — the reproduction writeup is a dry run.
- **Duke CBB**: published prereq list (diff eq, lin alg, stats,
  algorithms, genetics, mol bio, biochem) is a literal self-audit
  checklist — each item maps to an artifact or transcript line.
- **Georgia Tech MS (PSM)**: GRA awards are competitive after semester
  1 — arriving with a public portfolio is direct leverage.
- **Umbrella applications**: UNC→BBSP, Michigan→PIBS, Harvard→HILS
  ("Biomedical Informatics" → BIG track).

## What this portfolio deliberately does not claim

- Wet-lab competence — no substitute for bench rotations (BU/CPCB/MIT
  require experimental rotations; the portfolio shows *literacy*, not
  bench skill).
- Clinical/protected-data experience — all work is on public data;
  IRB-regulated work is post-admission territory.
