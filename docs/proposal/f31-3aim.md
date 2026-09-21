# Convergent verification of cancer-driver biology through scripted, multi-arm re-analysis of public cohorts

**Format:** NIH F31 / general-exam style — Specific Aims (1 p.) +
Research Strategy (~5 pp). Modeled on UW Genome Sciences general
exam (6-page, 3-aim proposal) and UNC BCB written qualifying
structure.
**Author:** Reuben P. Bowlby · **Date:** 2026-09-21
**Preliminary data:** `analysis/flagship-reproduction/` — TCGA-BRCA
3-arm reproduction (stats + ML + CNA), all scripts committed.

---

## Specific Aims

Computational cancer genomics has a structural verification problem:
published driver–subtype associations are re-derived ad hoc in each
lab's private pipeline, and negative results of re-analysis are
almost never reported. The result is a literature whose canonical
findings are cited far more often than they are re-tested.

Our preliminary work demonstrates that this need not be true. Using
only the cBioPortal public API and ~500 lines of dependency-light
Python, we reproduced the TCGA 2012 subtype-enrichment associations
on the later PanCancer Atlas cohort (TP53-basal OR 27.7; PIK3CA-LumA
OR 3.25; GATA3 absent in basal), showed a supervised classifier
independently recovers exactly the subtypes the statistics identify —
and fails exactly where biology predicts mutation features should be
uninformative (Her2, F1 0.18) — and confirmed that copy-number data
closes the gap (ERBB2-AMP in Her2, OR 33.6). Three independent arms,
one biological story.

We propose to generalize this **convergent-verification framework**
into a reusable methodology and apply it across TCGA:

- **Aim 1 — Multi-cohort replication atlas.** Apply the scripted
  association pipeline to 10 additional TCGA PanCancer Atlas studies
  (LUAD, LUSC, COAD, STAD, HNSC, BLCA, OV, UCEC, KIRC, GBM) for each
  cancer type's published driver–subtype associations; produce a
  replication scorecard quantifying which canonical claims hold under
  uniform re-analysis.

- **Aim 2 — Convergent ML verification.** Extend the classifier arm to
  joint mutation + CNA + (where available) expression features; test
  the hypothesis that out-of-fold class-separability per subtype
  predicts which alteration *type* carries subtype signal — turning
  the Her2 observation into a general diagnostic.

- **Aim 3 — A reproducibility standard and open toolkit.** Package the
  pipeline as a versioned, containerized toolkit (one command →
  verified association tables for any cBioPortal study), with a
  reporting standard that captures cohort provenance, alteration-type
  definitions, and negative results. Success = an external user
  regenerates our TCGA-BRCA numbers from the tool alone.

**Impact.** A positive result reframes replication from a labor
problem to an infrastructure problem — and provides the field a
checklist for "has this association been independently re-derived?"

---

## Research Strategy

### A. Significance

Replication rates in computational biology are poorly characterized
precisely because re-analysis is labor-intensive and unrewarded. The
claims we target — driver–subtype associations — underpin subtype-
directed clinical reasoning (e.g., basal-like BRCA → platinum/PARP
sensitivity rationale via TP53/homologous-recombination context;
luminal PIK3CA → alpelisib eligibility). Verifying them cheaply and
continuously has direct translational value.

The significance is methodological as much as biological: if a
~500-line stdlib pipeline reproduces a *Nature* paper's core
associations, the barrier to routine replication is tooling, not
effort — and the remedy is a toolkit, not an exhortation.

### B. Innovation

1. **Convergent verification as method**: running statistics and ML
   as *independent arms* and treating their agreement — including
   agreement about *where signal should be absent* — as the
   verification criterion. To our knowledge this triangulation has
   not been formalized as a replication standard.
2. **Negative-result capture as a first-class output**: the scorecard
   reports non-replicating claims with the same rigor as replicating
   ones.
3. **Zero-credential, zero-infrastructure replication**: everything
   runs against public APIs in stdlib-light code — auditable by a
   single reviewer in one sitting.

### C. Approach

#### Aim 1 — Multi-cohort replication atlas

*Design.* For each of 10 Atlas studies, pre-register the specific
published associations being tested (source paper, claimed direction
and approximate magnitude) before running the pipeline — preventing
post-hoc claim selection. Associations drawn from each cancer type's
TCGA marker paper plus Atlas follow-ups.

*Pipeline.* Generalize the BRCA scripts: study discovery → clinical
attribute resolution (subtype variable names differ by study) →
alteration fetch (mutation + GISTIC) → Fisher/BH battery →
scorecard JSON + auto-generated table.

*Analysis.* Each claim scored: replicated (direction + q<0.05),
direction-consistent (direction right, underpowered), or failed.
Inter-rater check: a subset re-derived via the raw GDC/TCGA API to
test cBioPortal-dependence.

*Pitfalls.* Subtype nomenclature varies across studies (mitigate:
attribute-value census first, explicit label maps committed to repo);
some claims reference proteomic/methylation subtypes outside the
API's mutation/CNA profiles (mitigate: pre-registration marks these
out-of-scope rather than silently dropping them).

*Expected outcome.* A published scorecard: for each cancer type,
which canonical associations survive uniform scripted re-analysis —
including, honestly reported, which do not.

#### Aim 2 — Convergent ML verification

*Design.* Extend the 46-gene panel model to a three-view feature
space (mutation, CNA, mRNA z-scores where the study provides them).
Per subtype, per cancer type: stratified out-of-fold macro-F1 per
feature view and for concatenated views.

*Hypothesis.* The view(s) that recover a subtype in silico are the
same views in which the subtype's defining alterations live (Her2 ↔
CNA; Basal ↔ mutation + CNA; luminal splits ↔ expression). If
confirmed, per-class classifier failure becomes a *diagnostic* for
which alteration type carries subtype signal — a general tool.

*Analysis.* Ablation per view; calibration curves; SHAP/permutation
importances compared against the Fisher-test gene ranking from Aim 1.
Biological readout validated against published subtype definitions.

*Pitfalls.* Expression profiles exist for fewer patients; joint-model
n will shrink (mitigate: report per-view n transparently, nested CV
for any hyperparameter choice). Multi-view classifiers can leak
through imputation (mitigate: no imputation; per-view models on
view-complete patients only).

*Expected outcome.* A rule: "subtype X is recoverable from view Y"
verified across cancer types — or a documented map of where the rule
breaks.

#### Aim 3 — Reproducibility standard + toolkit

*Design.* Package Aims 1–2 into `cgverify` (working name): one CLI
command takes a study ID + claim file → runs fetch/stats/ML → emits
a signed JSON scorecard + human-readable report. Versioned,
containerized (Dockerfile already in repo), CI-tested.

*Standard.* A `REPLICATION.md` template: cohort provenance,
alteration-type definitions, denominator rules, negative results,
deviations from the source paper — the fields our BRCA
REPRODUCTION.md already uses, formalized.

*Validation.* External-user test: a collaborator unfamiliar with the
codebase regenerates our BRCA numbers using only the tool + docs.
Failure to reproduce from the tool alone is treated as a defect in
the tool, not the user.

*Expected outcome.* A citable toolkit + reporting standard submitted
as a software paper (JOSS/F1000) accompanying the atlas results.

### Timeline (36 months, F31-shaped)

- Yr 1: Aim 1 pipeline generalization + 5/10 cohorts; Aim 3 skeleton
- Yr 2: Aim 1 complete + scorecard manuscript; Aim 2 model runs
- Yr 3: Aim 2 manuscript; Aim 3 toolkit + software paper; dissertation

### Preliminary data summary (this repository)

| Claim tested | Result | Arm |
|---|---|---|
| TP53 enriched in basal | 89.5%, OR 27.7, q<10⁻⁵⁰ ✓ | stats |
| PIK3CA enriched in LumA | 47.5%, OR 3.25 ✓ | stats |
| GATA3 absent in basal | 0/171 ✓ | stats |
| Her2 learnable from mutations | F1 0.18 — correctly unlearnable | ML |
| ERBB2-AMP defines Her2 | 70.5%, OR 33.6 ✓ | CNA |
| Flat learning curve → feature-limited | confirmed n=196→981 | ML |
