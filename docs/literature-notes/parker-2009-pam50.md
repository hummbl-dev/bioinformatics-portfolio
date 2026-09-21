# Note — Parker et al. (2009), *J Clin Oncol* 27:1160–1167

**Why this paper:** the subtype labels my entire flagship analysis
inherits come from this 50-gene assay.

## What they did

Derived a minimal 50-gene qPCR expression signature (PAM50) that
recovers the intrinsic subtypes (Luminal A/B, HER2-enriched,
Basal-like, Normal-like) from microarray-era centroids; validated
prognostic value (risk-of-relapse score) on independent cohorts.

## Key ideas

- Subtype is an *expression-space* object — a nearest-centroid
  classification, not a mutation-defined category.
- ROR score adds proliferation weighting → subtype + risk in one assay.
- Clinical adoption (Prosigna) made it the de facto subtype standard.

## Strengths

- Small, assayable gene set → clinically deployable.
- Consistent centroid structure across platforms/cohorts.

## Caveats relevant to my work

- Subtype calls are sample-conditional: purity, batch, and platform
  shift the centroids. My pipeline takes calls as-is — a limitation
  documented in the manuscript.
- "Normal-like" frequently reflects normal-tissue contamination —
  matches my observation that it's the weakest learnable class.
- Mutation→subtype is fundamentally many-to-one: mutations
  statistically associate with expression-defined states but don't
  determine them — which is *why* my classifier plateaus at macro-F1
  0.40. This paper explains the ceiling better than any ML paper could.

## Takeaway for portfolio work

Any claim "subtype X has property Y" implicitly depends on the
expression-defined identity of X. Testing mutation→subtype prediction
tests the *strength of the association*, not a labeling error.
