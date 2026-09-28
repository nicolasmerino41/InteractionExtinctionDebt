# Single-site removal: support fragility

All ten datasets use the validated site-interaction incidence exported by analysis 42. Support k is the number of distinct sites recording an interaction, not observation frequency. No 50% removal scenario or null model is used.

- A: sum of 1/k - 1/(k+1) over interactions at the removed site. This reciprocal-support index weights scarce support more heavily. Interactions with k = 1 are included: their contribution is 1 - 1/2 = 0.5, so interaction extinctions remain represented in the fragility index.
- B: number of regional interactions lost immediately (k = 1).
- C: number of interactions newly left at one supporting site (k = 2). Existing singleton interactions elsewhere are not included.

All interactions follow the same weighting rule. Raw sums measure total site impact; fragility_per_interaction in the table separates average contribution from richness. Richness associations are descriptive Pearson/Spearman correlations, not proof that richness explains all variation. Constant outcomes have undefined correlations (NaN).

The output root contains three faceted scatter figures; each dataset folder contains corresponding individual figures and site-level tables. Actual richness is used on x, with no jitter or aggregation; exact duplicate points overlap. Each scatter plot includes an ordinary least-squares linear regression line as a descriptive summary.

These figures measure removal effects, not nestedness or containment.

Validation passed: toy examples, distinct incidences, support recalculated from site identities, and identities summing full loss, new singleton counts and fragility across sites.
