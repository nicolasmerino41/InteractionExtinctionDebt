# Single-site removal: support fragility

All ten datasets use the validated site-interaction incidence exported by analysis 42. Support k is the number of distinct sites recording an interaction, not observation frequency. No 50% removal scenario or null model is used.

- A: sum of 1/(k-1) - 1/k over interactions at the removed site with k > 1. This reciprocal-support index weights scarce support more heavily. It is not an extinction probability. Interactions lost completely are excluded and reported in B; thus A alone is not total damage.
- B: number of regional interactions lost immediately (k = 1).
- C: number of interactions newly left at one supporting site (k = 2). Existing singleton interactions elsewhere are not included.

All interactions follow the same weighting rule. Raw sums measure total site impact; fragility_per_interaction in the table separates average contribution from richness. Richness associations are descriptive Pearson/Spearman correlations, not proof that richness explains all variation. Constant outcomes have undefined correlations (NaN).

The output root contains three faceted figures; each dataset folder contains corresponding two-panel figures (richness order and actual richness scatter) and site-level tables. Rank spacing is uniform and does not represent differences in richness. Despite the requested folder name, these figures measure removal effects, not nestedness or containment.

Validation passed: toy examples, distinct incidences, support recalculated from site identities, and identities summing full loss, new singleton counts and fragility across sites.
