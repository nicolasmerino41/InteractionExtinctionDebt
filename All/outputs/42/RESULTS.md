# Analysis 42: exploratory results

All ten datasets were analysed. The three main figures have individual versions in each dataset folder, alongside site-similarity heatmaps, null-chain traces and secondary removal curves. All figures are supplied as PNG and PDF. See METHODS.md for definitions, controls, diagnostics and interpretation limits.

## What the results support

1. **Strong site nestedness is not a general description of these recorded networks.** Mean containment ranges from 0.037 (Salix_Galpar) to 0.579 (Garraf_PP2). Containment is the fraction of a poorer site's interactions shared with a richer site, averaged over unequal-richness nonempty pairs. Turnover exceeds the nestedness-resultant component in every dataset. These are descriptive comparisons, not independent tests of every pair.
2. **No dataset shows significantly greater containment than the fixed-margin null after the exploratory multiple-test correction.** Garraf_PP2 has slightly *lower* containment than the unconstrained fixed-margin expectation (0.579 versus 0.589; approximate BH q = 0.038). Its opportunity-constrained comparison does not show a clear departure. Null agreement does not establish a causal explanation or prove absence of structure.
3. **Some site-importance departures occur under the unconstrained null.** At approximately 50% site removal, 19 sites pass the within-dataset/model BH threshold: 11 in Garraf_PP2, three each in Gottin_PP and Montseny, and two in Olot. None passes under the opportunity-constrained null among sites with adequate diagnostics. Many opportunity-constrained site statistics are invariant or poorly explored; those cases are unresolved, not evidence of ecological equivalence.
4. **Salix_Galpar requires special caution.** Its observed mean containment is 0.0372, mean turnover 0.9645 and nestedness-resultant component 0.0197. Its fixed-margin comparison is unexceptional. The opportunity-constrained comparison fails the present mixing criteria (containment split R-hat 1.194, estimated effective sample size 45), so no inferential probability is reported. More sampling and potentially a sampler supporting longer alternating cycles would be needed before interpreting this null. Passing scalar diagnostics in other datasets does not prove the structural-zero state space is fully reachable.

## Consequences for the ecological interpretation

The richness-importance relationship in analysis 41 does not demonstrate site nestedness. Rich sites can contain many interactions and cause substantial regional loss even when sites contain largely different interaction sets. This analysis measures that distinction directly.

The secondary removal curves describe rich-first, poor-first and random site removal. Under fixed interaction support, expected regional retention under uniform random removal is identical by construction; arrangement can instead affect targeted removal and variability. These exploratory curves should not be promoted as evidence of an additional mechanism without an informative primary comparison.

This analysis neither identifies geographic barriers nor distinguishes environmental effects from sampling and detection. Opportunity constraints use the occurrence definitions supplied by the existing loader; they do not establish independent measurements of species presence. Existing interaction-network nestedness findings are not equivalent to nestedness of this site-by-interaction matrix.

## Verification

All 102 recorded checks passed: known nested/disjoint examples, analytical loss probabilities against exhaustive small examples, observed opportunity validity, preservation of both margins and allowed cells at every saved null draw, and zero average focal-site excess importance. There are 1,000 saved draws from four chains per dataset and null model. Mixing flags and approximate probabilities are retained in the CSV tables. These checks verify implementation constraints, not ecological assumptions or convergence.
