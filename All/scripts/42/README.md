# Analysis 42 — site nestedness, complementarity and site importance

Run `Rscript All/scripts/42/run42.R` from the repository root. Set `IED_PYTHON` to a Python executable with numpy if the default bundled runtime is unavailable. R uses the existing shared loader and ggplot2/dplyr/tidyr/jsonlite. Output: `All/outputs/42`.

## Questions and definitions

Each column is an initially realised regional interaction, each row a site. Containment is shared interactions divided by the poorer site's richness; its overall mean excludes equal-richness and empty-site pairs. Pairwise Sørensen dissimilarity is partitioned into Simpson turnover and nestedness-resultant dissimilarity following Baselga (2010), https://doi.org/10.1111/j.1466-8238.2009.00490.x. The latter is not a nestedness index. Turnover and nestedness-component summaries use all nonempty site pairs. Site pairs are descriptive observations, not independent inferential replicates.

## Null models

Margins preserves every row sum (site richness), every column sum (interaction support) and the regional link set. Opportunity preserves those same quantities and forbids cells lacking recorded pair co-occurrence. The first model is deliberately a structural baseline and may place interactions outside species ranges; the second is the biologically constrained comparison.

Four chains start from the observed matrix. Each update selects two distinct rows uniformly. Their shared links stay fixed; exclusive links eligible in both rows are pooled and reassigned uniformly while retaining row sizes. Links forbidden in the other row stay fixed. This conditional row-trade kernel is symmetric and has the uniform distribution as a stationary distribution within its reachable state space. It does NOT guarantee irreducibility when structural zeros constrain movement.

Each chain burns 100 sweeps and saves 250 draws at five sweeps between draws (one sweep = number of sites row-pair update attempts). Total 1,000 draws per dataset/model. Reported diagnostics: changed-trade frequency, relocation from observed incidence, scalar split R-hat and positive-autocorrelation ESS. These are diagnostics, not convergence proofs. Approximate tail probabilities are suppressed for invariant metrics or split R-hat >=1.05 or ESS <100. Global BH families are each metric across eligible dataset/model tests; site-level BH families are dataset × model across eligible sites. Markov draws are not independent permutations; probabilities remain exploratory. Invariant metrics are not evidence for a lack of ecological structure.

## Site importance

Exactly as the link endpoint in script 41: expected fraction of links lost when focal site s is forced into an m-site removal set minus the expectation for uniform m-site removal, with m approximately N/2. All other removed sites are uniform without replacement. This is total excess damage, not solely future vulnerability after immediate damage. Null-adjusted importance subtracts the null mean for that same site. Single-site links are shown separately. Support, richness and amount of loss are controlled. Null uncertainty and diagnostics are saved for each site.

## Secondary diagnostics and limitations

Individual folders include similarity heatmaps (no formally inferred groups), chain traces and exploratory richness-ordered/removal comparisons. They are computed consistently for all datasets to avoid selective execution, but only interpret removal differences where primary structure comparisons pass diagnostics and support a departure. Random-removal expected means are analytically identical under fixed supports; simulations show only finite-sample deviations. Variance is saved separately. Removal comparisons randomise ties and use 40 snapshots per null model, so do not treat their precision as equivalent to the 1,000-draw primary analysis.

Fixed margins do not remove sampling/detection bias. Opportunity constraints inherit the shared loader's occurrence definitions. No environmental drivers, geographic barriers, demographic extinction, or causal habitat-loss response are identified. Site similarity and rare-link patterns alone do not demonstrate ecological mechanisms.

Three combined PNG/PDF figures are written at the output root, with matching individual figures plus supplementary diagnostic figures in each dataset folder. CSV files preserve underlying observations, null comparisons, diagnostics and settings. Main figures are exploratory and should be read alongside the diagnostic flags.
