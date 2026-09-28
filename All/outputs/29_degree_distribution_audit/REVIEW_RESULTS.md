# Nestedness, Figure 31 and degree-distribution audit

## Overall site nestedness

The Julia implementation reports one site NODF per complete site-by-interaction matrix. All ten results were independently checked against vegan::nestednodf (both site-row and full matrix scores agree to 1e-8). Main site scores on the 0-100 scale: Salix_Galpar 3.44; Garraf_HP 14.98; Gottin_PP 24.84; Garraf_PP 27.30; Montseny 27.43; Gottin_HP 33.56; Nahuel 34.95; Quercus 35.28; Olot 35.86; Garraf_PP2 57.60. These are descriptive scores, not significance tests. Whole-matrix NODF combining sites and interaction columns is also in the CSV; it answers a broader question and should not replace the site-row score for this question.

NODF is an aggregate of ordered overlaps. A whole-dataset score still uses pairwise overlap internally. Equal-richness pairs contribute zero under standard NODF, explaining the difference from the earlier mean containment restricted to unequal-richness pairs. The earlier nestedness-resultant dissimilarity component was not a general nestedness score.

## Figure 31 update

Original entry point: All/scripts/31_figure3B_link_states_by_dataset.R; prior version backed up beside it with suffix .before_annotations.R. The wrapper exports fresh pair supports from the existing loader and invokes 31_figure3B_render.jl for exact analysis and Makie plotting. The requested PNG and PDF are replaced in their original output folder.

The original removal levels and uniform random removal design are retained. Means are now exact hypergeometric expectations rather than Monte Carlo estimates over 500 subsets. For any original interaction, the pink-state probability equals the probability that all interaction-support sites disappear minus the probability that all co-occurrence sites disappear. State fractions sum to one. Pink area is integrated by trapezoids across the original six plotted levels, from 0 to 0.8 removal.

Labels show p_emp = total realised pair-site occurrences / total co-occurring pair-site opportunities, including pairs never observed interacting, and f = regional realised pairs / regional co-occurring pairs. Neither is the fitted homogeneous-model p. The table also supplies unweighted pair-average conversion for transparency.

Across ten datasets, pink area versus p_emp has Pearson r=-0.797883 and Spearman rho=-0.636364. Versus f: r=-0.650475, rho=-0.442424. These are descriptive associations, not independent mechanistic evidence. The strongest association is for conversion averaged over realised pairs only, a different conditional denominator; it must not be confused with all opportunities.

## Script 29 and Galiana

The original local Galiana script data-analysis_occurrence.R (lines 177-218) computes a reverse cumulative frequency table over observed degree values. Script 29's cdf_from_degrees computes mean(degree >= x) at every integer x. Both are complementary cumulative distributions of distinct regional partners, with consumer/resource guilds separate and both axes logarithmic. At zero removal the all-original and active-only denominators coincide. Script 29 does not fit the distribution families shown as dashed fitted curves in the paper.

The numerical CCDF values coincide at observed thresholds when degree data coincide. Connecting only observed degree values bridges missing degrees with sloping segments; inserting all integer thresholds reveals horizontal runs and drops. Script 29 additionally shares y ranges within guild rows through facet_grid; Galiana's separately drawn panels need not share those limits. These display choices affect appearance without changing the underlying distribution.

The audit compared current degree counts with every original network_real_pred/prey CSV stored in the repository. Seven datasets match exactly for both guilds. Exceptions:

- Gottin_HP: 108 original versus 113 current links; current has an extra consumer V1 with five partners and two extra resource nodes.
- Garraf_HP: 90 original versus 89 current links; degree differences affect Sapyga_quinquepunctata and Ancistrocerus longispinosus.
- Olot: 93 original versus 92 current links; differences affect Trichr..Cyanea and Passaloecus spp (suma2).

All zero-removal script-29 exports reproduce current incidence-based degree distributions to numerical precision. Thus script 29's CCDF calculation itself agrees with its inputs.

### Confirmed loader issue and scope

Gottin_HP V1 was traced to an unnamed consumer column produced by frame2webs. It contains positive entries at sites 7,18,20,29. clean_web_matrix converts the matrix to a data frame before testing for blank column names; that conversion renames the blank column V1, so it evades the intended blank-name filter. This is not evidence for a named consumer species. The shared loader was not changed during this audit, since doing so changes the inputs to all prior analyses. Gottin_HP results in this update and prior shared-loader analyses therefore remain provisional pending that coordinated correction. The two one-link discrepancies in Garraf_HP and Olot are established, but their preprocessing origin has not been attributed.

Reference: Galiana et al., Nature Ecology & Evolution 8, 209-217 (2024), https://doi.org/10.1038/s41559-023-02254-y . NODF: Almeida-Neto et al. (2008), https://doi.org/10.1111/j.0030-1299.2008.16644.x .
