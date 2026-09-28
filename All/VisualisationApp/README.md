# Spatial interaction explorer

Open **index.html** in Edge, Chrome or Firefox. The app works locally with embedded data; no server or ChatGPT session is required. The standalone renderer may attempt to load optional interface assets from a CDN, but the matrix and controls have no external data or JavaScript dependencies.

## Explore

- Choose any of the ten datasets (Salix–Galpar opens first).
- Select a consumer or resource to view its partner set.
- Order by site richness / link support, or by shared-site similarity.
- Select a site with the menu or by clicking a matrix row.
- Use “Only links at selected site” to inspect its interactions, and “Remove selected site” to highlight immediate losses.
- Select an interaction to list every supporting site. The selected-site table lists all matching pairs, support counts and single-site-removal consequences.
- Reset returns to the full current dataset.

## Interpretation

Rows are sites; columns are regional consumer–resource pairs; filled cells are observations. Top bars show regional support; left bars show site richness in the current view. Filtering never changes the regional support denominator. Removal affects one selected site only, and changing that site replaces the scenario rather than accumulating removals.

Shared-site similarity uses average-linkage clustering with binary/Jaccard dissimilarity, computed on the complete dataset. It is exploratory ordering, not a significance test or geographic distance. Species-filtered views retain the global ordering.

“Loses every recorded partner” includes disappearance from interaction records. It is not demographic extinction or independently verified presence without interactions. Absence of an observed interaction is not proof of ecological absence. Data come from the existing shared project loader, with its sampling and occurrence limitations.

## Rebuild

From the repository root, run:

```text
Rscript All/VisualisationApp/build_data.R
python All/VisualisationApp/build_app.py
```

The first step uses the same R dependencies as `All/scripts/00_dataset_loaders_and_helpers_all.R`, plus jsonlite. It reads original data and writes `datasets.json` and `data_checks.csv` here. Python uses its standard library and the included renderer. The template and renderer are retained so the app can be rebuilt without a Codex installation.
