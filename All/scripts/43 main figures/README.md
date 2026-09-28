# Figures 2 and 3

Edit `main_figures.jl` for data summaries, colours, fonts, axes and panel layout. All analyses and figures use Julia/CairoMakie. `export_inputs.R` is only an input adapter to the established dataset loaders.

Run from the InteractionExtinctionDebt repository root:

```powershell
Rscript "All/scripts/43 main figures/export_inputs.R"
julia "All/scripts/43 main figures/main_figures.jl"
```

To edit styling and redraw from existing calculation tables:

```powershell
julia "All/scripts/43 main figures/main_figures.jl" . --plot-only
```

Figures are saved under `All/outputs/43 main figures` as PNG, PDF and SVG. Figure captions and analysis definitions are in `METHODS_AND_CAPTIONS.md`, outside the figures. Reproducible quantities are in `tables/`.

Dependencies: Julia packages CSV, DataFrames and CairoMakie plus standard libraries; the existing shared R loader's dependencies for input export. The main Julia script has an optional first argument specifying the repository root. Palette, dataset order, representative dataset and replicate count are defined at the top. Label offsets for Figure 2d are editable near its plot code.

The input adapter corrects the confirmed blank-name/V1 issue locally for this run. Earlier analyses are not overwritten. The previously identified one-link discrepancies in Garraf_HP and Olot remain documented in the output methods and should be reconciled before submission.

Figure 3 revision: edit figure3_revised.jl. The main script includes this file. To redraw Figure 3 only, preserving Figure 2, run: julia "All/scripts/43 main figures/figure3_revised.jl" . --plot-only. Omit --plot-only to recompute the group analysis. The revised caption is FIGURE3_REVISED_CAPTION.md; Figure 3 is exported as PNG, PDF and SVG.
