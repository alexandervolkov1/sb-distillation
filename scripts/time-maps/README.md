# Distillation-time map checks

Run the notebook and the automated suite from the repository root.

```powershell
clojure -M:test
clojure -M:clay -r notebooks/sb_distillation/analysis.clj
```

The notebook writes `target/clay/sb_distillation.analysis.html`. The first render for a study and calculation revision runs 15,006 binary solves and 2,501 analytical pure Sb₄ calculations. Subsequent renders reuse the scalar cache in `target/distillation-maps/`. The key covers the operating conditions, grid, transport parameters, solver options, Clojure version, and sources used by the calculation. Changes to the notebook presentation and contour styling reuse the calculated data. To time a fresh calculation without deleting the notebook's cache, use a new directory under `target/`.

```powershell
clojure -M -e "(require '[sb-distillation.experiments.time-grid :as g]) (prn (select-keys (g/cached-grid g/default-study `"target/distillation-maps-timing`") [:calculation-seconds]))"
Measure-Command { clojure -M:clay -r notebooks/sb_distillation/analysis.clj | Out-Host }
```

## Data and contours

Generate or reuse `cpim-dense.csv` and `dgm-dense.csv` under `target/cpim-verification/` using the [numerical verification commands](../cpim-verification/DGM.md#reference-and-reproduction). Compare the map data and run the contour refinement from the repository root.

```powershell
clojure -M scripts/time-maps/check.clj
```

All 15,006 binary values of mass flow and elapsed time match the earlier dense-grid data exactly. The two pure-vapor panels have identical data. The map grid has 2,501 distinct coordinates per panel. The check also exercises cache key changes for operating inputs and solver tolerances.

The refinement checks three contours per map in cool and narrow, central, and hot and wide regions. It also checks the one-day contour wherever it crosses a map. At each of those 28 sections, it recalculates a local 0.01 mm diameter grid at a temperature between the main grid rows. It compares the interpolated level crossing with the main plot. The largest measured displacement is 0.0002192 mm, about 0.22% of one diameter grid interval. This supports the 0.1 mm and 1 °C spacing for these levels and display sizes. The plotted range is 0.8893107–38.1392226 days, so the 6, 12, and 18 hour levels fall outside every map.

## Browser check

With Microsoft Edge available, install Playwright in a Python environment and run the check against the generated HTML. Python is only used for this optional browser check and the independent reference generation.

```powershell
python -m pip install playwright
python -X utf8 scripts/time-maps/check_browser.py target/clay/sb_distillation.analysis.html target/clay/screenshots
```

The check opens the HTML in Edge and inspects the rendered charts, common logarithmic scale, physical legend labels, contour levels, selectable point tooltip, and embedded Old Standard font. It verifies that neither the 1440 px layout nor the 390 px layout overflows, each chart renders 2,501 cells, and the zero-feed panels have the same rendered colors. It checks for visible source code and JavaScript or Vega warnings and errors. The output directory contains the full notebook, an overview of all eight maps, a tooltip example, and the stacked mobile pair.

On Windows with an AMD Ryzen 5 6600H and Temurin 25.0.4.1, two fresh grid calculations took 277.15 s and 332.58 s. These measurements exclude cache-file serialization. The final `clojure -M:clay` run with the grid already cached took 10.61 s including JVM startup and HTML generation. Clay reported 1.28 s to evaluate the notebook itself.
