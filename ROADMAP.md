# Development roadmap

## Current state

- The CPIM mathematical core is complete and tested. It assembles constant and local coefficients and evaluates both partial-pressure ODE right-hand sides. The seven core functions and seven basic transport functions agree with the pinned report. [Implemented functions](ARCHITECTURE.md#software-architecture)
- `clojure -M:test` runs independent reference and property checks for the core and transport functions alongside the legacy regressions. Core checks include pure-component states and the continuum and free-molecular ODE limits. [Model tests](test/sb_distillation/cpim/model_test.clj)

## In progress

- Between the mathematical core and solver stages.

## Next steps

1. Choose the numerical integration method and agree on error control. Implement integration from the receiver at `z = L` toward the source at `z = 0`.
2. Implement hydrogen flux and the three pressure boundary conditions. Use a consistent throughput convention for the trap and enforce physical bounds on trial profiles.
3. Implement shooting with the inlet saturation residual. Find a sign bracket for the nonnegative Sb₄ molar flux, solve for the root, and check the boundary residuals and accepted profile. [Theory summary](ARCHITECTURE.md#mathematical-model)

## Planned

1. Check solver behavior in both limiting regimes and convergence of flux and pressure profiles under integration refinement. Compute Sb₄ mass flow from the accepted molar flux, tube area, and molar mass.
2. Use Clay for computational experiments and plots of transfer rates, pressures, and composition. Compare predictions with reference calculations and available experimental data.

## Deferred

- Model spatial deposition inside the receiving ampoule after the tube model is established. Define receiver geometry, thermal conditions, and the deposition model with the user.

## Open decisions

- Choose the CPIM integration and bracketed root-finding algorithms, libraries if needed, tolerances, and handling of failed brackets or rejected trials. The report gives examples. The legacy solver belongs to the DGM reference.
- Confirm the inlet saturation-pressure correlation and its applicable temperature range, the trap conductance and throughput calibration, and the operating range and data for validation.
- Justify the effective Sb₄ collision diameter and decide how to assess sensitivity of the transfer rate to uncertain transport inputs.
