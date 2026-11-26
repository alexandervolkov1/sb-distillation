# Development roadmap

## Current state

- The CPIM mathematical core is complete and tested. It assembles constant and local coefficients and evaluates both partial-pressure ODE right-hand sides. The seven core functions and seven basic transport functions agree with the pinned report. [Implemented functions](ARCHITECTURE.md#software-architecture)
- The generic numerical core is complete and tested. It provides adaptive RK4 with step doubling and normalized local error control, plus bracketed bisection with exact residual comparisons across numeric types. Tests cover analytical solutions, RK4 convergence, failure diagnostics, and floating-point boundary cases. [Numerical core](ARCHITECTURE.md#numerical-core)
- `clojure -M:test` runs independent reference and property checks for the core and transport functions alongside the legacy regressions. Core checks include pure-component states and the continuum and free-molecular ODE limits. [Model tests](test/sb_distillation/cpim/model_test.clj)

## In progress

- Between the numerical core and physical solver stages.

## Next steps

1. Implement the physical CPIM solver's hydrogen flux and three pressure boundary conditions, with a consistent throughput convention for the trap.
2. Use adaptive RK4 to integrate trial profiles from the receiver at `z = L` toward the source at `z = 0`. Enforce physical bounds and retain diagnostics from failed trials.
3. Implement shooting with the inlet saturation residual. Bracket the nonnegative Sb₄ molar flux and use bisection to find the root. Check the boundary residuals and accepted profile. [Theory summary](ARCHITECTURE.md#mathematical-model)

## Planned

1. Check solver behavior in both limiting regimes and convergence of flux and pressure profiles under integration refinement. Compute Sb₄ mass flow from the accepted molar flux, tube area, and molar mass.
2. Use Clay for computational experiments and plots of transfer rates, pressures, and composition. Compare predictions with reference calculations and available experimental data.

## Deferred

- Model spatial deposition inside the receiving ampoule after the tube model is established. Define receiver geometry, thermal conditions, and the deposition model with the user.

## Open decisions

- Agree on CPIM-specific integration and shooting tolerances, and on how the physical solver reports failed brackets and rejected profiles.
- Confirm the inlet saturation-pressure correlation, its temperature range, and the trap conductance and throughput calibration.
- Justify the effective Sb₄ collision diameter and decide how to assess sensitivity of the transfer rate to uncertain transport inputs.
- Select operating conditions and measurements for experimental validation. Validation remains open.
