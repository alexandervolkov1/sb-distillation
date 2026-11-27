# Development roadmap

## Current state

- The CPIM mathematical core is complete and tested. It assembles constant and local coefficients and evaluates both partial-pressure ODE right-hand sides. The seven core functions and seven basic transport functions agree with the pinned report. [Implemented functions](ARCHITECTURE.md#software-architecture)
- The generic numerical core is complete and tested. It provides adaptive RK4 with step doubling and normalized local error control, plus bracketed bisection with exact residual comparisons across numeric types. Tests cover analytical solutions, RK4 convergence, failure diagnostics, and floating-point boundary cases. [Numerical core](ARCHITECTURE.md#numerical-core)
- CPIM is numerically verified over the provisional 680–720 °C, 6–12 mm, 100–300 ml/h operating range. All 189 reference-grid cases passed the boundary and pressure checks. The full CPIM boundary-value solution approaches both transport limits in synthetic cases. [Verification study](scripts/cpim-verification/README.md)
- The physical solver validates required inputs before calculation and checks derived values for positivity and finite double representation. Regressions cover conversion and calculation overflow and underflow, mixed numeric types, and zero flow. The 189-point reference comparison and 7503-point dense grid passed after the validation fixes. [Input contract](ARCHITECTURE.md#physical-solver) · [Solver tests](test/sb_distillation/cpim/solver_test.clj)
- CPIM and DGM share transport, thermodynamics, the total-pressure gradient, boundary preparation, and the shooting solver. CPIM remains the default. The DGM implementation follows theory revision `e6a744f3a66f7a58dcb99a81f7a8655fa8bdacf6` and passes the baseline solve and model selection checks. The original independent CPIM verification used revision `93c4a9e09dae96427635c2e8130a058c982b6931`. [Accepted theory](https://github.com/alexandervolkov1/sb-distillation-theory/blob/e6a744f3a66f7a58dcb99a81f7a8655fa8bdacf6/main.typ#L1156-L1364) · [Model selection tests](test/sb_distillation/solver_test.clj) · [CPIM verification study](scripts/cpim-verification/README.md)
- Independent mathematical tests verify the DGM partial-pressure equations and their relation to CPIM. Both models pass the 189-point SciPy comparison and the 7503-point dense-grid checks. DGM approaches the molecular and continuum approximations in the full boundary-value problem. The maximum relative DGM mass-flow error against SciPy is `2.509e-5` in the reference grid. [DGM model tests](test/sb_distillation/dgm/model_test.clj) · [Numerical verification and benchmarks](scripts/cpim-verification/DGM.md)
- On the paired dense grids over the provisional operating range, DGM predicts 1.115–10.541% higher Sb₄ mass flow and 1.103–9.535% shorter distillation time than CPIM. The comparison includes the baseline, all corners, and temperature and diameter trends at each feed. [Model comparison](scripts/cpim-verification/COMPARISON.md)
- The pure Sb₄ vacuum case has a shared analytical solution for CPIM and DGM. It uses the same saturation pressure, Knudsen diffusivity and pure viscosity correlations as the binary models. The current theory revision is `cd291cd1ee7689b3cff53d99da2017ff9243e3d0`. [Pure-vapor theory](https://github.com/alexandervolkov1/sb-distillation-theory/blob/cd291cd1ee7689b3cff53d99da2017ff9243e3d0/main.typ#L1366-L1399) · [Tests](test/sb_distillation/physics/pure_vapor_test.clj)
- The Clay notebook displays eight preliminary time maps. The two zero-feed maps reuse the analytical results. The six positive-feed maps use `1e-8 Pa` absolute and `1e-6` relative RK4 tolerances and a `1e-3 Pa` absolute shooting tolerance. All binary map scalars match the verification grid. [Notebook](notebooks/sb_distillation/analysis.clj) · [Map verification](scripts/time-maps/README.md) · [Convergence and tolerances](scripts/cpim-verification/README.md#convergence-and-tolerances)

## Next steps

- Compare predictions with measured antimony transfer rates and uncertainty estimates when suitable data become available.

## Deferred

- Model spatial deposition inside the receiving ampoule after the tube model is established. Define receiver geometry, thermal conditions, and the deposition model with the user.
- Revisit a relative shooting tolerance and recovery from failed upper trials if future operating ranges require them.

## Open decisions

- Confirm the inlet saturation-pressure correlation, its temperature range, and the trap conductance and throughput calibration.
- Confirm the hydrogen-flow reference conditions and the interpretation of vacuum-gauge readings.
- Justify the effective Sb₄ collision diameter and decide how to assess sensitivity of the transfer rate to uncertain transport inputs.
- Select operating conditions and measurements for experimental validation. Validation remains open.
