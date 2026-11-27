# Development roadmap

## Current state

- The CPIM mathematical core is complete and tested. It assembles constant and local coefficients and evaluates both partial-pressure ODE right-hand sides. The seven core functions and seven basic transport functions agree with the pinned report. [Implemented functions](ARCHITECTURE.md#software-architecture)
- The generic numerical core is complete and tested. It provides adaptive RK4 with step doubling and normalized local error control, plus bracketed bisection with exact residual comparisons across numeric types. Tests cover analytical solutions, RK4 convergence, failure diagnostics, and floating-point boundary cases. [Numerical core](ARCHITECTURE.md#numerical-core)
- The physical solver is implemented and numerically verified over the provisional 680–720 °C, 6–12 mm, 100–300 ml/h operating range. All 189 reference-grid cases passed the boundary and pressure checks. The full boundary-value solution approaches both transport limits in synthetic cases. [Verification study](scripts/cpim-verification/README.md)
- The physical solver validates required inputs before calculation and checks derived values for positivity and finite double representation. Regressions cover conversion and calculation overflow and underflow, mixed numeric types, and zero flow. The 189-point reference comparison and 7503-point dense grid passed after the validation fixes. [Input contract](ARCHITECTURE.md#physical-solver) · [Solver tests](test/sb_distillation/cpim/solver_test.clj)
- Preliminary level maps will use the current defaults of `1e-8 Pa` absolute and `1e-6` relative for RK4 and `1e-3 Pa` absolute for the shooting residual. [Convergence and tolerances](scripts/cpim-verification/README.md#convergence-and-tolerances)

## Next steps

1. Use Clay for preliminary level maps with the accepted solver tolerances.
2. Compare predictions with measured transfer rates and uncertainty estimates when suitable data become available.

## Deferred

- Model spatial deposition inside the receiving ampoule after the tube model is established. Define receiver geometry, thermal conditions, and the deposition model with the user.
- Revisit a relative shooting tolerance and recovery from failed upper trials if future operating ranges require them.

## Open decisions

- Confirm the inlet saturation-pressure correlation, its temperature range, and the trap conductance and throughput calibration.
- Confirm the hydrogen-flow reference conditions and the interpretation of vacuum-gauge readings.
- Justify the effective Sb₄ collision diameter and decide how to assess sensitivity of the transfer rate to uncertain transport inputs.
- Select operating conditions and measurements for experimental validation. Validation remains open.
