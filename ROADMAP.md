# Development roadmap

## Current state

- The basic CPIM transport coefficient stage is complete. All seven functions for Knudsen and binary diffusion, pure viscosities, Wilke interaction factors, and mixture viscosity agree with the formulas and SI dimensions in the pinned report. [Implemented functions](ARCHITECTURE.md#software-architecture)
- `clojure -M:test` runs direct checks of all seven functions against independent numerical values and mathematical properties, alongside the legacy regressions. The Clay analysis notebook builds as a scaffold. [CPIM tests](test/sb_distillation/cpim/transport_test.clj)

## In progress

- Between coefficient stages.

## Next steps

1. Prepare the quantities constant along the isothermal tube for given geometry and parameters, including `D_{K,i}`, `K_12`, pure viscosities, Wilke interaction factors, and `A_K`.
2. Compute local pressure, composition, diffusivity, and mixture viscosity. Implement the Bosanquet coefficients `𝒟_i` and pressure resistances `A_C` and `A`. Check the continuum and free-molecular limits. [Theory summary](ARCHITECTURE.md#mathematical-model)

## Planned

1. Implement the two partial-pressure ODE right-hand sides and check their signs, units, and limiting forms.
2. Integrate from the known outlet pressures toward the source. Set up hydrogen input and the trap pressure boundary using a consistent throughput convention. Reject nonphysical trial profiles.
3. Implement shooting with the inlet saturation residual. Find a sign bracket for the nonnegative Sb₄ molar flux and solve for the root. Check all three boundary conditions and the physical bounds along the accepted profile.
4. Check both limiting regimes and numerical convergence of flux and pressure profiles under integration refinement. Compute Sb₄ mass flow from the accepted molar flux, tube area, and molar mass.
5. Use Clay for computational experiments and plots of transfer rates, pressures, and composition. Compare predictions with reference calculations and available experimental data.

## Deferred

- Model spatial deposition inside the receiving ampoule after the tube model is established. Define receiver geometry, thermal conditions, and the deposition model with the user.

## Open decisions

- Choose the CPIM integration and bracketed root-finding algorithms, libraries if needed, tolerances, and handling of failed brackets or rejected trials. The report gives examples. The legacy solver belongs to the DGM reference.
- Confirm the inlet saturation-pressure correlation and its applicable temperature range, the trap conductance and throughput calibration, and the operating range and data for validation.
- Justify the effective Sb₄ collision diameter and decide how to assess sensitivity of the transfer rate to uncertain transport inputs.
