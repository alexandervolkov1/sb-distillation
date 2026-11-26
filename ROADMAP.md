# Development roadmap

## Current state

- [CPIM transport](src/sb_distillation/cpim/transport.clj) implements `knudsen-diffusivity` (`D_{K,i}`), `binary-diffusivity-factor` (`K_12`), `binary-diffusivity` (`D_12`), and `sb4-viscosity` (`μ_1`). Their expressions agree with the pinned theory. The default parameters already include the H₂ Sutherland constants.
- [The automated suite](test/sb_distillation_test.clj) exercises legacy geometry, vapor pressure, and transport regressions. All four CPIM functions are outside its coverage. The Clay analysis notebook builds as a scaffold.

## In progress

- Implement CPIM transport coefficients one block at a time, following the coefficient dependencies in the [theory summary](ARCHITECTURE.md#mathematical-model). Establish direct checks for the four existing functions and verify each subsequent block as it is implemented.

## Next steps

1. Implement H₂ dynamic viscosity using the documented Sutherland correlation. Implement Wilke interaction factors and the composition-dependent mixture viscosity. Check the reference temperature and pure-component limits.
2. Implement the remaining CPIM coefficients in dependency order. Prepare the quantities constant along the tube, then compute local pressure, composition, diffusivity and viscosity, the Bosanquet coefficients `𝒟_i`, and the pressure resistances `A_C`, `A_K`, and `A`.
3. Implement the two partial-pressure ODE right-hand sides. Check signs, units, and the continuum and free-molecular reductions before using the equations in a solver.

## Planned

1. Integrate from the known outlet pressures toward the source. Set up hydrogen input and the trap pressure boundary using a consistent throughput convention. Reject nonphysical trial profiles.
2. Implement shooting with the inlet saturation residual. Find a sign bracket for the nonnegative Sb₄ molar flux and solve for the root. Check all three boundary conditions and the physical bounds along the accepted profile.
3. Check both limiting regimes and numerical convergence of flux and pressure profiles under integration refinement. Compute Sb₄ mass flow from the accepted molar flux, tube area, and molar mass.
4. Use Clay for computational experiments and plots of transfer rates, pressures, and composition. Compare predictions with reference calculations and available experimental data.

## Deferred

- Model spatial deposition inside the receiving ampoule after the tube model is established. Define receiver geometry, thermal conditions, and the deposition model with the user.

## Open decisions

- Choose the CPIM integration and bracketed root-finding algorithms, libraries if needed, tolerances, and handling of failed brackets or rejected trials. The report gives examples. The legacy solver belongs to the DGM reference.
- Confirm the inlet saturation-pressure correlation and its applicable temperature range, the trap conductance and throughput calibration, and the operating range and data for validation.
- Justify the effective Sb₄ collision diameter and decide how to assess sensitivity of the transfer rate to uncertain transport inputs.
