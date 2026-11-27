# Project architecture

The mathematical reference is the Typst report at revision `93c4a9e09dae96427635c2e8130a058c982b6931`. [sb-distillation-theory/main.typ][theory]

## Physical problem

The apparatus carries antimony vapor from a heated source ampoule through a heated quartz tube to a receiving ampoule in a continuous H₂ stream. The receiver connects through a trap to a vacuum pump, with a gauge downstream of the trap. The current goal is the tube transfer rate. A later model will address spatial deposition inside the receiver.

The model assumes steady, one-dimensional transport with cross-section averages in a straight circular tube at the uniform source temperature. The ideal gas contains only Sb₄ and H₂. Equilibrium with liquid antimony sets inlet saturation, and evaporation is assumed fast enough to sustain it. The receiver is an ideal Sb₄ sink. Both component flows are conserved in the tube, with wall deposition and chemical reactions excluded. The regulator fixes H₂ molar flow, the trap has constant effective conductance, and other line resistances are neglected. These are model assumptions. Experimental validation remains open. [Problem statement and assumptions][assumptions]

## Mathematical model

Young and Todd's Cylindrical Pore Interpolation Model (CPIM) links two capillary transport limits. At small Knudsen numbers, Maxwell–Stefan mutual diffusion describes relative species motion and the pressure gradient drives laminar viscous flow. Unequal tangential momentum transfer causes diffusion slip at the stationary wall, with equal diffuse reflection coefficients for both species. The continuum limit neglects viscous slip. At large Knudsen numbers, molecule–wall collisions dominate, and each partial-pressure gradient drives its species' Knudsen flux with fully diffuse reflection. [Limiting regimes][limits] · [Young and Todd (2005)][young-todd]

CPIM interpolates composition and total-pressure transport separately. Bosanquet reciprocal addition combines mutual diffusivity `D_12` and each Knudsen diffusivity `D_{K,i}` into the effective species coefficient `𝒟_i`. Reciprocal addition of the continuum and molecular resistances `A_C` and `A_K` gives `A`. The pressure gradient uses `A` and the species fluxes weighted by the square roots of their molar masses. The interpolation reproduces the two limits and gives an approximation in the transition region. [Interpolation and final equations][interpolation]

The independent coordinate `z` increases from source (`0`) to receiver (`L`). Two first-order ODEs govern `p_1(z), p_2(z)`, with index 1 for Sb₄ and 2 for H₂, `p = p_1 + p_2`, and `x_i = p_i/p`. The constant tube-frame fluxes `N_i` use mol/(m²·s), positive toward the receiver. Hydrogen flux `N_2` follows from the regulator's reference volumetric flow, pressure and temperature, and tube area `S`. The unknowns are `N_1` and the pressure profiles, including inlet H₂ pressure. [Boundary conditions and unknowns][boundaries]

The inlet imposes `p_1(0) = p_{1,sat}(T_s)`. The outlet has `p_1(L) = 0` and `p_2(L) = p_g + Q_2/C_trap`, using a common throughput convention for `Q_2` and `C_trap`. The shooting method integrates backward from the two known outlet pressures for each trial `N_1 ≥ 0` and seeks a sign-bracketed root of `F(N_1) = p_1(0; N_1) - p_{1,sat}(T_s)`. It checks physical pressure and composition bounds along the tube, the inlet residual against absolute and relative tolerances, and numerical convergence. Antimony mass flow is `M_1 N_1 S` in kg/s. The profiles also give total pressure and composition along the tube. [Shooting, checks, and mass flow][numerics]

The Knudsen diffusivities, hard-sphere Chapman–Enskog factor `K_12 = p D_12`, and pure-component viscosities are constant along the tube. Sb₄ viscosity uses hard spheres and H₂ viscosity uses Sutherland's correlation. Local pressure sets `D_12` and local composition sets Wilke mixture viscosity. The provisional 0.600 nm Sb₄ collision diameter affects Sb₄ viscosity and binary diffusion through inverse-square factors, with diffusion using the mean collision diameter. Hard-sphere and Sutherland approximations also contribute uncertainty. [Transport coefficients, uncertainties, and SI units][coefficients]

## Software architecture

- `sb-distillation-theory` owns the physical theory and mathematical definitions. `src/` holds the Clojure model and numerical code, `test/` holds automated checks, and `notebooks/` holds Clay experiments and visualization.
- `sb-distillation.cpim.transport` implements the basic coefficients in SI. `default-parameters` supplies the physical constants, species molar masses and collision diameters, and H₂ Sutherland parameters. Functions receive the parameters and operating inputs explicitly. [CPIM transport][transport]

| Implemented CPIM function | Quantity and dependencies |
| --- | --- |
| `knudsen-diffusivity` | `D_{K,i}` in m²/s from gas constant, temperature, tube inner diameter, and species molar mass. |
| `binary-diffusivity-factor` | `K_12` in Pa·m²/s from temperature, molecular reduced mass, and mean collision diameter. |
| `binary-diffusivity` | `D_12` in m²/s from `K_12` and local total pressure. |
| `sb4-viscosity` | `μ_1` in Pa·s from temperature and hard-sphere Sb₄ parameters. |
| `h2-viscosity` | `μ_2` in Pa·s from temperature and Sutherland parameters. |
| `wilke-interaction-factor` | Dimensionless `Φ_ij` from ordered viscosities and molar masses. |
| `mixture-viscosity` | `μ` in Pa·s from mole fractions, pure viscosities, and cross interaction factors. |

`sb-distillation.cpim.model` assembles the isothermal tube coefficients and evaluates the partial-pressure ODE right-hand side. [CPIM model][model]

- `molecular-resistance` computes `A_K`. `prepare-coefficients` prepares `RT`, diameter, square roots of molar masses, `D_{K,i}`, `K_12`, pure viscosities, Wilke factors, and `A_K` once per tube setup.
- `effective-diffusivity` computes `𝒟_i` (`D_eff,i` in code). `inverse-continuum-resistance` computes `1/A_C`, and `pressure-resistance` combines it with `A_K` to obtain `A`. `local-coefficients` uses the current partial pressures to compute `p`, mole fractions, mixture viscosity, both effective diffusivities, `1/A_C`, and `A`. It requires positive total pressure.
- `pressure-derivatives` returns `[dp_1/dz dp_2/dz]` for the given state and fluxes. The derivatives follow increasing `z` from source to receiver.

### Numerical core

The numerical namespaces operate on caller-supplied functions. The physical solver uses them to integrate the CPIM boundary-value problem and find the inlet saturation residual's root. [Development roadmap](ROADMAP.md)

**`sb-distillation.numerics.ode`** integrates vector initial-value problems in either direction using adaptive classical RK4. The caller supplies `f(z, state)` and can supply a `valid-state?` predicate. States are checked before each RHS evaluation and at the end of each RK4 step. Derivative vectors must match the state dimension. [ODE source][ode] · [ODE tests][ode-tests]

- `rk4-step` advances one signed step. `rk4-step-doubling` compares one full step with two half steps and returns the fine state as `:state` and the signed Richardson correction `(fine - coarse)/15` as `:error`.
- `normalized-error` takes the largest componentwise ratio of the absolute correction to `abs-tol + rel-tol * max(|initial|, |fine|)`. `next-step-size` scales the signed step by a bounded factor based on `safety * E^(-1/5)`, using the maximum growth factor at zero error.
- `attempt-step` accepts the fine state at finite normalized error `E <= 1`. Its map reports `:accepted` with `:state` and `:error`, or `:rejected` with `:reason`. Excessive error includes the error estimate. Invalid states and nonfinite derivatives include stage data in `:details`. Nonfinite error also rejects the trial. Invalid derivative shape and unexpected RHS exceptions propagate.
- `integrate` applies `default-options` and caller overrides for tolerances, step control, and attempt limits. It returns `:status`, `:z`, `:state`, `:attempts`, and an accepted `:profile` of coordinate and state maps. Attempts include rejections. Success reaches the supplied end coordinate. Failure preserves the accepted profile and reports `:max-attempts`, `:step-too-small`, `:coordinate-stagnation`, or `:step-size-unrepresentable` in `:reason`. A step-too-small result also includes `:last-rejection`. Invalid arguments, options, or initial states throw.

`integrate` uses the displacement between representable double coordinates for each trial. A final step may be shorter than `:min-step`, and stage coordinates can coincide near the resolution limit. States and derivatives must be finite in double precision. Wider arithmetic protects the RK4 weighted increment and the error scale at extreme magnitudes. Scaling protects the Richardson subtraction when its intermediate difference overflows.

**`sb-distillation.numerics.root`** provides `bisect` for a continuous scalar residual on an ordered finite interval. An endpoint can satisfy the absolute tolerance directly. Otherwise, the endpoints must have opposite signs. Each new midpoint counts as one iteration, and endpoint values are reused. The result contains `:status`, `:x`, `:fx`, `:bracket`, and `:iterations`. A `:success` satisfies `abs(f(x)) <= f-tol`. A `:failed` result includes `:reason`. The reasons `:not-bracketed`, `:max-iterations`, and `:midpoint-stagnation` return the current endpoint with smaller absolute residual. `:non-finite-value` identifies the offending evaluation. Invalid arguments or nonnumeric residuals throw, and user exceptions propagate. [Root source][root] · [Bisection tests][root-tests]

Bisection coordinates use double precision and can exhaust the available midpoints. The residual and tolerance retain their Clojure numeric types. Comparison of two floating-point values uses the floating-point path. Mixed comparisons use exact cross products with binary floating-point operands converted to their exact values. This preserves the residual sign and tolerance comparison even outside the double range.

### Physical solver

- `sb-distillation.cpim.thermodynamics` provides the liquid and solid antimony vapor-pressure correlations in Pa with Kelvin inputs. The liquid correlation supplies the inlet boundary condition. The working temperatures use its extrapolation below the 1058–1373 K fit interval. [Thermodynamics source][thermodynamics]
- `sb-distillation.cpim.solver` prepares the SI boundary-value problem from the regulator reference conditions, tube geometry, trap pressure drop, and transport parameters. `prepare-problem` accepts Clojure numbers, checks their original signs, and normalizes them to finite IEEE 754 doubles. Representable subnormal values are allowed. A nonzero input that converts to zero is rejected. Required inputs are strictly positive except for the Sutherland constant, volumetric H₂ feed, gauge pressure, and trap pressure drop, which can be zero. The converted source temperature must exceed 152.352 K, the Antoine correlation's singularity. [Solver source][solver]
- Preparation requires positive finite area, reference `RT`, saturation and outlet pressures, and constant transport coefficients. Hydrogen throughput, molar flow, and molar flux must remain positive for positive feed. All three are zero for zero feed, with positive outlet pressure still required. `ExceptionInfo` identifies invalid inputs by `:parameter` and failed calculations by `:quantity`. The double range applies to this physical boundary. The generic ODE and root modules have their own input contracts.
- The solver integrates trial profiles from the receiver to the source and enforces finite, physically admissible pressure states at RK4 stages. The bracket search doubles the molecular flux estimate and stops with diagnostics if a trial fails. Each solve reuses the bracket endpoint profiles and caches subsequent trials for that call. A successful solve returns the selected profile and Sb₄ mass flow in kg/s. Failure identifies the bracketing or shooting stage and retains the corresponding diagnostics.
- Preliminary level maps use the default RK4 tolerances of `1e-8 Pa` absolute and `1e-6` relative and the absolute shooting tolerance of `1e-3 Pa`. A relative shooting tolerance and recovery from a failed trial are deferred. [Verification study][verification]

### Reference calculations and verification

- Legacy references preserve ampoule geometry, vapor pressures, conductance and Knudsen corrections, viscosity, diffusion, and pure and binary Dusty Gas Model (DGM) rates. Binary DGM uses RK4 and bisection. Some legacy APIs use °C, Torr, or g/h. The early empirical estimate sets receiver pressure by solid-Sb saturation. DGM and CPIM theory use zero outlet Sb₄ pressure. [Legacy code][legacy] · [Legacy units and assumptions][legacy-docs]
- [Plot helpers][plot] prepare a Vega-Lite theme and grid layers and extract contour points by interpolating the first crossing in log time. [Clay analysis notebook][notebook] is currently a namespace and heading scaffold for future experiments.
- Automated tests cover all seven basic transport functions using independent references, scaling laws, and viscosity limits and symmetries. The model tests cover all seven core functions, including coefficient maps, pure-component states, flux linearity, pressure-gradient consistency, and both ODE limits at finite positive pressures. [Transport tests][cpim-tests] · [Model tests][model-tests]
- Solver and thermodynamics tests cover boundary preparation, both pressure correlations, trial and full profiles against an independent SciPy reference, and numerical failure propagation. The separate verification study compares all 189 operating points against DOP853 and Brent, studies convergence and both limiting regimes, and measures serial JVM performance. [Solver tests][solver-tests] · [Thermodynamics tests][thermodynamics-tests] · [Verification study][verification]
- Legacy tests cover geometry, saturation pressures, and four historical transport outputs. The liquid-Sb test at 700 °C checks the correlation below its 1058–1373 K Antoine fit interval. [Legacy tests][tests] · [NIST fit data][nist]

Run from the project root using the aliases in [deps.edn](deps.edn).

```powershell
clojure -M:test
clojure -M:clay -r notebooks/sb_distillation/analysis.clj
```

Clay writes `target/clay/sb_distillation.analysis.html`. [Clay configuration](clay.edn)

[theory]: https://github.com/alexandervolkov1/sb-distillation-theory/blob/93c4a9e09dae96427635c2e8130a058c982b6931/main.typ
[assumptions]: https://github.com/alexandervolkov1/sb-distillation-theory/blob/93c4a9e09dae96427635c2e8130a058c982b6931/main.typ#L32-L88
[boundaries]: https://github.com/alexandervolkov1/sb-distillation-theory/blob/93c4a9e09dae96427635c2e8130a058c982b6931/main.typ#L103-L208
[limits]: https://github.com/alexandervolkov1/sb-distillation-theory/blob/93c4a9e09dae96427635c2e8130a058c982b6931/main.typ#L210-L554
[interpolation]: https://github.com/alexandervolkov1/sb-distillation-theory/blob/93c4a9e09dae96427635c2e8130a058c982b6931/main.typ#L556-L752
[numerics]: https://github.com/alexandervolkov1/sb-distillation-theory/blob/93c4a9e09dae96427635c2e8130a058c982b6931/main.typ#L754-L937
[coefficients]: https://github.com/alexandervolkov1/sb-distillation-theory/blob/93c4a9e09dae96427635c2e8130a058c982b6931/main.typ#L940-L1255
[young-todd]: https://www.sciencedirect.com/science/article/pii/S0017931005004904
[transport]: src/sb_distillation/cpim/transport.clj
[model]: src/sb_distillation/cpim/model.clj
[solver]: src/sb_distillation/cpim/solver.clj
[thermodynamics]: src/sb_distillation/cpim/thermodynamics.clj
[ode]: src/sb_distillation/numerics/ode.clj
[root]: src/sb_distillation/numerics/root.clj
[legacy]: src/sb_distillation/reference/legacy.clj
[legacy-docs]: https://github.com/alexandervolkov1/sb-distillation/blob/4784a763ff7d0c385fdeb27f44311637af7b2950/README.md#справочные-расчёты
[plot]: src/sb_distillation/reference/plot.clj
[notebook]: notebooks/sb_distillation/analysis.clj
[cpim-tests]: test/sb_distillation/cpim/transport_test.clj
[model-tests]: test/sb_distillation/cpim/model_test.clj
[solver-tests]: test/sb_distillation/cpim/solver_test.clj
[thermodynamics-tests]: test/sb_distillation/cpim/thermodynamics_test.clj
[verification]: scripts/cpim-verification/README.md
[ode-tests]: test/sb_distillation/numerics/ode_test.clj
[root-tests]: test/sb_distillation/numerics/root_test.clj
[tests]: test/sb_distillation_test.clj
[nist]: https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440360&Mask=4
