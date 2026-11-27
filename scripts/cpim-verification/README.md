# CPIM solver verification

The follow-up study verifies DGM and reruns the CPIM reference and dense grids. It includes a controlled serial benchmark of both models and the pre-refactor CPIM revision. [DGM verification](DGM.md)

## Reference and reproduction

The reference integrates total pressure and Sb₄ mole fraction in distance from the outlet using SciPy DOP853. A trial that reaches equilibrium early supplies a positive objective from the remaining tube length. A trial that reaches the inlet below equilibrium supplies the negative pressure deficit. Brent's method finds the zero of this continuous objective. A hydrogen-depletion event or integrator failure raises an error. Each final flux is integrated across the full tube to produce the pressure and composition profiles. The transport and boundary equations come from the pinned theory report. The two vapor-pressure correlations use the legacy coefficients with Kelvin inputs and Pa outputs. [Theory revision](https://github.com/alexandervolkov1/sb-distillation-theory/blob/93c4a9e09dae96427635c2e8130a058c982b6931/main.typ) · [Legacy correlations](../../src/sb_distillation/reference/legacy.clj)

The reference refinement used `(rtol, atol_p, xtol_N1)` of `(2e-11, 2e-12 Pa, 2e-14 mol/(m²·s))` and `(2e-13, 2e-14 Pa, 2e-16 mol/(m²·s))`. The absolute mole-fraction tolerance is `atol_p / p1_sat`. Over all 189 solutions the change was at most `7.83e-13` relative in Sb₄ flux and `1.29e-9` and `1.21e-9 Pa` in the two pressure profiles at 11 common stations. The tighter reference produces the EDN fixture with full double precision. [Reference source](reference.py) · [Fixture](../../test/fixtures/cpim-reference.edn)

Run the following commands from the repository root. The reported run used Python 3.14.6, SciPy 1.18.1, and NumPy 2.5.5. SciPy is needed to regenerate the reference or analyze convergence. The ordinary Clojure tests read the committed fixture.

```powershell
New-Item -ItemType Directory -Force target/cpim-verification | Out-Null
python -m pip install scipy==1.18.1 numpy==2.5.5
python scripts/cpim-verification/reference.py --fixture test/fixtures/cpim-reference.edn --study target/cpim-verification/reference.edn
clojure -M scripts/cpim-verification/study.clj target/cpim-verification/reference.edn target/cpim-verification/pressures.csv --full-grid
python scripts/cpim-verification/reference.py --pressures target/cpim-verification/pressures.csv
clojure -M:test
```

Omit `--full-grid` for the 189-point study. The current benchmark runs separately from the grid checks. The validation regressions run with `clojure -M:test`. [Solver tests](../../test/sb_distillation/cpim/solver_test.clj) · [Benchmark instructions](DGM.md#performance)

## Baseline at 700 °C, 10 mm, 200 ml/h

The tube length is 0.500 m. The flow reference is 293.15 K and 101325 Pa. The outlet uses 0.0368 Torr directly as hydrogen pressure, with zero trap drop. The mass is 0.500 kg. Both implementations use the same physical parameter values.

| Quantity | SciPy reference | Clojure defaults |
| --- | ---: | ---: |
| H₂ molar flow, mol/s | 2.309510939589e-6 | 2.309510939589e-6 |
| H₂ molar flux, mol/(m²·s) | 0.02940560657283 | 0.02940560657283 |
| Liquid saturation pressure, Pa | 64.24681994631 | 64.24681994631 |
| Sb₄ molar flux, mol/(m²·s) | 0.05748751772423 | 0.05748803510203 |
| Inlet Sb₄ partial pressure, Pa | 64.24681994631 | 64.24722190616 |
| Inlet H₂ partial pressure, Pa | 24.99996037695 | 24.99992207720 |
| Inlet total pressure, Pa | 89.24678032326 | 89.24714398335 |
| Sb₄ mass flow, kg/s | 2.199014376217e-6 | 2.199034166970e-6 |
| Distillation time, h | 63.15960931906 | 63.15904089853 |

| Source-side coordinate z, m | Reference p₁, Pa | Reference p₂, Pa | Reference x₁ | Reference x₂ |
| ---: | ---: | ---: | ---: | ---: |
| 0.00 | 64.246820 | 24.999960 | 0.719878 | 0.280122 |
| 0.10 | 56.972593 | 21.437776 | 0.726595 | 0.273405 |
| 0.20 | 48.786166 | 17.473611 | 0.736286 | 0.263714 |
| 0.30 | 39.214572 | 12.929415 | 0.752044 | 0.247956 |
| 0.40 | 27.098057 | 7.464976 | 0.784018 | 0.215982 |
| 0.45 | 18.746234 | 4.387586 | 0.810339 | 0.189661 |
| 0.50 | 0.000000 | 4.906263 | 0.000000 | 1.000000 |

The fixture includes 11 stations and four fixed-flux trials at zero, half, one, and eight times the molecular flux estimate. Tests compare partial pressures, mole fractions, and initial derivative signs and magnitudes. The zero-flux Clojure trial retains zero Sb₄ pressure throughout. A separate zero-H₂-feed case with positive outlet H₂ pressure also solves.

## Operating grid

The grid spans 680–720 °C by 5 °C, 6–12 mm by 1 mm, and 100, 200, and 300 ml/h of H₂. The corresponding gauge readings are 0.02010, 0.03680, and 0.05110 Torr at source and receiver temperatures of 695 °C and 510 °C. The 150 and 250 ml/h readings are outside this grid. Reference flow conditions and direct gauge interpretation are provisional. The liquid-Sb Antoine fit covers 1058–1373 K, above the entire 953.15–993.15 K operating interval.

All 189 cases passed with default options. The study checked all accepted pressure states and compositions, strict coordinate ordering, the exact outlet boundary, the inlet residual, and positive finite mass flow. It compared each case with the independent reference, including the geometric center and all corners. The largest relative error in flux and mass flow was `2.886e-5` at 680 °C, 10 mm, and 100 ml/h. Maximum inlet H₂ error was `1.099e-4 Pa` at 720 °C, 9 mm, and 300 ml/h. The largest absolute inlet Sb₄ residual was `9.984e-4 Pa`. Maximum errors across the 11-station profiles were `9.984e-4 Pa` for Sb₄ and `1.099e-4 Pa` for H₂. No doubling trial aborted the bracket search.

| Metric | Minimum | Maximum |
| --- | ---: | ---: |
| Sb₄ mass flow, kg/s | 3.228204e-7 | 6.364076e-6 |
| Distillation time, h | 21.8239 | 430.2358 |
| RK4 attempts in returned profile | 36 | 52 |
| RK4 attempts across each solve | 334 | 877 |
| Profile integrations per solve | 9 | 20 |

The slowest transfer occurred at 680 °C, 6 mm, and 100 ml/h. The fastest occurred at 720 °C, 12 mm, and 300 ml/h. An RK4 attempt includes accepted and rejected steps.

| Temperature, °C | Diameter, mm | H₂ feed, ml/h | Reference mass flow, kg/s | Relative Clojure mass-flow error | Absolute inlet H₂ error, Pa |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 680 | 6 | 100 | 3.22822143e-7 | 5.432e-6 | 1.080e-6 |
| 680 | 6 | 300 | 5.05384800e-7 | 9.516e-6 | 5.828e-5 |
| 680 | 12 | 100 | 2.19863297e-6 | 1.342e-6 | 4.217e-6 |
| 680 | 12 | 300 | 2.67432261e-6 | 9.559e-6 | 4.030e-5 |
| 700 | 9 | 200 | 1.62877013e-6 | 1.662e-5 | 8.687e-5 |
| 720 | 6 | 100 | 6.46341928e-7 | 7.916e-6 | 5.700e-5 |
| 720 | 6 | 300 | 9.40606882e-7 | 1.726e-6 | 8.935e-8 |
| 720 | 12 | 100 | 5.70662357e-6 | 6.578e-6 | 1.293e-5 |
| 720 | 12 | 300 | 6.36418781e-6 | 1.757e-5 | 9.890e-5 |

## Convergence and tolerances

The study compares 680 °C / 6 mm / 300 ml/h, 700 °C / 10 mm / 200 ml/h, and 720 °C / 12 mm / 100 ml/h. First it fixes the SciPy flux to isolate integration error. Each native Clojure node is compared directly with SciPy dense output at the same coordinate. At common stations, Clojure uses cubic Hermite interpolation with CPIM derivatives. Interpolating exact reference pressures and slopes on each Clojure mesh gives a separate interpolation-only estimate. For the default integration at fixed flux, the largest interpolation-only errors were `3.38e-5 Pa` and `2.29e-5 Pa`. The reported station errors can reflect cancellation between interpolation and integration errors.

| Fixed-flux RK4 tolerances (absolute Pa, relative) | Max inlet H₂ error, Pa | Max native-node errors p₁ / p₂, Pa | Max common-station errors p₁ / p₂, Pa |
| --- | ---: | ---: | ---: |
| 1e-8, 1e-6 | 4.46e-6 | 1.44e-5 / 1.44e-5 | 2.13e-5 / 2.24e-5 |
| 1e-9, 1e-6 | 4.35e-6 | 1.44e-5 / 1.44e-5 | 2.15e-5 / 2.25e-5 |
| 1e-10, 1e-6 | 4.34e-6 | 1.44e-5 / 1.44e-5 | 2.15e-5 / 2.26e-5 |
| 1e-8, 1e-7 | 1.32e-6 | 1.95e-6 / 1.70e-6 | 4.79e-6 / 3.41e-6 |
| 1e-8, 1e-8 | 1.59e-7 | 3.47e-7 / 2.56e-7 | 9.85e-7 / 5.97e-7 |
| 1e-9, 1e-7 | 1.13e-6 | 1.94e-6 / 1.70e-6 | 4.95e-6 / 3.51e-6 |
| 1e-10, 1e-8 | 1.92e-7 | 3.24e-7 / 2.54e-7 | 5.79e-7 / 5.30e-7 |
| 1e-12, 1e-10 | 4.68e-9 | 8.78e-9 / 5.99e-9 | 1.33e-8 / 1.43e-8 |

The next sweep solves for flux at each setting. The middle five rows hold RK4 at `1e-12 Pa` absolute and `1e-10` relative to isolate shooting error. Relative mass-flow error equals relative flux error for fixed area and molar mass.

| Run | Residual tolerance, Pa | Max relative flux and mass-flow error | Max inlet H₂ error, Pa | Max common-station errors p₁ / p₂, Pa |
| --- | ---: | ---: | ---: | ---: |
| All defaults | 1e-3 | 9.52e-6 | 5.83e-5 | 4.51e-4 / 6.25e-5 |
| Tight RK4 | 1e-3 | 9.52e-6 | 6.27e-5 | 4.41e-4 / 7.05e-5 |
| Tight RK4 | 1e-4 | 1.49e-6 | 9.84e-6 | 7.27e-5 / 1.10e-5 |
| Tight RK4 | 1e-5 | 1.58e-7 | 3.34e-7 | 7.98e-6 / 3.34e-7 |
| Tight RK4 | 1e-6 | 1.39e-8 | 7.98e-8 | 6.27e-7 / 9.80e-8 |
| Tight RK4 | 1e-9 | 1.27e-10 | 5.52e-9 | 1.50e-8 / 1.52e-8 |
| Candidate RK4 (1e-9 Pa, 1e-7) | 1e-4 | 1.49e-6 | 8.71e-6 | 7.25e-5 / 9.06e-6 |

The default shooting tolerance dominates the transfer-rate error. The accepted settings for preliminary mass-flow and time maps are the current defaults of `1e-8 Pa` absolute and `1e-6` relative for RK4 and `1e-3 Pa` absolute for shooting. These defaults yielded a maximum relative mass-flow error of about `3e-5` on the reference grid. The tighter candidate in the table was tested at the three convergence points.

## Limiting regimes

The molecular comparison uses the analytic linear partial-pressure profiles and `N1 = DK1 * p1_sat / (RT * L)`. At 973.15 K, the outlet is 0.2 Pa and the hydrogen volumetric feed scales with the cube of diameter to keep its limiting pressure rise fixed. As the diameter decreases from 10 to 1 to 0.1 µm, the relative CPIM flux difference from the molecular limit falls from `5.59e-3` to `5.58e-4` to `5.58e-5`. Maximum profile differences fall from `0.0898 / 0.1414 Pa` to `0.000900 / 0.00142 Pa` for Sb₄ and H₂.

The continuum comparison removes Knudsen diffusion resistance and molecular pressure conductance from the independent equations. The tube diameter stays 10 mm and the feed stays 200 ml/h. At outlet pressures of 100, 1000, and 10000 Pa the relative flux difference decreases from `5.26e-2` to `2.69e-2` to `2.78e-3`. Maximum pressure differences decrease from `4.22 / 9.76 Pa` to `0.00434 / 0.00446 Pa`. Across both series Clojure agrees with the full CPIM reference within `3.15e-9` relative in flux and `1.43e-6 Pa` in common-station profiles. These are synthetic limit tests of the full boundary-value problem.

## Findings from the initial verification

The table records behavior at the initial verification commit `ff2b67dc5727b008129a9e055986f67988b37023`. The examples use the baseline conditions unless overridden. The original optional defect suite had seven assertion failures and five unexpected exception types. The validation follow-up addresses each case and moves the regressions into the standard suite.

| Reproduction | Observed behavior | Proposed fix |
| --- | --- | --- |
| Set `:sb4-collision-diameter` to `-0.600e-9` | Preparation succeeds and changes `K_12` from about 20.42 to 166.88 Pa·m²/s | Validate transport parameters before arithmetic, with field-specific errors and positive domains for masses, diameters, constants and reference quantities. Allow a nonnegative Sutherland constant. |
| Set `:h2-reference-viscosity` to NaN or infinity | Preparation returns nonfinite viscosity and Wilke coefficients | Validate the supplied constants and require finite positive derived coefficients. |
| Set `:gas-const` to zero, omit it, or pass a string | A raw arithmetic, null-pointer, or class-cast exception escapes | Reject the input with `ExceptionInfo` identifying the parameter and value. |
| Set `:diameter` to `1e200` or `1e-200` m | The area overflows to infinity and the H₂ flux becomes zero, or area underflow causes a raw division exception | Validate representable positive area before dividing. |
| Set `:ref-volumetric-flow` to `Double/MAX_VALUE` | Throughput, molar flow, and molar flux are infinite | Validate derived values before constructing the problem. Allow exact zero flow. |
| Set `:temp` to 153 K | Saturation underflows to zero and the bracket rejects the zero residual at its lower endpoint | Reject saturation underflow at preparation without changing the correlation. If zero saturation is physically intended, define its boundary behavior separately. |
| Pass `Long/MAX_VALUE` for both gauge and trap pressures | Checked integer addition throws an arithmetic exception despite a finite double sum | Define the numerical type contract at the SI boundary. Normalize supported numbers before arithmetic and check the result. |
| Call `(valid-pressure-state? [-1e-500M 1.0])` | A negative decimal passes the sign check after rounding to negative zero | Check the original signs with type-aware predicates before verifying representability and the summed pressure. |

The physical-constant defects could occur at ordinary operating conditions with invalid caller-supplied parameters. The overflow, underflow, and mixed-type cases were outside the tested operating domain. All working-range cases used the default transport parameters. The provisional gauge convention, trap calibration, saturation extrapolation, and effective collision diameter remain scientific uncertainties.

No physically invalid doubling trial occurred on the 189-point reference grid. Recovery from a failed doubling trial and the report's relative shooting criterion are deferred. A candidate recovery for future study is to reduce the trial increment between the last admissible negative residual and the failed flux, keeping only successful trials as bracket endpoints.

## Validation follow-up

`prepare-problem` checks signs before conversion and accepts finite doubles that preserve whether an input is zero. It validates the derived area, reference thermal energy, pressures, hydrogen flow quantities, and transport coefficients. Zero hydrogen feed stays zero, and a positive feed that rounds to zero at any flow conversion is rejected. Gauge and trap pressures are converted before addition. Input errors identify the parameter and calculation errors identify the quantity. [Input contract](../../ARCHITECTURE.md#physical-solver)

`clojure -M:test` passed 83 tests and 1635 assertions. All 189 reference-grid solutions passed. The numerical summary matched the initial verification, including the maximum relative mass-flow error of `2.886e-5` and maximum station pressure errors of `9.984e-4 Pa` for Sb₄ and `1.099e-4 Pa` for H₂. The convergence CSV matched byte for byte. All 7503 dense-grid solutions passed the profile and boundary checks. The measured serial solve time was 126.11 s, with a mean of 16.81 ms per call.

## Performance

The initial study warmed the JVM with five serial passes over the 189-point grid. It timed only `solve` with counters around trial evaluations. The run used Windows, Temurin 25.0.4.1, and an AMD Ryzen 5 6600H. The full 7503-point grid was timed serially with default tolerances and no contour computation. Each solve owns its cache. The retained profile state count includes the final bracket endpoints and all subsequent midpoint profiles. Earlier negative trials are replaced during bracketing.

| Grid | Successful calls | Mean solve, ms | Max solve, ms | Mean / max integrations | Mean / max RK4 attempts | Max cached states | Serial solve time, s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Reference | 189 | 16.79 | 22.58 | 18.10 / 20 | 698.76 / 877 | 731 | 3.17 |
| Full | 7503 | 16.88 | 27.09 | 18.13 / 20 | 695.90 / 913 | 761 | 126.62 |

The mean from the 189-point timing predicts 126 seconds for 7503 calls. All 7503 measured solutions also passed the accepted-profile and boundary checks.

Profile retention grows with shooting evaluations and accepted ODE steps. For serial batches the live trial cache is per solve. Concurrent solves multiply that cost, and retaining all returned profiles adds memory proportional to the batch size. Keep scalar results for level maps if profiles are only needed for selected points. The study counted retained states and did not measure heap bytes.

## Verification scope

Mathematical checks cover the pinned CPIM pressure and composition equations, boundaries, units, and mass conversion. Numerical verification compares Clojure with a separate formulation and integrator, including refinement and both limits. Experimental validation of Sb transfer rates remains open. The gauge readings are operating inputs for this study.
