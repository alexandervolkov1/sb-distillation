# Working on sb-distillation

Read [ARCHITECTURE.md](ARCHITECTURE.md) for the model and code boundaries and [ROADMAP.md](ROADMAP.md) for the current stage. Inspect the relevant functions and tests before changing behavior. Use the Typst report linked from the architecture document as the primary source for formulas, notation, and units. Report conflicts among theory, implementation, and tests to the user for resolution.

## Scientific collaboration

- The user and ChatGPT develop the physics, mathematical model, numerical methods, and implementation together. Develop the core computational functions in small, sequential blocks. Explain the mathematical meaning, inputs, and place of each block in the model.
- The coding agent primarily tests, checks implementations, refactors, and carries out other explicitly assigned work. Implement a later stage of the mathematical model only when the user assigns it. Obtain the user's agreement before changing physical assumptions, equations, or numerical methods.
- Distinguish agreement with a formula, automated test coverage, and validation against experimental data when reporting results. Keep unresolved scientific choices visible in the roadmap.

## Implementation

- Compute internally in SI. Convert external quantities at the model boundary and make physical constants and model parameters explicit. Prefer small, pure functions and pass parameters explicitly.
- Write function documentation and comments in English. Each computational docstring gives the report symbol, the physical meaning, the result units, and the argument units where needed. Follow this example from [transport.clj](src/sb_distillation/cpim/transport.clj).

```clojure
(defn knudsen-diffusivity
  "D_{K,i} - Knudsen diffusivity in m^2/s.
  Temperature in K, diameter in m, molar mass [kg/mol]."
  [{:keys [gas-const]} temp diameter molar-mass]
  (* (/ diameter 3.0)
     (Math/sqrt
      (/ (* 8.0 gas-const temp)
         (* Math/PI molar-mass)))))
```

- Keep the physics and numerical algorithms in `src/`. Use Clay notebooks for computational experiments and visualization. Keep automated tests in `test/`, separate from the notebooks. Run `clojure -M:test` for the automated suite and report the scope of its coverage.

## Git and documentation

- Use a separate branch for significant changes. Start from the requested base, normally the current `main`, and preserve existing user work. Keep the diff focused and verify it before committing or publishing.
- Commit and publish branches when the user directs those actions. Merge into `main` only with explicit permission. Prefer a fast-forward merge when possible. Force-push requires explicit permission.
- Update `ARCHITECTURE.md` when the architecture changes and `ROADMAP.md` when a development stage is completed. Editing `AGENTS.md` requires explicit user permission.
