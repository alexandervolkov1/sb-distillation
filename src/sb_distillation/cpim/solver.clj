(ns sb-distillation.cpim.solver
  "CPIM entry points for the shared tube boundary-value solver."
  (:require [sb-distillation.cpim.model :as model]
            [sb-distillation.solver :as solver]))

(def valid-pressure-state? solver/valid-pressure-state?)
(def cross-sectional-area solver/cross-sectional-area)
(def hydrogen-throughput solver/hydrogen-throughput)
(def hydrogen-molar-flow solver/hydrogen-molar-flow)
(def hydrogen-molar-flux solver/hydrogen-molar-flux)
(def trap-pressure-drop solver/trap-pressure-drop)
(def outlet-hydrogen-pressure solver/outlet-hydrogen-pressure)
(def prepare-problem solver/prepare-problem)
(def molecular-flux-estimate solver/molecular-flux-estimate)
(def default-solver-options solver/default-solver-options)
(def sb-mass-flow solver/sb-mass-flow)
(def distillation-time solver/distillation-time)

(defn integrate-trial
  "Integrate a CPIM trial from z = L to z = 0.
  N_1 in mol/(m^2*s), coordinates in m, pressures in Pa."
  ([problem n1]
   (solver/integrate-trial model/pressure-derivatives problem n1))
  ([problem n1 options]
   (solver/integrate-trial model/pressure-derivatives problem n1 options)))

(defn evaluate-trial
  "Integrate the CPIM trial and return F(N_1) in Pa on success.
  N_1 in mol/(m^2*s)."
  ([problem n1]
   (solver/evaluate-trial model/pressure-derivatives problem n1))
  ([problem n1 options]
   (solver/evaluate-trial model/pressure-derivatives problem n1 options)))

(defn find-shooting-bracket
  "Find a CPIM bracket for F(N_1), with fluxes in mol/(m^2*s)
  and residuals in Pa. max-trials includes the zero-flux trial."
  ([problem]
   (solver/find-shooting-bracket model/pressure-derivatives problem))
  ([problem ode-options max-trials]
   (solver/find-shooting-bracket model/pressure-derivatives problem ode-options max-trials)))

(defn solve
  "Solve the CPIM boundary-value problem by shooting.
  N_1 in mol/(m^2*s), pressures and residual tolerance in Pa."
  ([problem]
   (solver/solve-with model/pressure-derivatives problem))
  ([problem options]
   (solver/solve-with model/pressure-derivatives problem options)))
