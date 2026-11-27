(ns sb-distillation.physics.pure-vapor
  (:require [sb-distillation.physics.common :as physics]
            [sb-distillation.physics.thermodynamics :as thermo]
            [sb-distillation.physics.transport :as transport]))

(defn molar-flux
  "N_1 - Steady pure Sb4 molar flux into vacuum in mol/(m^2*s).
  The inlet is saturated. Temperature in K, inner diameter and length in m."
  [{:keys [gas-const sb4-molar-mass] :as params} temp diameter length]
  (let [p-sat (thermo/sb-liquid-vapor-pressure temp)
        dk (transport/knudsen-diffusivity params temp diameter sb4-molar-mass)
        mu (transport/sb4-viscosity params temp)
        b0 (physics/viscous-permeability diameter)]
    (/ (+ (* dk p-sat)
          (/ (* b0 p-sat p-sat) (* 2.0 mu)))
       (* gas-const temp length))))
