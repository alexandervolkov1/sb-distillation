(ns sb-distillation.cpim.model
  (:require [sb-distillation.cpim.transport :as tr]))

(defn prepare-coefficients
  "Constant CPIM coefficients for an isothermal circular tube.
  Temperature in K, inner diameter in m.
  All returned quantities are in SI units."
  [{:keys [gas-const sb4-molar-mass h2-molar-mass]
    :as params}
   temp diameter]
  (let [m1 sb4-molar-mass
        m2 h2-molar-mass
        mu1 (tr/sb4-viscosity params temp)
        mu2 (tr/h2-viscosity params temp)]
    {:rt       (* gas-const temp)
     :diameter diameter
     :sqrt-m1  (Math/sqrt m1)
     :sqrt-m2  (Math/sqrt m2)
     :dk1      (tr/knudsen-diffusivity params temp diameter m1)
     :dk2      (tr/knudsen-diffusivity params temp diameter m2)
     :k12      (tr/binary-diffusivity-factor params temp)
     :mu1      mu1
     :mu2      mu2
     :phi12    (tr/wilke-interaction-factor mu1 mu2 m1 m2)
     :phi21    (tr/wilke-interaction-factor mu2 mu1 m2 m1)}))

(defn effective-diffusivity
  "D_eff,i - Effective species diffusivity in m^2/s,
  using Bosanquet's rule.
  Pressure in Pa, K_12 in Pa*m^2/s, D_{K,i} in m^2/s."
  [pressure k12 dk]
  (/ 1.0
     (+ (/ pressure k12)
        (/ 1.0 dk))))
