(ns sb-distillation.cpim.model
  (:require [sb-distillation.physics.common :as physics]))

(def molecular-resistance physics/molecular-resistance)
(def prepare-coefficients physics/prepare-coefficients)
(def inverse-continuum-resistance physics/inverse-continuum-resistance)
(def pressure-resistance physics/pressure-resistance)

(defn effective-diffusivity
  "D_eff,i - Effective species diffusivity in m^2/s,
  using Bosanquet's rule.
  Pressure in Pa, K_12 in Pa*m^2/s, D_{K,i} in m^2/s."
  [pressure k12 dk]
  (/ 1.0
     (+ (/ pressure k12)
        (/ 1.0 dk))))

(defn local-coefficients
  "Local CPIM coefficients at partial pressures p_1 and p_2.
  Pressures in Pa, all returned quantities in SI units.
  Partial pressures must be nonnegative,
  total pressure must be positive."
  [{:keys [dk1 dk2 k12] :as coefficients} p1 p2]
  (let [{:keys [p] :as shared} (physics/local-coefficients coefficients p1 p2)]
    (assoc shared
           :d-eff1 (effective-diffusivity p k12 dk1)
           :d-eff2 (effective-diffusivity p k12 dk2))))

(defn pressure-derivatives
  "dp_1/dz, dp_2/dz - CPIM partial pressure derivatives in Pa/m.
  Partial pressures in Pa, molar fluxes in mol/(m^2*s)."
  [{:keys [rt] :as coefficients} p1 p2 n1 n2]
  (let [{:keys [x1 x2 d-eff1 d-eff2 a]}
        (local-coefficients coefficients p1 p2)

        h (- (/ (* x1 n2) d-eff2)
             (/ (* x2 n1) d-eff1))

        g (physics/total-pressure-gradient coefficients a n1 n2)
        rt-h (* rt h)]
    [(+ (* x1 g) rt-h)
     (- (* x2 g) rt-h)]))
