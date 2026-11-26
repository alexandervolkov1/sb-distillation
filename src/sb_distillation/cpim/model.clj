(ns sb-distillation.cpim.model
  (:require [sb-distillation.cpim.transport :as tr]))

(defn molecular-resistance
  "A_K - Molecular pressure resistance coefficient
  in (kg/mol)^0.5/s.
  RT in J/mol, inner diameter in m."
  [rt diameter]
  (* (/ 3.0 diameter)
     (Math/sqrt (/ (* Math/PI rt) 8.0))))

(defn prepare-coefficients
  "Constant CPIM coefficients for an isothermal circular tube.
  Temperature in K, inner diameter in m.
  All returned quantities are in SI units."
  [{:keys [gas-const sb4-molar-mass h2-molar-mass]
    :as params}
   temp diameter]
  (let [rt (* gas-const temp)
        m1 sb4-molar-mass
        m2 h2-molar-mass
        mu1 (tr/sb4-viscosity params temp)
        mu2 (tr/h2-viscosity params temp)]
    {:rt       rt
     :diameter diameter
     :ak       (molecular-resistance rt diameter)
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

(defn inverse-continuum-resistance
  "1/A_C - Inverse continuum pressure resistance coefficient
  in s*(mol/kg)^0.5.
  Pressure in Pa, viscosity in Pa*s, RT in J/mol,
  diameter in m, mole fractions dimensionless,
  square roots of molar masses in (kg/mol)^0.5."
  [pressure viscosity rt diameter x1 x2 sqrt-m1 sqrt-m2]
  (/ (* pressure diameter diameter
        (+ (* x1 sqrt-m1)
           (* x2 sqrt-m2)))
     (* 32.0 viscosity rt)))

(defn pressure-resistance
  "A - Interpolated CPIM pressure resistance coefficient
  in (kg/mol)^0.5/s.
  Inverse A_C in s*(mol/kg)^0.5,
  A_K in (kg/mol)^0.5/s."
  [inverse-ac ak]
  (/ 1.0
     (+ inverse-ac
        (/ 1.0 ak))))

(defn local-coefficients
  "Local CPIM coefficients at partial pressures p_1 and p_2.
  Pressures in Pa, all returned quantities in SI units.
  Partial pressures must be nonnegative,
  total pressure must be positive."
  [{:keys [rt diameter sqrt-m1 sqrt-m2
           dk1 dk2 k12 mu1 mu2 phi12 phi21 ak]}
   p1 p2]
  (let [p (+ p1 p2)
        x1 (/ p1 p)
        x2 (/ p2 p)
        mu (tr/mixture-viscosity x1 x2 mu1 mu2 phi12 phi21)
        inverse-ac (inverse-continuum-resistance
                    p mu rt diameter x1 x2 sqrt-m1 sqrt-m2)]
    {:p          p
     :x1         x1
     :x2         x2
     :mu         mu
     :d-eff1     (effective-diffusivity p k12 dk1)
     :d-eff2     (effective-diffusivity p k12 dk2)
     :inverse-ac inverse-ac
     :a          (pressure-resistance inverse-ac ak)}))

(defn pressure-derivatives
  "dp_1/dz, dp_2/dz - CPIM partial pressure derivatives
  in Pa/m.
  Partial pressures in Pa, molar fluxes in mol/(m^2*s)."
  [{:keys [rt sqrt-m1 sqrt-m2] :as coefficients}
   p1 p2 n1 n2]
  (let [{:keys [x1 x2 d-eff1 d-eff2 a]}
        (local-coefficients coefficients p1 p2)

        h (- (/ (* x1 n2) d-eff2)
             (/ (* x2 n1) d-eff1))

        g (- (* a (+ (* sqrt-m1 n1)
                     (* sqrt-m2 n2))))

        rt-h (* rt h)]
    [(+ (* x1 g) rt-h)
     (- (* x2 g) rt-h)]))
