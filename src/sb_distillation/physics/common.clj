(ns sb-distillation.physics.common
  (:require [sb-distillation.physics.transport :as tr]))

(defn molecular-resistance
  "A_K - Molecular pressure resistance coefficient
  in (kg/mol)^0.5/s.
  RT in J/mol, inner diameter in m."
  [rt diameter]
  (* (/ 3.0 diameter)
     (Math/sqrt (/ (* Math/PI rt) 8.0))))

(defn- checked-coefficient [quantity compute]
  (let [value (try
                (compute)
                (catch ArithmeticException cause
                  (throw (ex-info (str "Cannot compute tube coefficient " (name quantity) ".")
                                  {:quantity quantity} cause))))
        x (double value)]
    (when-not (and (pos? x) (Double/isFinite x))
      (throw (ex-info (str "Expected a finite positive tube coefficient " (name quantity) ".")
                      {:quantity quantity :value value})))
    value))

(defn prepare-coefficients
  "Shared constant transport coefficients for an isothermal circular tube.
  Temperature in K, inner diameter in m.
  All returned quantities are in SI units.
  Throws if a computed coefficient is outside the positive finite double range."
  [{:keys [gas-const sb4-molar-mass h2-molar-mass]
    :as params}
   temp diameter]
  (let [rt (checked-coefficient :rt #(* gas-const temp))
        m1 sb4-molar-mass
        m2 h2-molar-mass
        mu1 (checked-coefficient :mu1 #(tr/sb4-viscosity params temp))
        mu2 (checked-coefficient :mu2 #(tr/h2-viscosity params temp))]
    {:rt       rt
     :diameter diameter
     :ak       (checked-coefficient :ak #(molecular-resistance rt diameter))
     :sqrt-m1  (checked-coefficient :sqrt-m1 #(Math/sqrt m1))
     :sqrt-m2  (checked-coefficient :sqrt-m2 #(Math/sqrt m2))
     :dk1      (checked-coefficient :dk1 #(tr/knudsen-diffusivity params temp diameter m1))
     :dk2      (checked-coefficient :dk2 #(tr/knudsen-diffusivity params temp diameter m2))
     :k12      (checked-coefficient :k12 #(tr/binary-diffusivity-factor params temp))
     :mu1      mu1
     :mu2      mu2
     :phi12    (checked-coefficient :phi12 #(tr/wilke-interaction-factor mu1 mu2 m1 m2))
     :phi21    (checked-coefficient :phi21 #(tr/wilke-interaction-factor mu2 mu1 m2 m1))}))

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
  "A - Combined continuum and molecular pressure resistance coefficient
  in (kg/mol)^0.5/s.
  Inverse A_C in s*(mol/kg)^0.5,
  A_K in (kg/mol)^0.5/s."
  [inverse-ac ak]
  (/ 1.0
     (+ inverse-ac
        (/ 1.0 ak))))

(defn local-coefficients
  "p, x_i, mu, 1/A_C, and A - Local pressure in Pa, dimensionless mole
  fractions, viscosity in Pa*s, and pressure resistances in SI units.
  Partial pressures in Pa.
  Partial pressures must be nonnegative,
  total pressure must be positive."
  [{:keys [rt diameter sqrt-m1 sqrt-m2 mu1 mu2 phi12 phi21 ak]}
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
     :inverse-ac inverse-ac
     :a          (pressure-resistance inverse-ac ak)}))

(defn total-pressure-gradient
  "G = dp/dz - Total pressure gradient in Pa/m, with z increasing
  from source to receiver. A in (kg/mol)^0.5/s,
  molar fluxes in mol/(m^2*s)."
  [{:keys [sqrt-m1 sqrt-m2]} a n1 n2]
  (- (* a (+ (* sqrt-m1 n1)
             (* sqrt-m2 n2)))))
