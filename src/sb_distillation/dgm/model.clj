(ns sb-distillation.dgm.model
  (:require [sb-distillation.physics.common :as physics]
            [sb-distillation.physics.transport :as transport]))

(defn prepare-coefficients
  "B_0 - Add the viscous permeability in m^2 to the shared tube coefficients.
  Inner diameter in m."
  [{:keys [diameter] :as coefficients}]
  (let [b0 (/ (* diameter diameter) 32.0)]
    (when-not (and (pos? b0) (Double/isFinite b0))
      (throw (ex-info "Expected finite positive DGM permeability."
                      {:quantity :b0 :value b0})))
    (assoc coefficients :b0 b0)))

(defn pressure-derivatives
  "dp_1/dz, dp_2/dz - DGM partial pressure derivatives in Pa/m for z
  increasing from source to receiver. Pressures in Pa, molar fluxes
  in mol/(m^2*s). Coefficients include B_0 in m^2 from prepare-coefficients.
  Both partial pressures are nonnegative and total pressure is positive."
  [{:keys [rt dk1 dk2 k12 b0] :as coefficients} p1 p2 n1 n2]
  (let [{:keys [p x1 x2 mu a]} (physics/local-coefficients coefficients p1 p2)
        d12 (transport/binary-diffusivity k12 p)
        g (physics/total-pressure-gradient coefficients a n1 n2)
        mutual (/ (- (* x2 n1) (* x1 n2)) d12)
        p-b0-over-mu (/ (* p b0) mu)
        beta1 (/ (* x1 p-b0-over-mu) dk1)
        beta2 (/ (* x2 p-b0-over-mu) dk2)]
    [(- (* (- rt) (+ (/ n1 dk1) mutual)) (* beta1 g))
     (- (* (- rt) (- (/ n2 dk2) mutual)) (* beta2 g))]))
