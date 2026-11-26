(ns sb-distillation.cpim.transport)

(def default-parameters
  {:gas-const           8.314462618
   :boltzmann-const     1.380649e-23
   :avogadro-const      6.02214076e23

   :sb4-molar-mass         0.48704
   :h2-molar-mass          0.00201588

   :sb4-collision-diameter 0.600e-9
   :h2-collision-diameter  0.289e-9

   :h2-reference-viscosity 8.76e-6
   :h2-reference-temp      293.85
   :h2-sutherland-constant 72.0})

(defn knudsen-diffusivity
  "Knudsen diffusivity in m^2/s.
  Temperature in K, diameter in m, molar mass in kg/mol."
  [{:keys [gas-const]} temp diameter molar-mass]
  (* (/ diameter 3.0)
     (Math/sqrt
      (/ (* 8.0 gas-const temp)
         (* Math/PI molar-mass)))))
