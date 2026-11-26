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
  "D_{K,i} - Knudsen diffusivity in m^2/s.
  Temperature in K, diameter in m, molar mass [kg/mol]."
  [{:keys [gas-const]} temp diameter molar-mass]
  (* (/ diameter 3.0)
     (Math/sqrt
      (/ (* 8.0 gas-const temp)
         (* Math/PI molar-mass)))))

(defn binary-diffusivity-factor
  "K_12 = p * D_12 - Pressure-independent binary diffusion factor
  [Pa*m^2/s], using the hard-sphere Chapman-Enskog model."
  [{:keys [bolzman-const
           avogadro-const
           sb4-molar-mass
           h2-molar-mass
           sb4-collision-diameter
           h2-collision-diameter]}
   temp]
  (let [sigma (/ (+ sb4-collision-diameter
                    h2-collision-diameter)
                 2.0)
        reduced-mass (/ (* sb4-molar-mass h2-molar-mass)
                        (* avogadro-const
                           (+ sb4-molar-mass h2-molar-mass)))
        kt (* bolzman-const temp)]
    (* (/ (* 3.0 kt)
          (* 16.0 sigma sigma))
       (Math/sqrt
        (/ (* 2.0 kt)
           (* Math/PI reduced-mass))))))

(defn binary-diffusivity
  "D_12 - Binary diffusivity in m^2/s.
  Pressure in Pa, diffusion factor [Pa*m^2/s]."
  [factor pressure]
  (/ factor pressure))
