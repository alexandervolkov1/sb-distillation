(ns sb-distillation.cpim.solver
  (:require
   [sb-distillation.cpim.model :as model]
   [sb-distillation.cpim.thermodynamics :as thermo]))

(defn cross-sectional-area
  "S - Tube cross-sectional area in m^2.
  Inner diameter in m."
  [diameter]
  (/ (* Math/PI diameter diameter) 4.0))

(defn hydrogen-throughput
  "Q_2 - Hydrogen throughput in Pa*m^3/s
  under the reference-condition convention.
  Reference pressure in Pa,
  reference volumetric flow in m^3/s."
  [ref-pressure ref-volumetric-flow]
  (* ref-pressure ref-volumetric-flow))

(defn hydrogen-molar-flow
  "n_dot_2 - Hydrogen molar flow in mol/s.
  Throughput in Pa*m^3/s,
  reference temperature in K."
  [{:keys [gas-const]} throughput ref-temp]
  (/ throughput (* gas-const ref-temp)))

(defn hydrogen-molar-flux
  "N_2 - Hydrogen molar flux in mol/(m^2*s).
  Molar flow in mol/s, inner diameter in m."
  [molar-flow diameter]
  (/ molar-flow (cross-sectional-area diameter)))

(defn trap-pressure-drop
  "Delta p_trap - Hydrogen pressure drop
  across the trap in Pa.
  Throughput in Pa*m^3/s, conductance in m^3/s.
  Both must use the same reference conditions."
  [throughput trap-conductance]
  (/ throughput trap-conductance))

(defn outlet-hydrogen-pressure
  "p_2(L) - Hydrogen partial pressure
  at the tube outlet in Pa.
  Gauge pressure and trap pressure drop in Pa."
  [gauge-pressure pressure-drop]
  (+ gauge-pressure pressure-drop))

(defn prepare-problem
  "Prepare constant coefficients and boundary conditions
  for the CPIM boundary-value problem.
  All input and output quantities are in SI units."
  [params
   {:keys [temp
           diameter
           length
           ref-pressure
           ref-temp
           ref-volumetric-flow
           gauge-pressure
           trap-pressure-drop]}]
  (let [area (cross-sectional-area diameter)

        throughput (hydrogen-throughput
                    ref-pressure
                    ref-volumetric-flow)

        molar-flow (hydrogen-molar-flow
                    params throughput ref-temp)

        n2 (hydrogen-molar-flux molar-flow diameter)

        p1-in (thermo/sb-liquid-vapor-pressure temp)

        p2-out (outlet-hydrogen-pressure
                gauge-pressure trap-pressure-drop)

        coefficients (model/prepare-coefficients
                      params temp diameter)]

    {:temp                temp
     :diameter            diameter
     :length              length
     :area                area
     :sb4-molar-mass       (:sb4-molar-mass params)
     :hydrogen-throughput throughput
     :hydrogen-molar-flow molar-flow
     :n2                  n2
     :p1-in               p1-in
     :outlet-state        [0.0 p2-out]
     :coefficients        coefficients}))
