(ns sb-distillation.cpim.solver)

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
