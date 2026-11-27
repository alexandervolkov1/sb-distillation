(ns sb-distillation.cpim.solver
  (:require
   [sb-distillation.cpim.model :as model]
   [sb-distillation.cpim.thermodynamics :as thermo]
   [sb-distillation.numerics.ode :as ode]))

(defn- finite-number?
  "Return true if x is a finite number
  representable as a double."
  [x]
  (and (number? x)
       (Double/isFinite (double x))))

(defn valid-pressure-state?
  "Check physical admissibility of the pressure
  state [p_1 p_2], with pressures in Pa.
  Both partial pressures and total pressure
  must be finite.
  Requires p_1 >= 0 and p_2 > 0."
  [state]
  (and (sequential? state)
       (= 2 (count state))
       (let [[p1 p2] state]
         (and (finite-number? p1)
              (finite-number? p2)
              (>= p1 0.0)
              (> p2 0.0)
              (Double/isFinite
               (+ (double p1) (double p2)))))))

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

(defn- validate-conditions!
  "Validate operating conditions of the CPIM
  boundary-value problem.
  All quantities are in SI units."
  [conditions]
  (doseq [k [:temp
             :diameter
             :length
             :ref-pressure
             :ref-temp]]
    (let [v (get conditions k)]
      (when-not (and (finite-number? v)
                     (pos? v))
        (throw
         (ex-info "Expected a finite positive parameter."
                  {:parameter k :value v})))))

  (doseq [k [:ref-volumetric-flow
             :gauge-pressure
             :trap-pressure-drop]]
    (let [v (get conditions k)]
      (when-not (and (finite-number? v)
                     (not (neg? v)))
        (throw
         (ex-info "Expected a finite nonnegative parameter."
                  {:parameter k :value v})))))

  (when-not (> (:temp conditions) 152.352)
    (throw
     (ex-info "Temperature outside Antoine correlation domain."
              {:temp (:temp conditions)})))

  (let [p2-out (outlet-hydrogen-pressure
                (:gauge-pressure conditions)
                (:trap-pressure-drop conditions))]
    (when-not (valid-pressure-state? [0.0 p2-out])
      (throw
       (ex-info "Invalid outlet pressure state."
                {:outlet-state [0.0 p2-out]}))))

  conditions)

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
           trap-pressure-drop]
    :as conditions}]
  (validate-conditions! conditions)
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

(defn integrate-trial
  "Integrate a trial CPIM pressure profile
  backward from z = L to z = 0.
  N_1 is the trial Sb4 molar flux in mol/(m^2*s).
  Coordinates in m, pressures in Pa."
  ([problem n1]
   (integrate-trial problem n1 {}))

  ([{:keys [length
            coefficients
            n2
            outlet-state]}
    n1 options]

   (when-not (and (finite-number? n1)
                  (>= n1 0.0))
     (throw
      (ex-info "Invalid trial Sb4 molar flux."
               {:n1 n1})))

   (let [rhs (fn [_z [p1 p2]]
               (model/pressure-derivatives
                coefficients p1 p2 n1 n2))]

     (assoc
      (ode/integrate
       rhs
       length
       0.0
       outlet-state
       (assoc options
              :valid-state? valid-pressure-state?))
      :n1 n1))))
