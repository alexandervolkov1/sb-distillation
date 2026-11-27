(ns sb-distillation.cpim.solver
  (:require
   [sb-distillation.cpim.model :as model]
   [sb-distillation.cpim.thermodynamics :as thermo]
   [sb-distillation.numerics.ode :as ode]
   [sb-distillation.numerics.root :as root]))

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

(defn evaluate-trial
  "Evaluate a trial Sb4 molar flux N_1 in mol/(m^2*s).
  F(N_1) = p_1(0; N_1) - p_{1,sat}, in Pa.
  The residual is available only after successful
  integration to z = 0."
  ([problem n1]
   (evaluate-trial problem n1 {}))

  ([{:keys [p1-in] :as problem} n1 options]
   (let [{:keys [status state] :as trial}
         (integrate-trial problem n1 options)]

     (if (= status :success)
       (assoc trial
              :residual (- (first state) p1-in))
       trial))))

(defn molecular-flux-estimate
  "N_{1,mol} - Sb4 molar flux estimate
  in the free-molecular limit, in mol/(m^2*s).
  Used as the initial shooting bracket scale."
  [{:keys [length p1-in coefficients]}]
  (let [{:keys [rt dk1]} coefficients]
    (/ (* dk1 p1-in)
       (* rt length))))

(defn find-shooting-bracket
  "Find a bracket for the shooting residual F(N_1).
  Start at zero flux and double the molecular estimate
  until a nonnegative residual is found.
  Failed integrations terminate the search.
  Fluxes in mol/(m^2*s), residuals in Pa.
  max-trials includes the initial zero-flux trial."
  ([problem]
   (find-shooting-bracket problem {} 32))

  ([problem ode-options max-trials]
   (when-not (and (integer? max-trials)
                  (>= max-trials 2))
     (throw
      (ex-info "Invalid bracket trial limit."
               {:max-trials max-trials})))

   (loop [lower (evaluate-trial problem 0.0 ode-options)
          upper-n1 (molecular-flux-estimate problem)
          trials 1]

     (cond
       (not= :success (:status lower))
       {:status :failed
        :reason :integration-failed
        :trial lower
        :trials trials}

       (not (finite-number? (:residual lower)))
       {:status :failed
        :reason :non-finite-residual
        :trial lower
        :trials trials}

       (not (neg? (:residual lower)))
       {:status :failed
        :reason :invalid-lower-residual
        :trial lower
        :trials trials}

       (>= trials max-trials)
       {:status :failed
        :reason :max-trials
        :lower-trial lower
        :trials trials}

       (not (and (finite-number? upper-n1)
                 (> upper-n1 (:n1 lower))))
       {:status :failed
        :reason :invalid-upper-flux
        :upper-n1 upper-n1
        :lower-trial lower
        :trials trials}

       :else
       (let [upper (evaluate-trial
                    problem upper-n1 ode-options)
             trials (inc trials)]

         (cond
           (not= :success (:status upper))
           {:status :failed
            :reason :integration-failed
            :lower-trial lower
            :trial upper
            :trials trials}

           (not (finite-number? (:residual upper)))
           {:status :failed
            :reason :non-finite-residual
            :trial upper
            :trials trials}

           (>= (:residual upper) 0.0)
           {:status :success
            :bracket [(:n1 lower) (:n1 upper)]
            :lower-trial lower
            :upper-trial upper
            :trials trials}

           :else
           (recur upper
                  (* 2.0 upper-n1)
                  trials)))))))

(def default-solver-options
  {:residual-tol      1.0e-3
   :max-bracket-trials 32
   :max-iterations    64
   :ode-options       {}})

(defn sb-mass-flow
  "m_dot_1 - Antimony mass flow in kg/s.
  Molar flux in mol/(m^2*s), area in m^2,
  Sb4 molar mass in kg/mol."
  [molar-flux area molar-mass]
  (* molar-mass molar-flux area))

(defn solve
  "Solve the CPIM boundary-value problem by shooting.
  Returns the Sb4 molar flux N_1 in mol/(m^2*s),
  the pressure profile in Pa, and solver diagnostics.
  Residual tolerance in Pa."
  ([problem]
   (solve problem {}))

  ([problem options]
   (let [{:keys [residual-tol
                 max-bracket-trials
                 max-iterations
                 ode-options]}
         (merge default-solver-options options)]

     (when-not (and (finite-number? residual-tol)
                    (pos? residual-tol)
                    (integer? max-iterations)
                    (pos? max-iterations))
       (throw
        (ex-info "Invalid shooting solver options."
                 {:options options})))

     (let [{:keys [status
                   bracket
                   lower-trial
                   upper-trial
                   trials]
            :as bracket-result}
           (find-shooting-bracket
            problem ode-options max-bracket-trials)]

       (if (not= status :success)

         (assoc bracket-result :stage :bracketing)

         (let [[a b] bracket

               cache (volatile!
                      {(double a) lower-trial
                       (double b) upper-trial})

               residual-fn
               (fn [n1]
                 (let [trial
                       (or (get @cache n1)
                           (let [result
                                 (evaluate-trial
                                  problem n1 ode-options)]
                             (vswap! cache assoc n1 result)
                             result))]

                   (if (and (= :success (:status trial))
                            (finite-number? (:residual trial)))
                     (:residual trial)
                     ##NaN)))

               root-result
               (root/bisect
                residual-fn a b
                residual-tol max-iterations)

               trial (get @cache (:x root-result))]

           (if (= :success (:status root-result))

             (assoc trial
                    :mass-flow (sb-mass-flow
                                (:n1 trial)
                                (:area problem)
                                (:sb4-molar-mass problem))
                    :root-result root-result
                    :bracket-trials trials)

             {:status :failed
              :stage :shooting
              :reason (cond
                        (not= :success (:status trial))
                        :integration-failed

                        (not (finite-number? (:residual trial)))
                        :non-finite-residual

                        :else
                        (:reason root-result))
              :trial trial
              :root-result root-result
              :bracket-trials trials})))))))

(defn distillation-time
  "t - Distillation time in seconds.
  Antimony mass in kg, mass flow in kg/s.
  Positive mass at zero flow requires infinite time."
  [mass mass-flow]
  (when-not (and (finite-number? mass)
                 (>= mass 0.0)
                 (finite-number? mass-flow)
                 (>= mass-flow 0.0))
    (throw
     (ex-info "Invalid distillation parameters."
              {:mass mass
               :mass-flow mass-flow})))

  (cond
    (zero? mass)      0.0
    (zero? mass-flow) ##Inf
    :else             (/ mass mass-flow)))
