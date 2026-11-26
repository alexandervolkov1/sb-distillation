(ns sb-distillation.numerics.ode)

(defn- finite-number?
  "Return true if x is a finite number."
  [x]
  (and (number? x)
       (Double/isFinite (double x))))

(defn- check-state!
  "Check the numerical and problem-specific validity
  of a state vector."
  [valid-state? z y]
  (when-not (and (every? finite-number? y)
                 (valid-state? y))
    (throw (ex-info "Invalid integration state"
                    {:type :invalid-state
                     :z z
                     :state y})))
  y)

(defn- checked-derivatives
  "Evaluate the ODE right-hand side after state validation
  and check the resulting derivative vector."
  [f valid-state? z y]
  (check-state! valid-state? z y)
  (let [dy (f z y)]
    (when-not (and (vector? dy)
                   (= (count y) (count dy))
                   (every? number? dy))
      (throw (ex-info "Invalid derivative vector"
                      {:type :invalid-derivative-vector
                       :z z
                       :state y
                       :derivatives dy})))
    (when-not (every? finite-number? dy)
      (throw (ex-info "Non-finite derivatives"
                      {:type :non-finite-derivatives
                       :z z
                       :state y
                       :derivatives dy})))
    dy))

(defn rk4-step
  "y_{n+1} - Advance the state by one classical fourth-order
  Runge-Kutta step.
  f takes the independent coordinate and state vector
  and returns a derivative vector of the same dimension.
  h is a signed integration step.
  valid-state? checks problem-specific state constraints."
  ([f z y h]
   (rk4-step f z y h (constantly true)))

  ([f z y h valid-state?]
   (let [rhs (partial checked-derivatives f valid-state?)
         half-h (* 0.5 h)
         k1 (rhs z y)
         k2 (rhs (+ z half-h)
                 (mapv #(+ %1 (* half-h %2)) y k1))
         k3 (rhs (+ z half-h)
                 (mapv #(+ %1 (* half-h %2)) y k2))
         k4 (rhs (+ z h)
                 (mapv #(+ %1 (* h %2)) y k3))
         next-y (mapv
                 (fn [yi a b c d]
                   (+ yi (* (/ h 6.0)
                            (+ a (* 2.0 b) (* 2.0 c) d))))
                 y k1 k2 k3 k4)]
     (check-state! valid-state? (+ z h) next-y))))

(defn rk4-step-doubling
  "y_f, e - RK4 step with local error estimation
  using step doubling.
  Returns the fine state and the estimated signed
  correction (y_f - y_c)/15.
  h is a signed integration step.
  valid-state? checks problem-specific state constraints."
  ([f z y h]
   (rk4-step-doubling f z y h (constantly true)))

  ([f z y h valid-state?]
   (let [half-h (* 0.5 h)
         coarse (rk4-step f z y h valid-state?)
         middle (rk4-step f z y half-h valid-state?)
         fine (rk4-step f (+ z half-h)
                        middle half-h valid-state?)
         error (mapv #(/ (- %1 %2) 15.0) fine coarse)]
     {:state fine
      :error error})))

(defn normalized-error
  "E - Dimensionless maximum normalized local error.
  y and fine are the initial and fine state vectors.
  error is the estimated local correction vector.
  abs-tol has state units, rel-tol is dimensionless."
  [y fine error abs-tol rel-tol]
  (reduce max 0.0
          (map (fn [yi fi ei]
                 (/ (Math/abs (double ei))
                    (+ abs-tol
                       (* rel-tol
                          (max (Math/abs (double yi))
                               (Math/abs (double fi)))))))
               y fine error)))

(defn next-step-size
  "h_new - Adjust the signed RK4 integration step
  using the normalized local error E.
  Safety and step adjustment factors are dimensionless.
  Error must be finite and nonnegative."
  [h error safety min-factor max-factor]
  (let [factor (if (zero? error)
                 max-factor
                 (* safety (Math/pow error -0.2)))]
    (* h (max min-factor
              (min max-factor factor)))))

(defn attempt-step
  "Attempt one adaptive RK4 integration step.
  Returns an accepted state or a rejection reason.
  Programming errors are propagated."
  [f z y h valid-state? abs-tol rel-tol]
  (try
    (let [{:keys [state error]}
          (rk4-step-doubling f z y h valid-state?)

          e (normalized-error
             y state error abs-tol rel-tol)]

      (cond
        (not (finite-number? e))
        {:status :rejected
         :reason :non-finite-error}

        (<= e 1.0)
        {:status :accepted
         :state state
         :error e}

        :else
        {:status :rejected
         :reason :error-too-large
         :error e}))

    (catch clojure.lang.ExceptionInfo ex
      (let [data (ex-data ex)]
        (if (#{:invalid-state
               :non-finite-derivatives} (:type data))
          {:status :rejected
           :reason (:type data)
           :details data}
          (throw ex))))))

(def default-options
  {:abs-tol      1.0e-8
   :rel-tol      1.0e-6
   :safety       0.9
   :min-factor   0.2
   :max-factor   5.0
   :max-attempts 100000
   :valid-state? (constantly true)})

(defn- prepare-options
  "Prepare and validate adaptive integration options.
  span is the positive integration interval length.
  Step sizes are positive magnitudes."
  [span options]
  (when-not (and (finite-number? span) (pos? span))
    (throw (ex-info "Invalid integration interval"
                    {:span span})))

  (let [{:keys [abs-tol rel-tol safety
                min-factor max-factor max-attempts
                initial-step min-step max-step
                valid-state?]
         :as opts}
        (merge default-options
               {:initial-step (/ span 100.0)
                :min-step     (* span 1.0e-12)
                :max-step     span}
               options)]

    (when-not
     (and
      (every? finite-number?
              [abs-tol rel-tol safety
               min-factor max-factor
               initial-step min-step max-step])
      (pos? abs-tol)
      (>= rel-tol 0.0)
      (< 0.0 safety 1.0)
      (< 0.0 min-factor 1.0)
      (>= max-factor 1.0)
      (< 0.0 min-step)
      (<= min-step initial-step max-step)
      (integer? max-attempts)
      (pos? max-attempts)
      (ifn? valid-state?))
      (throw (ex-info "Invalid integration options"
                      {:options opts})))

    opts))
