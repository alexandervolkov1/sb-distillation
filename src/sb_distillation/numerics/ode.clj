(ns sb-distillation.numerics.ode
  (:import [java.math BigDecimal MathContext RoundingMode]))

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

(defn- wide-quotient
  [^BigDecimal numerator ^BigDecimal denominator]
  (let [result (.doubleValue (.divide numerator denominator MathContext/DECIMAL128))]
    (if (<= (Math/abs result) Double/MIN_NORMAL)
      ;; Round directly in units of the smallest double to preserve ties.
      (let [unit-denominator (.multiply denominator (BigDecimal. Double/MIN_VALUE))
            units (.divide numerator unit-denominator 0 RoundingMode/HALF_EVEN)]
        (* Double/MIN_VALUE (.doubleValue units)))
      result)))

(defn- rk4-increment
  [h fraction a b c d]
  (let [weight (/ (* h fraction) 6.0)
        weighted-slope (+ a (* 2.0 b) (* 2.0 c) d)]
    (if (and (Double/isFinite weighted-slope)
             (>= (Math/abs (double weight)) Double/MIN_NORMAL))
      (* weight weighted-slope)
      ;; Preserve the product when the slope overflows or the weight underflows.
      (let [inner (.add (BigDecimal. (double b)) (BigDecimal. (double c)))
            outer (.add (BigDecimal. (double a)) (BigDecimal. (double d)))
            slope (.add outer (.multiply inner 2M))
            length (.multiply (BigDecimal. (double h)) (BigDecimal. (double fraction)))]
        (wide-quotient (.multiply length slope) 6M)))))

(defn- rk4-step-at
  [f [z z-mid z-end] y h fraction valid-state?]
  (let [rhs (partial checked-derivatives f valid-state?)
        step (* h fraction)
        half-step (* 0.5 step)
        ;; For tiny steps, apply the fractional weights before scaling by h.
        subnormal-half? (< (Math/abs (double half-step)) Double/MIN_NORMAL)
        half-increment (if subnormal-half?
                         #(* h (* 0.5 fraction %))
                         #(* half-step %))
        full-increment (if subnormal-half?
                         #(* h (* fraction %))
                         #(* step %))
        k1 (rhs z y)
        k2 (rhs z-mid (mapv #(+ %1 (half-increment %2)) y k1))
        k3 (rhs z-mid (mapv #(+ %1 (half-increment %2)) y k2))
        k4 (rhs z-end (mapv #(+ %1 (full-increment %2)) y k3))
        next-y (mapv (fn [yi a b c d] (+ yi (rk4-increment h fraction a b c d)))
                     y k1 k2 k3 k4)]
    (check-state! valid-state? z-end next-y)))

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
   (rk4-step-at f [z (+ z (* 0.5 h)) (+ z h)] y h 1.0 valid-state?)))

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
   (let [z-mid (+ z (* 0.5 h))
         z-end (+ z h)
         coarse (rk4-step-at f [z z-mid z-end] y h 1.0 valid-state?)
         ;; Anchor both half steps to the original coordinate and share z-end.
         middle (rk4-step-at f [z (+ z (* 0.25 h)) z-mid] y h 0.5 valid-state?)
         fine (rk4-step-at f [z-mid (+ z (* 0.75 h)) z-end] middle h 0.5 valid-state?)
         error (mapv (fn [fine-value coarse-value]
                       (let [difference (- fine-value coarse-value)]
                         (if (Double/isInfinite difference)
                           (- (/ fine-value 15.0) (/ coarse-value 15.0))
                           (/ difference 15.0))))
                     fine coarse)]
     {:state fine
      :error error})))

(defn- normalized-component-error
  [scale error abs-tol rel-tol]
  (let [denominator (+ abs-tol (* rel-tol scale))
        wide? (or (Double/isInfinite denominator)
                  (and (< 0.0 denominator Double/MIN_NORMAL)
                       (pos? rel-tol) (pos? scale)))]
    (cond
      (not (Double/isFinite error)) error
      (not wide?) (/ error (double denominator))
      (and (Double/isFinite scale) (Double/isFinite abs-tol) (Double/isFinite rel-tol))
      (let [relative (.multiply (BigDecimal. (double rel-tol)) (BigDecimal. (double scale)))
            wide-denominator (.add (BigDecimal. (double abs-tol)) relative)]
        (wide-quotient (BigDecimal. (double error)) wide-denominator))
      :else Double/NaN)))

(defn normalized-error
  "E - Dimensionless maximum normalized local error.
  y and fine are the initial and fine state vectors.
  error is the estimated local correction vector.
  abs-tol has state units, rel-tol is dimensionless."
  [y fine error abs-tol rel-tol]
  (reduce max 0.0
          (map (fn [yi fi ei]
                 (normalized-component-error
                  (max (Math/abs (double yi)) (Math/abs (double fi)))
                  (Math/abs (double ei)) abs-tol rel-tol))
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

(defn- step-coordinate
  [z z-end step min-step max-step retry?]
  (let [direction (if (< z z-end) 1.0 -1.0)
        remaining (Math/abs (double (- z-end z)))
        toward (if (pos? direction) Double/POSITIVE_INFINITY Double/NEGATIVE_INFINITY)
        proposed (if (>= step remaining) z-end (+ z (* direction step)))
        ;; A reduced proposal can round to zero even when a smaller step exists.
        target (if (and retry? (= proposed z))
                 (Math/nextAfter (double z) (double z-end))
                 proposed)
        h (- target z)
        admissible (fn [candidate]
                     (let [actual-h (- candidate z)
                           magnitude (* direction actual-h)]
                       (when (and (pos? magnitude)
                                  (<= (* direction candidate) (* direction z-end))
                                  (<= magnitude max-step)
                                  (or (= candidate z-end) (>= magnitude min-step))
                                  (= candidate (+ z actual-h)))
                         candidate)))]
    (if (= target z)
      z
      ;; Try neighboring coordinates and displacements if rounding crosses a
      ;; bound or if adding the computed displacement cannot reach the target.
      (or (admissible target)
          (admissible (+ z h))
          (admissible (Math/nextAfter (double target) (double z)))
          (admissible (+ z (Math/nextAfter (double h) 0.0)))
          (admissible (Math/nextAfter (double target) (double toward)))
          (admissible (+ z (Math/nextAfter (double h) (double toward))))))))

(defn integrate
  "Integrate a system of ODEs using adaptive RK4
  with step doubling and local error control.
  Supports forward and backward integration.
  Returns the final state, solution profile,
  attempt count and integration status. A step that cannot
  be represented within the bounds fails with
  :step-size-unrepresentable."
  ([f z-start z-end y-start]
   (integrate f z-start z-end y-start {}))

  ([f z-start z-end y-start options]
   (when-not (and (ifn? f)
                  (finite-number? z-start)
                  (finite-number? z-end)
                  (vector? y-start)
                  (seq y-start)
                  (map? options))
     (throw (ex-info "Invalid integration arguments"
                     {:z-start z-start
                      :z-end z-end
                      :state y-start})))

   (let [z-start (double z-start)
         z-end (double z-end)
         span (Math/abs (- z-end z-start))
         {:keys [abs-tol rel-tol safety
                 min-factor max-factor max-attempts
                 initial-step min-step max-step
                 valid-state?]}
         (prepare-options span options)]

     (check-state! valid-state? z-start y-start)

     (loop [z z-start
            y y-start
            next-z (step-coordinate z-start z-end initial-step min-step max-step false)
            attempts 0
            profile [{:z z-start :state y-start}]]

       (let [result {:z z
                     :state y
                     :attempts attempts
                     :profile profile}]

         (cond
           (= z z-end)
           (assoc result :status :success)

           (>= attempts max-attempts)
           (assoc result
                  :status :failed
                  :reason :max-attempts)

           (= next-z z)
           (assoc result :status :failed :reason :coordinate-stagnation)

           (nil? next-z)
           (assoc result :status :failed :reason :step-size-unrepresentable)

           :else
           (let [h (- next-z z)
                 step (Math/abs (double h))
                 {:keys [status state error] :as trial}
                 (attempt-step f z y h valid-state? abs-tol rel-tol)]
             (if (= status :accepted)
               (let [next-step (next-step-size step error safety min-factor max-factor)
                     proposal (max min-step (min max-step next-step))]
                 (recur next-z
                        state
                        (step-coordinate next-z z-end proposal min-step max-step false)
                        (inc attempts)
                        (conj profile {:z next-z :state state})))
               (let [reduced-step (if (= (:reason trial) :error-too-large)
                                    (next-step-size step error safety min-factor max-factor)
                                    (* 0.5 step))
                     retry-limit (Math/nextDown step)
                     proposal (min retry-limit (max min-step reduced-step))
                     retry-z (step-coordinate z z-end proposal min-step retry-limit true)]
                 (if (or (nil? retry-z) (= z retry-z))
                   (assoc result :status :failed :reason :step-too-small
                          :attempts (inc attempts) :last-rejection trial)
                   (recur z y retry-z (inc attempts) profile)))))))))))
