(ns sb-distillation.numerics.ode-test
  (:require [clojure.test :refer [deftest is testing]]
            [sb-distillation.numerics.ode :as ode]))

(defn- close?
  [expected actual tolerance]
  (and (number? actual)
       (<= (Math/abs (- (double actual) (double expected))) tolerance)))

(defn- state-close?
  [expected actual tolerance]
  (and (vector? actual)
       (= (count expected) (count actual))
       (every? true? (map #(close? %1 %2 tolerance) expected actual))))

(defn- exponential [_ [y]] [y])
(defn- decay [_ [y]] [(- y)])
(defn- oscillator [_ [x v]] [v (- x)])
(defn- oscillator-solution [z] [(Math/cos z) (- (Math/sin z))])
(defn- variable-rate [z [y]] [(* 2.0 z y)])
(defn- quartic-rhs [z _] [(* z z z z)])

(defn- observed-orders
  [errors]
  (mapv (fn [[coarse fine]]
          (/ (Math/log (/ coarse fine)) (Math/log 2.0)))
        (partition 2 1 errors)))

(defn- exception-info-data
  [call]
  (try
    (call)
    nil
    (catch clojure.lang.ExceptionInfo ex
      (ex-data ex))))

(deftest rk4-step-analytical-solutions-test
  (testing "Exponential growth and decay in both coordinate directions"
    (doseq [[rhs rate] [[exponential 1.0] [decay -1.0]]
            h [0.125 -0.125]]
      (is (state-close? [(Math/exp (* rate h))]
                        (ode/rk4-step rhs 0.5 [1.0] h) 3.0e-7))))
  (testing "Two-dimensional harmonic oscillator"
    (doseq [h [0.125 -0.125]]
      (is (state-close? (oscillator-solution (+ 0.5 h))
                        (ode/rk4-step oscillator 0.5 (oscillator-solution 0.5) h)
                        3.0e-7))))
  (testing "Non-autonomous equation y' = 2zy with y(0.5) = 1"
    (let [target 0.5625
          exact (Math/exp (- (* target target) 0.25))]
      (is (state-close? [exact] (ode/rk4-step variable-rate 0.5 [1.0] 0.0625)
                        1.0e-7)))))

(deftest rk4-step-exact-quadrature-test
  (doseq [h [0.5 -0.5]]
    (testing "Constant vector derivative"
      (is (state-close? [(+ 1.0 (* 2.0 h)) (- (* 3.0 h))]
                        (ode/rk4-step (constantly [2.0 -3.0]) 1.0 [1.0 0.0] h)
                        1.0e-14)))
    (testing "A cubic derivative integrates exactly to a quartic"
      (let [rhs (fn [z _] [(* 4.0 z z z)])]
        (is (state-close? [(Math/pow (+ 1.0 h) 4.0)]
                          (ode/rk4-step rhs 1.0 [1.0] h) 1.0e-14))))))

(deftest rk4-convergence-test
  (testing "One-step error on y' = y is fifth order"
    (let [errors (for [h [0.25 0.125 0.0625 0.03125]]
                   (Math/abs (- (Math/exp h)
                                (first (ode/rk4-step exponential 0.0 [1.0] h)))))]
      (doseq [order (observed-orders errors)]
        (is (< 4.9 order 5.15) (str "Local order " order)))))
  (testing "Fixed-step error at z = 1 is fourth order"
    (let [errors (for [steps [8 16 32 64]
                       :let [h (/ 1.0 steps)
                             state (reduce (fn [y i] (ode/rk4-step exponential (* i h) y h))
                                           [1.0] (range steps))]]
                   (Math/abs (- Math/E (first state))))]
      (doseq [order (observed-orders errors)]
        (is (< 3.85 order 4.1) (str "Global order " order))))))

(deftest rk4-step-doubling-quadrature-test
  ;; For y' = z^4, each classical RK4 step is Simpson quadrature.
  ;; On the two half intervals the quadratures are 5/768 and 149/768.
  ;; The single full step is 5/24 and the exact integral is 1/5.
  (let [fine (+ 5/768 149/768)
        correction (/ (- fine 5/24) 15)
        {:keys [state error]} (ode/rk4-step-doubling quartic-rhs 0.0 [0.0] 1.0)
        backward (ode/rk4-step-doubling quartic-rhs 1.0 [0.2] -1.0)]
    (is (state-close? [fine] state 1.0e-14))
    (is (state-close? [correction] error 1.0e-15))
    (is (neg? (first error)))
    (is (close? (- 0.2 (first state)) (first error) 1.0e-15))
    (is (state-close? [(- 1/5 fine)] (:state backward) 1.0e-14))
    (is (state-close? [(- correction)] (:error backward) 1.0e-15))
    (is (pos? (first (:error backward))))))

(deftest rk4-step-doubling-convergence-test
  (let [estimates (for [h [0.5 0.25 0.125 0.0625]
                        :let [{:keys [state error]} (ode/rk4-step-doubling exponential 0.0 [1.0] h)
                              actual-correction (- (Math/exp h) (first state))]]
                    [(first error) actual-correction])
        residuals (map (fn [[estimated actual]] (Math/abs (- estimated actual))) estimates)]
    (doseq [[estimated actual] estimates]
      (is (and (pos? actual) (< 0.79 (/ estimated actual) 1.01))))
    (doseq [order (observed-orders (map first estimates))]
      (is (< 4.9 order 5.15) (str "Estimator order " order)))
    (doseq [order (observed-orders residuals)]
      (is (< 5.85 order 6.25) (str "Correction discrepancy order " order)))))

(deftest normalized-error-test
  (testing "Absolute tolerance sets the scale when both state components vanish"
    (is (close? 0.75 (ode/normalized-error [0.0] [0.0] [-0.375] 0.5 0.1) 1.0e-14)))
  (testing "The relative term uses the larger magnitude of the initial and fine components"
    (is (close? 0.5 (ode/normalized-error [-10.0] [2.0] [1.5] 1.0 0.2) 1.0e-14))
    (is (close? 0.5 (ode/normalized-error [2.0] [-10.0] [-1.5] 1.0 0.2) 1.0e-14)))
  (testing "Maximum component error and exact zero error"
    (is (close? 2.0 (ode/normalized-error [0.0 10.0 2.0] [0.0 8.0 3.0]
                                          [0.25 -3.0 0.0] 0.5 0.1)
                1.0e-14))
    (is (zero? (ode/normalized-error [0.0 4.0] [1.0 -5.0] [0.0 -0.0] 0.5 0.1)))))

(deftest next-step-size-test
  (let [h 2.0]
    (testing "Fifth-root controller below, at and above the error threshold"
      (doseq [[error expected] [[0.03125 3.6] [1.0 1.8] [32.0 0.9]]]
        (is (close? expected (ode/next-step-size h error 0.9 0.2 5.0) 1.0e-14))))
    (testing "Zero error and both factor clamps, with a signed backward step"
      (doseq [[error expected] [[0.0 8.0] [1.0e-10 8.0] [1.0e10 0.5]]]
        (is (close? expected (ode/next-step-size h error 0.8 0.25 4.0) 1.0e-14))
        (is (close? (- expected) (ode/next-step-size (- h) error 0.8 0.25 4.0) 1.0e-14))))
    (testing "The custom safety factor applies away from the clamps"
      (is (close? -1.4 (ode/next-step-size (- h) 1.0 0.7 0.3 3.0) 1.0e-14)))))

(deftest rk4-state-validation-test
  (let [positive-state? (fn [[y]] (pos? y))]
    (doseq [[label initial derivative expected-z calls]
            [["Initial state" [-1.0] (constantly -3.0) 0.0 0]
             ["Second stage" [1.0] (constantly -3.0) 0.5 1]
             ["Third stage" [1.0] #(if (zero? %) 0.0 -3.0) 0.5 2]
             ["Fourth stage" [1.0] #(if (zero? %) 0.0 -1.5) 1.0 3]
             ["Final state" [1.0] #(if (= % 1.0) -12.0 0.0) 1.0 4]]]
      (testing label
        (let [visited (atom [])
              rhs (fn [z y] (swap! visited conj y) [(derivative z)])
              data (exception-info-data #(ode/rk4-step rhs 0.0 initial 1.0 positive-state?))]
          (is (= :invalid-state (:type data)))
          (is (= expected-z (:z data)))
          (is (neg? (first (:state data))))
          (is (= calls (count @visited)))
          (is (every? positive-state? @visited)))))))

(deftest rk4-non-finite-state-and-derivatives-test
  (doseq [value [##NaN ##Inf ##-Inf]]
    (let [called? (atom false)
          data (exception-info-data
                #(ode/rk4-step (fn [_ _] (reset! called? true) [0.0]) 0.0 [value] 1.0))]
      (is (= :invalid-state (:type data)))
      (is (false? @called?)))
    (let [data (exception-info-data
                #(ode/rk4-step (constantly [value]) 0.0 [0.0] 1.0))]
      (is (= :non-finite-derivatives (:type data)))
      (is (= 0.0 (:z data)))
      (is (= [0.0] (:state data)))
      (is (false? (Double/isFinite (double (first (:derivatives data))))))))
  (testing "Overflow in an intermediate state is detected before the next RHS call"
    (let [calls (atom 0)
          rhs (fn [_ _] (swap! calls inc) [1.0e308])
          data (exception-info-data #(ode/rk4-step rhs 0.0 [1.0e308] 2.0))]
      (is (= {:type :invalid-state :z 1.0 :state [##Inf]} data))
      (is (= 1 @calls)))))

(deftest derivative-vector-and-exception-propagation-test
  (doseq [derivatives [nil 1.0 '(1.0) [] [1.0 2.0] ["1.0"]]]
    (let [rhs (constantly derivatives)
          expected {:type :invalid-derivative-vector :z 0.0 :state [1.0] :derivatives derivatives}]
      (is (= expected (exception-info-data #(ode/rk4-step rhs 0.0 [1.0] 1.0))))
      (is (= expected (exception-info-data #(ode/attempt-step rhs 0.0 [1.0] 1.0
                                                               (constantly true) 1.0e-8 0.0))))))
  (doseq [ex [(ex-info "RHS programming error" {:type :rhs-bug})
              (IllegalStateException. "Unexpected RHS failure")]]
    (let [rhs (fn [_ _] (throw ex))]
      (doseq [call [#(ode/attempt-step rhs 0.0 [1.0] 1.0 (constantly true) 1.0e-8 0.0)
                    #(ode/integrate rhs 0.0 1.0 [1.0])]]
        (is (identical? ex (try (call) nil (catch Exception actual actual))))))))

(deftest attempt-step-test
  (let [attempt (fn [rhs state step valid? abs-tol]
                  (ode/attempt-step rhs 0.0 state step valid? abs-tol 0.0))
        fine (+ 5/768 149/768)
        magnitude 1/1920
        accepted (attempt quartic-rhs [0.0] 1.0 (constantly true) 1.0e-3)
        rejected (attempt quartic-rhs [0.0] 1.0 (constantly true) 1.0e-4)]
    (testing "Accepted fine state and normalized correction, followed by a tolerance rejection"
      (is (= :accepted (:status accepted)))
      (is (state-close? [fine] (:state accepted) 1.0e-14))
      (is (close? (/ magnitude 1.0e-3) (:error accepted) 1.0e-13))
      (is (= {:status :rejected :reason :error-too-large} (dissoc rejected :error)))
      (is (close? (/ magnitude 1.0e-4) (:error rejected) 1.0e-13)))
    (testing "Numerical rejection reasons retain the failing stage diagnostics"
      (is (= {:status :rejected :reason :invalid-state
              :details {:type :invalid-state :z 0.5 :state [-0.5]}}
             (attempt (constantly [-3.0]) [1.0] 1.0 #(pos? (first %)) 1.0e-8)))
      (is (= {:status :rejected :reason :non-finite-derivatives
              :details {:type :non-finite-derivatives :z 0.0 :state [1.0] :derivatives [##Inf]}}
             (attempt (constantly [##Inf]) [1.0] 1.0 (constantly true) 1.0e-8))))
    (testing "A finite correction divided by a subnormal tolerance can overflow"
      (is (= {:status :rejected :reason :non-finite-error}
             (attempt quartic-rhs [0.0] 1.0 (constantly true) Double/MIN_VALUE))))))

(deftest adaptive-analytical-solutions-test
  (let [options {:abs-tol 1.0e-11 :rel-tol 1.0e-9}]
    (doseq [[label rhs exact] [["Exponential" exponential (fn [z] [(Math/exp z)])]
                              ["Harmonic oscillator" oscillator oscillator-solution]
                              ["Variable growth rate" variable-rate (fn [z] [(Math/exp (* z z))])]]
            [start end] [[0.0 1.03] [1.03 0.0]]]
      (testing (str label " from " start " to " end)
        (let [initial (exact start)
              {:keys [status z state profile attempts]} (ode/integrate rhs start end initial options)
              ordered? (if (< start end) < >)]
          (is (= :success status))
          (is (= end z))
          (is (= {:z start :state initial} (first profile)))
          (is (= {:z end :state state} (peek profile)))
          (is (apply ordered? (map :z profile)))
          (is (>= attempts (dec (count profile))))
          (doseq [{:keys [z state]} profile]
            (is (state-close? (exact z) state 2.0e-7) (str "Analytical state at " z))))))))

(deftest adaptive-tolerance-and-step-size-test
  (testing "Tighter tolerances improve agreement with exp(1) and require more trials"
    (let [results (mapv (fn [tol]
                          (ode/integrate exponential 0.0 1.0 [1.0]
                                         {:initial-step 1.0 :abs-tol (* 0.1 tol) :rel-tol tol}))
                        [1.0e-4 1.0e-6 1.0e-8 1.0e-10])
          errors (mapv #(Math/abs (- Math/E (first (:state %)))) results)]
      (is (every? #(= :success (:status %)) results))
      (is (apply > errors) (str "Endpoint errors " errors))
      (is (apply < (map :attempts results)))
      (is (< (peek errors) 1.0e-8))))
  (testing "Accepted step magnitudes change before the final endpoint clipping"
    (let [{:keys [status profile]} (ode/integrate exponential 0.0 1.0 [1.0]
                                                  {:initial-step 0.01 :abs-tol 1.0e-7 :rel-tol 0.0})
          interior-steps (map (fn [[a b]] (- (:z b) (:z a))) (partition 2 1 (pop profile)))]
      (is (= :success status))
      (is (<= 2 (count interior-steps)))
      (is (> (apply max interior-steps) (* 2.0 (apply min interior-steps)))))))

(deftest adaptive-rejection-and-attempt-count-test
  ;; Simpson's local correction has magnitude h^5/1920 for y' = z^4.
  ;; At h=1 the normalized error is 16. At h=1/2 it is 1/2.
  (let [options {:initial-step 1.0 :min-step 0.5
                 :abs-tol (/ 1.0 30720.0) :rel-tol 0.0
                 :safety 0.5 :min-factor 0.5 :max-factor 1.0}
        run (fn [limit] (ode/integrate quartic-rhs 0.0 1.0 [0.0] (assoc options :max-attempts limit)))
        rejected (run 1)
        partial (run 2)
        complete (run 3)]
    (testing "The first rejected full step consumes one attempt and adds no profile entry"
      (is (= {:status :failed :reason :max-attempts :attempts 1 :z 0.0 :state [0.0]
              :profile [{:z 0.0 :state [0.0]}]} rejected)))
    (testing "The next attempt succeeds and is retained when the attempt limit is reached"
      (is (= :failed (:status partial)))
      (is (= :max-attempts (:reason partial)))
      (is (= 2 (:attempts partial)))
      (is (= [0.0 0.5] (mapv :z (:profile partial))))
      (is (= (select-keys partial [:z :state]) (peek (:profile partial))))
      (is (state-close? [(+ 1/160 1/61440)] (:state partial) 1.0e-14)))
    (testing "Two accepted steps plus one rejection reach the target on the final allowed attempt"
      (is (= :success (:status complete)))
      (is (= 3 (:attempts complete)))
      (is (= [0.0 0.5 1.0] (mapv :z (:profile complete))))
      (is (state-close? [1/5] (:state complete) 4.0e-5)))))

(deftest adaptive-invalid-stage-recovery-test
  ;; The first coarse midpoint of y'=-2y at h=1 is zero.
  ;; The exact trajectory is strictly positive throughout the interval.
  (let [invalid-states (atom [])
        rhs-states (atom [])
        valid? (fn [[y :as state]]
                 (when-not (pos? y) (swap! invalid-states conj state))
                 (pos? y))
        rhs (fn [_ [y :as state]]
              (swap! rhs-states conj state)
              [(* -2.0 y)])
        {:keys [status state profile attempts]}
        (ode/integrate rhs 0.0 1.0 [1.0] {:initial-step 1.0 :valid-state? valid?})]
    (is (= :success status))
    (is (some #(zero? (first %)) @invalid-states))
    (is (seq @rhs-states))
    (is (every? #(pos? (first %)) @rhs-states))
    (is (> attempts (dec (count profile))))
    (is (state-close? [(Math/exp -2.0)] state 1.0e-6))))

(deftest adaptive-step-bounds-test
  (testing "Growth of zero-error steps is capped by max-step"
    (let [{:keys [status state profile attempts]}
          (ode/integrate (constantly [1.0]) 0.0 1.0 [0.0]
                         {:initial-step 0.125 :max-step 0.25})]
      (is (= :success status))
      (is (= [0.0 0.125 0.375 0.625 0.875 1.0] (mapv :z profile)))
      (is (= 5 attempts))
      (is (state-close? [1.0] state 1.0e-14))))
  (testing "An endpoint remainder smaller than min-step is allowed in both directions"
    (doseq [[start end] [[0.0 1.0] [1.0 0.0]]]
      (let [{:keys [status z state profile attempts]}
            (ode/integrate (constantly [2.0]) start end [(* 2.0 start)]
                           {:initial-step 0.375 :min-step 0.375 :max-step 0.375})
            penultimate (:z (nth profile (- (count profile) 2)))]
        (is (= :success status))
        (is (= end z))
        (is (= 3 attempts))
        (is (close? 0.25 (Math/abs (- z penultimate)) 1.0e-14))
        (is (state-close? [(* 2.0 end)] state 1.0e-14))))))

(deftest adaptive-minimum-step-failure-test
  (testing "A rejected step at min-step fails and retains the rejection diagnostics"
    (let [{:keys [status reason attempts profile last-rejection]}
          (ode/integrate exponential 0.0 1.0 [1.0]
                         {:initial-step 1.0 :min-step 1.0 :abs-tol 1.0e-8 :rel-tol 0.0})]
      (is (= :failed status))
      (is (= :step-too-small reason))
      (is (= 1 attempts))
      (is (= [{:z 0.0 :state [1.0]}] profile))
      (is (= :rejected (:status last-rejection)))
      (is (= :error-too-large (:reason last-rejection)))
      (is (> (:error last-rejection) 1.0))))
  (testing "A rejected final remainder keeps the profile accepted before the endpoint"
    (let [rhs (fn [z _] [(Math/pow (max 0.0 (- z 0.75)) 4.0)])
          {:keys [status reason z state attempts profile last-rejection]}
          (ode/integrate rhs 0.0 1.0 [0.0] {:initial-step 0.375 :min-step 0.375 :max-step 0.375
                                           :abs-tol 1.0e-8 :rel-tol 0.0})]
      (is (= :failed status))
      (is (= :step-too-small reason))
      (is (= 3 attempts))
      (is (= [0.0 0.375 0.75] (mapv :z profile)))
      (is (= {:z z :state state} (peek profile)))
      (is (state-close? [0.0] state 1.0e-14))
      (is (= :error-too-large (:reason last-rejection))))))

(deftest adaptive-coordinate-stagnation-test
  (let [start 1.0e16
        end (+ start 8.0)
        calls (atom 0)
        rhs (fn [_ _] (swap! calls inc) [1.0])
        result (ode/integrate rhs start end [0.0] {:initial-step 0.5})]
    (is (= {:status :failed :reason :coordinate-stagnation :attempts 0
            :z start :state [0.0] :profile [{:z start :state [0.0]}]}
           result))
    (is (zero? @calls))))

(deftest integration-input-validation-test
  (let [calls (atom 0)
        rhs (fn [_ _] (swap! calls inc) [1.0])]
    (testing "Invalid arguments fail before evaluating the RHS"
      (doseq [[f start end y options] [[nil 0.0 1.0 [1.0] {}]
                                       [rhs "0" 1.0 [1.0] {}]
                                       [rhs ##NaN 1.0 [1.0] {}]
                                       [rhs 0.0 ##Inf [1.0] {}]
                                       [rhs 0.0 1.0 '(1.0) {}]
                                       [rhs 0.0 1.0 [] {}]
                                       [rhs 0.0 1.0 [1.0] nil]]]
        (is (thrown-with-msg? clojure.lang.ExceptionInfo #"Invalid integration arguments"
                              (ode/integrate f start end y options)))))
    (testing "Zero-length and overflowing spans are rejected"
      (doseq [[start end] [[0.0 0.0] [-1.0e308 1.0e308]]]
        (is (thrown-with-msg? clojure.lang.ExceptionInfo #"Invalid integration interval"
                              (ode/integrate rhs start end [1.0])))))
    (testing "Every controller and integration option is validated"
      (doseq [options [{:abs-tol 0.0} {:abs-tol ##NaN} {:rel-tol -1.0} {:rel-tol ##Inf}
                       {:safety 0.0} {:safety 1.0} {:min-factor 0.0} {:min-factor 1.0}
                       {:max-factor 0.5} {:initial-step ##Inf} {:initial-step -1.0}
                       {:min-step 0.0} {:min-step 0.4 :initial-step 0.2}
                       {:max-step 0.2 :initial-step 0.4} {:max-attempts 0} {:max-attempts 1.5}
                       {:valid-state? nil}]]
        (is (thrown-with-msg? clojure.lang.ExceptionInfo #"Invalid integration options"
                              (ode/integrate rhs 0.0 1.0 [1.0] options))
            (str options))))
    (testing "Numerical and problem-specific initial state validation precede every attempt"
      (let [options {:valid-state? (fn [[y]] (pos? y))}]
        (doseq [initial [[##NaN] [##Inf] ["1.0"] [-1.0]]]
          (is (= :invalid-state
                 (:type (exception-info-data #(ode/integrate rhs 0.0 1.0 initial options))))))))
    (is (zero? @calls))))

(deftest rk4-finite-increment-overflow-regression-test
  ;; Integrating the finite constant slope over this short interval adds one.
  (is (state-close? [1.0] (ode/rk4-step (constantly [1.0e308]) 0.0 [0.0] 1.0e-308)
                    1.0e-14)))

(deftest richardson-finite-correction-overflow-regression-test
  ;; In relative coordinates the quartic takes [1 -2 1 -2 1] at the
  ;; quadrature nodes. Simpson's full and composite steps are +h and -h.
  (let [h 1.0e308
        rhs (fn [z _]
              (let [u (/ z h)
                    w (* u (- 1.0 u))]
                [(- 1.0 (* 256.0 w (- 0.25 w)))]))
        {:keys [state error]} (ode/rk4-step-doubling rhs 0.0 [0.0] h)]
    (is (state-close? [(- h)] state 1.0e294))
    (is (state-close? [(* (- (/ 2.0 15.0)) h)] error 1.0e294))))

(deftest error-scale-overflow-regression-test
  (testing "The normalized error is one half even when the scale exceeds Double/MAX_VALUE"
    (doseq [[abs-tol rel-tol] [[1.0e308 1.0] [1.0 2.0]]]
      (is (close? 0.5 (ode/normalized-error [1.0e308] [1.0e308] [1.0e308] abs-tol rel-tol)
                  1.0e-14)))))

(deftest adaptive-coordinate-rounding-regression-test
  ;; A requested step of 3 rounds to a displacement of 4 at this coordinate.
  ;; The solution of y' = 1 must still agree with the returned coordinate.
  (let [start 1.0e16
        end (+ start 8.0)
        {:keys [status z state]} (ode/integrate (constantly [1.0]) start end [0.0]
                                               {:initial-step 3.0})]
    (is (= :success status))
    (is (= end z))
    (is (state-close? [(- end start)] state 1.0e-14))))
