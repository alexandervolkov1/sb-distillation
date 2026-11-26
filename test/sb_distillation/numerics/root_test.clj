(ns sb-distillation.numerics.root-test
  (:require [clojure.test :refer [deftest is testing]]
            [sb-distillation.numerics.root :as root]))

(defn- close?
  [expected actual tolerance]
  (and (number? actual)
       (<= (Math/abs (- (double expected) (double actual))) tolerance)))

(defn- observed-bisection
  [f a b tolerance limit]
  (let [calls (atom [])
        result (root/bisect (fn [x] (swap! calls conj x) (f x)) a b tolerance limit)]
    {:result result :calls @calls}))

(defn- exception-data
  [call]
  (try
    (call)
    nil
    (catch clojure.lang.ExceptionInfo ex
      (ex-data ex))))

(deftest bisection-analytical-solutions-test
  ;; The cosine fixed point was calculated independently with Newton's method
  ;; and Taylor series for sine and cosine in 75-digit decimal arithmetic.
  (doseq [[label f a b expected]
          [["Increasing linear" #(- % 2.25) 0.0 5.0 2.25]
           ["Decreasing on the negative branch of a quadratic" #(- (* % %) 2.0)
            -2.0 -1.0 (- (Math/sqrt 2.0))]
           ["Increasing cubic" #(- (* % % %) 2.0) 0.0 2.0 (Math/cbrt 2.0)]
           ["Decreasing exponential" #(- 2.0 (Math/exp %)) 0.0 2.0 (Math/log 2.0)]
           ["Cosine fixed point" #(- (Math/cos %) %) 0.0 1.0 0.739085133215160641655312087674]
           ["Zero root" identity -0.25 1.0 0.0]]]
    (testing label
      (let [tolerance 1.0e-12
            {:keys [result calls]} (observed-bisection f a b tolerance 60)
            {:keys [status x fx bracket iterations]} result
            [left right] bracket]
        (is (= :success status))
        (is (nil? (:reason result)))
        (is (close? expected x 2.0e-12))
        (is (<= a left x right b))
        (is (= (f x) fx))
        (is (<= (Math/abs (double fx)) tolerance))
        (is (<= 1 iterations 60))
        (is (= (+ 2 iterations) (count calls) (count (distinct calls))))))))

(deftest bisection-bracket-and-evaluation-count-test
  (let [f #(- % 0.3)
        midpoints [0.5 0.25 0.375 0.3125 0.28125]
        brackets [[0.0 0.5] [0.25 0.5] [0.25 0.375]
                  [0.25 0.3125] [0.28125 0.3125]]]
    (doseq [[index expected-bracket] (map-indexed vector brackets)]
      (let [n (inc index)
            {:keys [result calls]} (observed-bisection f 0.0 1.0 0.0 n)
            {:keys [status reason x fx bracket iterations]} result
            [left right] bracket]
        (testing (str "After " n " midpoint evaluations")
          (is (= :failed status))
          (is (= :max-iterations reason))
          (is (= n iterations))
          (is (= expected-bracket bracket))
          (is (= (Math/scalb 1.0 (int (- n))) (- right left)))
          (is (< (f left) 0.0 (f right)))
          (is (<= left x right))
          (is (= (f x) fx))
          (is (= (+ 2 n) (count calls) (count (distinct calls))))
          (is (= (into [0.0 1.0] (take n midpoints)) calls)))))))

(deftest endpoint-convergence-test
  (doseq [[label f tolerance expected-x expected-fx]
          [["Exact left endpoint" identity 0.0 0.0 0.0]
           ["Exact right endpoint" #(- % 1.0) 0.0 1.0 0.0]
           ["Left residual equals tolerance" #(- % 0.125) 0.125 0.0 -0.125]
           ["Right residual equals tolerance" #(- % 0.875) 0.125 1.0 0.125]]]
    (testing label
      (let [{:keys [result calls]} (observed-bisection f 0.0 1.0 tolerance 2)]
        (is (= {:status :success :x expected-x :fx expected-fx
                :bracket [0.0 1.0] :iterations 0} result))
        (is (= [0.0 1.0] calls))))))

(deftest zero-tolerance-and-final-iteration-test
  (testing "An exact midpoint satisfies zero tolerance"
    (let [{:keys [result calls]} (observed-bisection identity -1.0 1.0 0.0 1)]
      (is (= {:status :success :x 0.0 :fx 0.0 :bracket [-1.0 1.0] :iterations 1} result))
      (is (= [-1.0 1.0 0.0] calls))))
  (testing "The last permitted midpoint can still succeed"
    (let [{:keys [result calls]} (observed-bisection #(- % 0.375) 0.0 1.0 0.0 3)]
      (is (= {:status :success :x 0.375 :fx 0.0 :bracket [0.25 0.5] :iterations 3} result))
      (is (= [0.0 1.0 0.5 0.25 0.375] calls)))))

(deftest unsuccessful-finite-brackets-test
  (testing "Same-sign endpoint values return the endpoint with the smaller residual"
    (doseq [[f expected] [[#(+ % 1.0) 0.0] [#(- 2.0 %) 1.0] [#(- % 2.0) 1.0]]]
      (let [{:keys [result calls]} (observed-bisection f 0.0 1.0 0.0 3)]
        (is (= {:status :failed :reason :not-bracketed :x expected :fx (f expected)
                :bracket [0.0 1.0] :iterations 0} result))
        (is (= [0.0 1.0] calls)))))
  (testing "An exhausted search returns the better endpoint of the updated bracket"
    (doseq [[zero-x expected] [[0.3 0.25] [0.35 0.375]]]
      (let [f #(- % zero-x)
            {:keys [result calls]} (observed-bisection f 0.0 1.0 0.0 3)]
        (is (= {:status :failed :reason :max-iterations :x expected :fx (f expected)
                :bracket [0.25 0.375] :iterations 3} result))
        (is (= [0.0 1.0 0.5 0.25 0.375] calls))))))

(deftest midpoint-with-extreme-endpoints-test
  (let [unit (Math/scalb 1.0 (int 1022))]
    (doseq [[label a b expected]
            [["Positive sum would overflow" (* 2.0 unit) (* 3.0 unit) (* 2.5 unit)]
             ["Negative sum would overflow" (* -3.0 unit) (* -2.0 unit) (* -2.5 unit)]
             ["Width would overflow" (* -3.0 unit) (* 2.0 unit) (* -0.5 unit)]
             ["Full symmetric double range" (- Double/MAX_VALUE) Double/MAX_VALUE 0.0]]]
      (testing label
        (let [f (fn [x] (- (/ x unit) (/ expected unit)))
              {:keys [result calls]} (observed-bisection f a b 0.0 1)]
          (is (= {:status :success :x expected :fx 0.0 :bracket [a b] :iterations 1} result))
          (is (= [a b expected] calls)))))))

(deftest midpoint-and-residuals-in-the-subnormal-range-test
  (let [unit Double/MIN_VALUE]
    (doseq [[a b expected] [[(* -4.0 unit) (* 6.0 unit) unit]
                            [unit (* 5.0 unit) (* 3.0 unit)]]]
      (let [f #(- % expected)
            {:keys [result calls]} (observed-bisection f a b 0.0 1)]
        (is (= {:status :success :x expected :fx 0.0 :bracket [a b] :iterations 1} result))
        (is (= [a b expected] calls))))))

(deftest adjacent-endpoint-stagnation-test
  (doseq [[a b scale root-offset expected]
          [[1.0 (Math/nextUp 1.0) (Math/ulp 1.0) 0.75 (Math/nextUp 1.0)]
           [Double/MIN_VALUE (* 2.0 Double/MIN_VALUE) Double/MIN_VALUE 0.25 Double/MIN_VALUE]]]
    (let [f (fn [x] (- (/ (- x a) scale) root-offset))
          {:keys [result calls]} (observed-bisection f a b 0.0 5)]
      (is (= (Math/nextUp (double a)) b))
      (is (= {:status :failed :reason :midpoint-stagnation
              :x expected :fx (f expected) :bracket [a b] :iterations 0} result))
      (is (= [a b] calls)))))

(deftest sign-check-with-unrepresentable-products-test
  (doseq [scale [1.0e308 1.0e-308]
          slope [scale (- scale)]]
    (let [f #(* slope %)
          {:keys [result calls]} (observed-bisection f -1.0 1.0 0.0 1)]
      (is (= :success (:status result)))
      (is (zero? (:fx result)))
      (is (= [-1.0 1.0 0.0] calls)))))

(deftest non-finite-function-values-test
  (doseq [value [##NaN ##Inf ##-Inf]
          [point iterations expected-calls] [[-1.0 0 [-1.0 1.0]]
                                             [1.0 0 [-1.0 1.0]]
                                             [0.0 1 [-1.0 1.0 0.0]]]]
    (testing (str "Non-finite value " value " at " point)
      (let [f (fn [x] (if (= point x) value x))
            {:keys [result calls]} (observed-bisection f -1.0 1.0 0.0 4)]
        (is (= {:status :failed :reason :non-finite-value :x point
                :bracket [-1.0 1.0] :iterations iterations}
               (dissoc result :fx)))
        (is (= (Double/doubleToLongBits (double value))
               (Double/doubleToLongBits (double (:fx result)))))
        (is (= expected-calls calls))))))

(deftest invalid-function-values-and-exceptions-test
  (testing "Nonnumeric function values report the offending point and value"
    (doseq [value [nil "invalid" [1.0]]
            point [-1.0 0.0 1.0]]
      (let [f (fn [x] (if (= point x) value x))]
        (is (= {:x point :fx value}
               (exception-data #(root/bisect f -1.0 1.0 0.0 4)))))))
  (testing "User exceptions propagate unchanged"
    (doseq [ex [(ex-info "Function evaluation failed" {:source :user})
                (IllegalStateException. "Unexpected exception")]
            point [-1.0 0.0 1.0]]
      (let [f (fn [x] (if (= point x) (throw ex) x))]
        (is (identical? ex (try (root/bisect f -1.0 1.0 0.0 4)
                               nil
                               (catch Exception actual actual))))))))

(deftest invalid-bisection-arguments-test
  (let [calls (atom 0)
        f (fn [x] (swap! calls inc) x)]
    (doseq [[a b tolerance limit]
            [[1.0 -1.0 0.0 1] [1.0 1.0 0.0 1] ["-1" 1.0 0.0 1]
             [##NaN 1.0 0.0 1] [##-Inf 1.0 0.0 1] [-1.0 ##Inf 0.0 1]
             [9007199254740992N 9007199254740993N 0.0 1]
             [-1.0 1.0 -1.0 1] [-1.0 1.0 ##NaN 1] [-1.0 1.0 ##Inf 1]
             [-1.0 1.0 ##-Inf 1] [-1.0 1.0 "0" 1]
             [-1.0 1.0 0.0 0] [-1.0 1.0 0.0 -1] [-1.0 1.0 0.0 1.5]
             [-1.0 1.0 0.0 ##Inf] [-1.0 1.0 0.0 nil]]]
      (is (thrown-with-msg? clojure.lang.ExceptionInfo #"Invalid bisection arguments"
                            (root/bisect f a b tolerance limit))
          (str [a b tolerance limit])))
    (is (thrown-with-msg? clojure.lang.ExceptionInfo #"Invalid bisection arguments"
                          (root/bisect nil -1.0 1.0 0.0 1)))
    (is (zero? @calls))))

(deftest residual-conversion-underflow-regression-test
  (testing "A nonzero endpoint residual below the double range cannot meet zero tolerance"
    (let [f (constantly 1.0e-400M)
          {:keys [status reason x fx]} (root/bisect f -1.0 1.0 0.0 1)]
      (is (= [:failed :not-bracketed] [status reason]))
      (is (= (f x) fx))))
  (testing "A nonzero midpoint residual below the double range still refines the bracket"
    (let [f (fn [x] (+ (bigdec x) 1.0e-400M))
          {:keys [status reason x fx bracket iterations]} (root/bisect f -1.0 1.0 0.0 1)]
      (is (= [:failed :max-iterations] [status reason]))
      (is (= [-1.0 0.0] bracket))
      (is (= 1 iterations))
      (is (= (f x) fx)))))

(deftest mixed-precision-residual-comparison-regression-test
  (doseq [[fx tolerance] [[1.00000000000000001M 1.0] [1.0 0.99999999999999999M]]]
    (let [{:keys [status reason] :as result} (root/bisect (constantly fx) -1.0 1.0 tolerance 1)]
      (is (= [:failed :not-bracketed] [status reason]))
      (is (= fx (:fx result))))))

(deftest negative-subnormal-tolerance-regression-test
  (is (thrown-with-msg? clojure.lang.ExceptionInfo #"Invalid bisection arguments"
                        (root/bisect (constantly 0.0) -1.0 1.0 -1.0e-400M 1))))

(deftest exact-decimal-residuals-outside-the-double-range-test
  (doseq [scale [1.0e400M 1.0e-400M]]
    (testing (str "Nonzero endpoint residuals at scale " scale)
      (doseq [fx [scale (- scale)]]
        (let [f (constantly fx)
              rejected (root/bisect f -1.0 1.0 0.0 1)
              accepted (root/bisect f -1.0 1.0 scale 1)]
          (is (= [:failed :not-bracketed] ((juxt :status :reason) rejected)))
          (is (identical? fx (:fx rejected)))
          (is (= :success (:status accepted)))
          (is (identical? fx (:fx accepted))))))
    (testing "Opposite endpoint signs allow a midpoint evaluation"
      (let [f (fn [x] (* scale (+ 1M (* 3M (bigdec x)))))
            {:keys [result calls]} (observed-bisection f -1.0 1.0 scale 1)]
        (is (= {:status :success :x 0.0 :fx scale :bracket [-1.0 1.0] :iterations 1}
               result))
        (is (instance? java.math.BigDecimal (:fx result)))
        (is (= [-1.0 1.0 0.0] calls))))))

(deftest mixed-numeric-magnitude-test
  (doseq [[label fx tolerance expected]
          [["An integer above the double precision limit" 9007199254740993N 9007199254740992.0 :failed]
           ["A large integer equals a decimal" 9007199254740993N 9007199254740993M :success]
           ["A ratio just below a large integer" 9007199254740993N 18014398509481985/2 :failed]
           ["A large integer just above a negative ratio's magnitude"
            -18014398509481985/2 9007199254740993N :success]
           ["The absolute value of the smallest long" Long/MIN_VALUE 9223372036854775808N :success]
           ["A terminating ratio equals a decimal" -3/4 0.75M :success]
           ["A decimal strictly below one third" -1/3 0.3333333333333333333333333333333333M :failed]
           ["One third is strictly above this decimal" 0.3333333333333333333333333333333333M 1/3 :success]
           ["Binary double 0.1 exceeds decimal 0.1" 0.1 0.1M :failed]
           ["Binary double 0.1 exceeds the ratio one tenth" 0.1 1/10 :failed]
           ["The ratio one tenth is smaller than binary double 0.1" 1/10 0.1 :success]
           ["Binary double 0.1 equals its full decimal expansion"
            0.1 0.1000000000000000055511151231257827021181583404541015625M :success]
           ["Binary float 0.1 exceeds binary double 0.1" (float 0.1) 0.1 :failed]
           ["Binary float 0.1 equals its exact ratio" (float -0.1) 13421773/134217728 :success]]]
    (testing label
      (let [{:keys [status reason iterations] :as result}
            (root/bisect (constantly fx) -1.0 1.0 tolerance 1)]
        (is (= expected status))
        (is (= (when (= expected :failed) :not-bracketed) reason))
        (is (zero? iterations))
        (is (identical? fx (:fx result)))))))

(deftest best-endpoint-with-mixed-precision-test
  (testing "Stagnation returns the boundary with the strictly smaller magnitude"
    (let [a 1.0
          b (Math/nextUp a)]
      (doseq [[fa fb] [[-1.00000000000000001M 1.0]
                       [-9007199254740993N 18014398509481985/2]]]
        (let [f (fn [x] (if (= x a) fa fb))
              {:keys [result calls]} (observed-bisection f a b 0 2)]
          (is (= {:status :failed :reason :midpoint-stagnation :x b :fx fb
                  :bracket [a b] :iterations 0} result))
          (is (identical? fb (:fx result)))
          (is (= [a b] calls))))))
  (testing "The best boundary is selected from the refined bracket after exhaustion"
    (let [f (fn [x] (cond (neg? x) -3.0 (zero? x) -1.00000000000000001M :else 1.0))
          {:keys [result calls]} (observed-bisection f -1.0 1.0 0.0 1)]
      (is (= {:status :failed :reason :max-iterations :x 1.0 :fx 1.0
              :bracket [0.0 1.0] :iterations 1} result))
      (is (= [-1.0 1.0 0.0] calls)))))

(deftest non-finite-float-residuals-and-tolerances-test
  (doseq [value [Float/NaN Float/POSITIVE_INFINITY Float/NEGATIVE_INFINITY]]
    (let [{:keys [status reason fx]} (root/bisect (constantly value) -1.0 1.0 0.0 1)]
      (is (= [:failed :non-finite-value] [status reason]))
      (is (identical? value fx)))
    (is (thrown-with-msg? clojure.lang.ExceptionInfo #"Invalid bisection arguments"
                          (root/bisect identity -1.0 1.0 value 1)))))
