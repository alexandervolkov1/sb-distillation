(ns sb-distillation.numerics.root
  (:import [java.math BigDecimal]))

(defn- finite-number?
  "Return true if x is a finite number."
  [x]
  (and (number? x)
       (or (not (float? x)) (Double/isFinite (double x)))))

(defn- magnitude-fraction
  [x]
  (cond
    (ratio? x) [(.abs ^BigDecimal (bigdec (numerator x))) (bigdec (denominator x))]
    ;; The constructor preserves the exact binary float value.
    (float? x) [(.abs (BigDecimal. (double x))) 1M]
    :else [(.abs ^BigDecimal (bigdec x)) 1M]))

(defn- magnitude-at-most?
  [a b]
  (if (and (float? a) (float? b))
    (<= (Math/abs (double a)) (Math/abs (double b)))
    (let [[^BigDecimal a-num ^BigDecimal a-den] (magnitude-fraction a)
          [^BigDecimal b-num ^BigDecimal b-den] (magnitude-fraction b)]
      ;; Compare cross products exactly, including nonterminating ratios.
      (not (pos? (.compareTo (.multiply a-num b-den) (.multiply b-num a-den)))))))

(defn- opposite-signs?
  "Return true if a and b have strictly opposite signs."
  [a b]
  (or (and (neg? a) (pos? b))
      (and (pos? a) (neg? b))))

(defn- midpoint
  "Compute the midpoint of a finite interval
  without intermediate overflow."
  [a b]
  (let [width (- b a)]
    (if (Double/isFinite (double width))
      (+ a (* 0.5 width))
      (* 0.5 (+ a b)))))

(defn bisect
  "Find a root of a continuous scalar function using bisection.
  a and b are finite ordered interval endpoints evaluated as doubles.
  f-tol is the absolute residual tolerance in function units.
  max-iterations limits the number of midpoint evaluations.
  Success requires abs(f(x)) <= f-tol. Residuals retain their original numeric type."
  [f a b f-tol max-iterations]

  (when-not
   (and (ifn? f)
        (every? #(and (number? %) (Double/isFinite (double %))) [a b])
        (finite-number? f-tol)
        (< (double a) (double b))
        (not (neg? f-tol))
        (integer? max-iterations)
        (pos? max-iterations))
    (throw (ex-info "Invalid bisection arguments"
                    {:a a
                     :b b
                     :f-tol f-tol
                     :max-iterations max-iterations})))

  (letfn [(evaluate [x]
            (let [fx (f x)]
              (when-not (number? fx)
                (throw (ex-info "Invalid function value"
                                {:x x :fx fx})))
              fx))

          (result [status reason x fx a b iterations]
            (cond-> {:status status
                     :x x
                     :fx fx
                     :bracket [a b]
                     :iterations iterations}
              reason (assoc :reason reason)))]

    (let [a (double a)
          b (double b)]
      (loop [a a
             b b
             fa (evaluate a)
             fb (evaluate b)
             iterations 0]

        (cond
          (not (finite-number? fa))
          (result :failed :non-finite-value
                  a fa a b iterations)

          (not (finite-number? fb))
          (result :failed :non-finite-value
                  b fb a b iterations)

          (magnitude-at-most? fa f-tol)
          (result :success nil a fa a b iterations)

          (magnitude-at-most? fb f-tol)
          (result :success nil b fb a b iterations)

          :else
          (let [[best-x best-fx]
                (if (magnitude-at-most? fa fb)
                  [a fa]
                  [b fb])]

            (cond
              (not (opposite-signs? fa fb))
              (result :failed :not-bracketed
                      best-x best-fx a b iterations)

              (>= iterations max-iterations)
              (result :failed :max-iterations
                      best-x best-fx a b iterations)

              :else
              (let [c (double (midpoint a b))]

                (if (or (= c a) (= c b))
                  (result :failed :midpoint-stagnation
                          best-x best-fx a b iterations)

                  (let [fc (evaluate c)
                        iterations (inc iterations)]

                    (cond
                      (not (finite-number? fc))
                      (result :failed :non-finite-value
                              c fc a b iterations)

                      (magnitude-at-most? fc f-tol)
                      (result :success nil
                              c fc a b iterations)

                      (opposite-signs? fa fc)
                      (recur a c fa fc iterations)

                      :else
                      (recur c b fc fb iterations))))))))))))
