(ns sb-distillation.numerics.root)

(defn- finite-number?
  "Return true if x is a finite number."
  [x]
  (and (number? x)
       (Double/isFinite (double x))))

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
