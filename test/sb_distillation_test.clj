
(ns sb-distillation-test
  (:require [clojure.test :refer [deftest is run-tests]]
            [sb-distillation :as model]))

(def radius 0.03)
(def straight-length 0.18)
(def eps 1.0e-12)

(defn close?
  [a b]
  (< (Math/abs (- a b)) eps))

(deftest filled-volume-test
  (let [total-volume (model/amp-volume radius straight-length)
        volume (fn [liquid-height]
                 (model/filled-volume radius straight-length liquid-height))]
    (is (close? (volume 0.0) 0.0))
    (is (close? (volume radius) (/ total-volume 2.0)))
    (is (close? (volume (* 2.0 radius)) total-volume))))

(deftest liquid-surface-area-test
  (let [liquid-area (fn [liquid-height]
                      (model/liquid-surface-area radius
                                                 straight-length
                                                 liquid-height))
        max-area (+ (* 2.0 radius straight-length)
                    (* Math/PI radius radius))]
    (is (close? (liquid-area 0.0) 0.0))
    (is (close? (liquid-area radius) max-area))
    (is (close? (liquid-area (* 2.0 radius)) 0.0))))

(defn -main
  [& _]
  (let [{:keys [fail error]} (run-tests 'sb-distillation-test)]
    (when (pos? (+ fail error))
      (System/exit 1))))
