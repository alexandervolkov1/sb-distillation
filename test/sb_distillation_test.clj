
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

(deftest liquid-height-from-mass-test
  (let [density 6400.0
        total-volume (model/amp-volume radius straight-length)
        height (fn [mass]
                 (model/liquid-height-from-mass
                  radius straight-length mass density))
        mass 0.5
        liquid-height (height mass)]

    (is (close? (height 0.0) 0.0))
    (is (close? (height (* density (/ total-volume 2.0)))
                radius))
    (is (close? (height (* density total-volume))
                (* 2.0 radius)))

    (is (close? (model/filled-volume
                 radius straight-length liquid-height)
                (/ mass density)))

    (is (thrown? clojure.lang.ExceptionInfo
                 (height -0.1)))

    (is (thrown? clojure.lang.ExceptionInfo
                 (height (* 2.0 density total-volume))))

    (is (thrown? clojure.lang.ExceptionInfo
                 (model/liquid-height-from-mass
                  radius straight-length mass 0.0)))))

(deftest vapor-pressure-test
  (let [liquid-pressure
        (fn [temp-c]
          (model/pa-to-torr
           (model/sb-liquid-vapor-pressure temp-c)))

        solid-pressure
        (fn [temp-c]
          (model/pa-to-torr
           (model/sb-solid-vapor-pressure temp-c)))]

    ;; Pressure increases with temperature.
    (is (< (liquid-pressure 650.0)
           (liquid-pressure 700.0)
           (liquid-pressure 750.0)))

    (is (< (solid-pressure 450.0)
           (solid-pressure 500.0)
           (solid-pressure 550.0)))

    ;; Reference values calculated from the correlations.
    (is (< (Math/abs (- (liquid-pressure 700.0)
                        0.48189078))
           1.0e-6))

    (is (< (Math/abs (- (solid-pressure 500.0)
                        0.001773684))
           1.0e-7))))

(defn -main
  [& _]
  (let [{:keys [fail error]} (run-tests 'sb-distillation-test)]
    (when (pos? (+ fail error))
      (System/exit 1))))
