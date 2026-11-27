(ns sb-distillation.cpim.thermodynamics-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :refer [deftest is testing]]
            [sb-distillation.cpim.thermodynamics :as thermo]
            [sb-distillation.reference.legacy :as legacy]))

(deftest vapor-pressure-in-pascals-test
  (let [reference (edn/read-string (slurp (io/resource "fixtures/cpim-reference.edn")))]
    (is (< (abs (- (:saturation reference) (thermo/sb-liquid-vapor-pressure 973.15))) 1e-12))
    (is (< (abs (- (:solid-pressure-at-773.15 reference)
                   (thermo/sb-solid-vapor-pressure 773.15))) 1e-13)))
  (is (apply < (map thermo/sb-liquid-vapor-pressure [953.15 973.15 993.15 1058.0 1373.0])))
  (is (apply < (map thermo/sb-solid-vapor-pressure [723.15 773.15 783.15 823.15])))
  (doseq [[current historical temperatures]
          [[thermo/sb-liquid-vapor-pressure legacy/sb-liquid-vapor-pressure [680 700 720 800 1050]]
           [thermo/sb-solid-vapor-pressure legacy/sb-solid-vapor-pressure [450 500 510 550]]]
          celsius temperatures]
    (testing (str "Legacy Celsius input at " celsius " degrees")
      (is (< (abs (- (/ (current (+ celsius 273.15)) (historical celsius)) 1.0)) 1e-13)))))
