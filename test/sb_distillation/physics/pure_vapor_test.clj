(ns sb-distillation.physics.pure-vapor-test
  (:require [clojure.test :refer [deftest is]]
            [sb-distillation.cpim.model :as cpim]
            [sb-distillation.dgm.model :as dgm]
            [sb-distillation.physics.common :as physics]
            [sb-distillation.physics.pure-vapor :as pure]
            [sb-distillation.physics.thermodynamics :as thermo]
            [sb-distillation.physics.transport :as transport]
            [sb-distillation.solver :as solver]))

(defn- relative-error [expected actual]
  (abs (- (/ actual expected) 1.0)))

(deftest pure-vapor-reference-test
  ;; SI values evaluated with 60-digit Decimal arithmetic, including the
  ;; Antoine, hard-sphere viscosity and Knudsen expressions.
  (doseq [[temp diameter flux mass-flow days]
          [[953.15 0.006 0.011018625238525551 1.5173453047467074e-7 38.13922262080664]
           [973.15 0.010 0.042114686606350275 1.6109723460578892e-6 3.592263424755948]
           [993.15 0.012 0.09669126426121122 5.326037782782578e-6 1.0865557611597736]]]
    (let [params transport/default-parameters
          actual (pure/molar-flux params temp diameter 0.5)
          flow (solver/sb-mass-flow actual (solver/cross-sectional-area diameter)
                                   (:sb4-molar-mass params))]
      (is (< (relative-error flux actual) 1e-13))
      (is (< (relative-error mass-flow flow) 1e-13))
      (is (< (relative-error days (/ (solver/distillation-time 0.5 flow) 86400.0)) 1e-13))
      (is (< (relative-error (/ flux 2.0) (pure/molar-flux params temp diameter 1.0)) 1e-13)))))

(deftest pure-vapor-limits-and-gradients-test
  (let [params transport/default-parameters
        temp 973.15
        p-sat (thermo/sb-liquid-vapor-pressure temp)]
    (let [small (pure/molar-flux params temp 1e-10 0.5)
          large (pure/molar-flux params temp 1e5 0.5)]
      ;; Knudsen flux grows as d, whereas the viscous term grows as d^2.
      (is (< (relative-error (* 2.0 small) (pure/molar-flux params temp 2e-10 0.5)) 1e-7))
      (is (< (relative-error (* 4.0 large) (pure/molar-flux params temp 2e5 0.5)) 1e-7)))
    (let [diameter 0.01
          n1 (pure/molar-flux params temp diameter 0.5)
          {:keys [rt dk1 mu1] :as tube} (physics/prepare-coefficients params temp diameter)]
      (doseq [p1 [1.0 (/ p-sat 2.0) p-sat 1e4]]
        (let [expected (/ (* (- rt) n1) (+ dk1 (/ (* diameter diameter p1) (* 32.0 mu1))))
              [cp1 cp2] (cpim/pressure-derivatives tube p1 0.0 n1 0.0)
              [dg1 dg2] (dgm/pressure-derivatives (dgm/prepare-coefficients tube) p1 0.0 n1 0.0)]
          (is (< (relative-error expected cp1) 1e-11))
          (is (< (relative-error expected dg1) 1e-11))
          (is (zero? cp2))
          (is (zero? dg2)))))))
