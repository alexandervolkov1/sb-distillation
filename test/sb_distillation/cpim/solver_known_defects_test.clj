(ns sb-distillation.cpim.solver-known-defects-test
  "Expected-behavior regressions for open solver defects. Run with -M:test --known-defects."
  (:require [clojure.test :refer [deftest is testing]]
            [sb-distillation.cpim.solver :as solver]
            [sb-distillation.cpim.solver-test :refer [conditions]]
            [sb-distillation.cpim.transport :as transport]))

(deftest negative-pressure-below-double-range-test
  (is (not (solver/valid-pressure-state? [-1e-500M 1.0]))))

(deftest invalid-physical-constants-test
  (doseq [[key value] [[:gas-const 0.0]
                       [:gas-const "8.314462618"]
                       [:sb4-collision-diameter -0.600e-9]
                       [:h2-reference-viscosity ##NaN]
                       [:h2-reference-viscosity ##Inf]]]
    (testing (str key " = " (pr-str value))
      (is (thrown? clojure.lang.ExceptionInfo
                   (solver/prepare-problem (assoc transport/default-parameters key value)
                                           conditions)))))
  (is (thrown? clojure.lang.ExceptionInfo
               (solver/prepare-problem (dissoc transport/default-parameters :gas-const)
                                       conditions))))

(deftest unrepresentable-derived-quantities-test
  (doseq [[label overrides]
          [["Area overflow" {:diameter 1e200}]
           ["Area underflow" {:diameter 1e-200}]
           ["Throughput overflow" {:ref-volumetric-flow Double/MAX_VALUE}]
           ["Antoine underflow" {:temp 153.0}]]]
    (testing label
      (is (thrown? clojure.lang.ExceptionInfo
                   (solver/prepare-problem transport/default-parameters
                                           (merge conditions overrides))))))
  (testing "Two valid numeric input pressures sum to a finite double"
    (is (== (* 2.0 Long/MAX_VALUE)
            (second (:outlet-state
                     (solver/prepare-problem transport/default-parameters
                                             (assoc conditions
                                                    :gauge-pressure Long/MAX_VALUE
                                                    :trap-pressure-drop Long/MAX_VALUE))))))))
