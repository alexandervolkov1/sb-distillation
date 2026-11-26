(ns sb-distillation.cpim.transport-test
  (:require [clojure.test :refer [deftest is testing]]
            [sb-distillation.cpim.transport :as transport]))

(def ^:private parameters transport/default-parameters)
(def ^:private hot-temperature 973.15)

;; Reference values were evaluated from the pinned Typst equations and SI
;; parameters with 50-digit Python Decimal arithmetic, then rounded to 16
;; significant digits. The viscosity and Wilke fixtures below use 973.15 K.
(def ^:private reference-viscosities [5.105179111305377e-5 1.848032814750365e-5])
(def ^:private reference-wilke-factors [0.04587204907763955 4.011869427733789])

(defn- relatively-close?
  [expected actual]
  ;; 1e-12 relative allows for fixture rounding and chains of double operations
  ;; across quantities with different units and magnitudes.
  (<= (Math/abs (- actual expected)) (* 1.0e-12 (Math/abs expected))))

(deftest knudsen-diffusivity-test
  (let [diameter 0.012
        {:keys [sb4-molar-mass h2-molar-mass]} parameters]
    (testing "Both species in a 12 mm tube at 973.15 K"
      (is (relatively-close? 0.8227250469449949
                             (transport/knudsen-diffusivity parameters hot-temperature
                                                            diameter sb4-molar-mass)))
      (is (relatively-close? 12.78806091762617
                             (transport/knudsen-diffusivity parameters hot-temperature
                                                            diameter h2-molar-mass))))
    (let [temperature 300.0
          base (transport/knudsen-diffusivity parameters temperature diameter sb4-molar-mass)]
      (testing "Linear dependence on tube diameter"
        (is (relatively-close? (* 2.0 base)
                               (transport/knudsen-diffusivity parameters temperature
                                                              (* 2.0 diameter) sb4-molar-mass))))
      (testing "Square-root dependence on absolute temperature"
        (is (relatively-close? (* 2.0 base)
                               (transport/knudsen-diffusivity parameters (* 4.0 temperature)
                                                              diameter sb4-molar-mass))))
      (testing "Inverse square-root dependence on molar mass"
        (is (relatively-close? (/ base 2.0)
                               (transport/knudsen-diffusivity parameters temperature
                                                              diameter (* 4.0 sb4-molar-mass))))))))

(deftest binary-diffusivity-factor-test
  (testing "Hard-sphere diffusion factor in Pa*m^2/s"
    (is (relatively-close? 20.42360712282102
                           (transport/binary-diffusivity-factor parameters hot-temperature))))
  (testing "Temperature scaling with exponent 3/2"
    (let [base (transport/binary-diffusivity-factor parameters 300.0)]
      (is (relatively-close? (* 8.0 base)
                             (transport/binary-diffusivity-factor parameters 1200.0))))))

(deftest binary-diffusivity-test
  (let [factor 20.4
        pressure 50.0
        base (transport/binary-diffusivity factor pressure)]
    (testing "Inverse dependence on total pressure"
      (is (relatively-close? (/ base 2.0)
                             (transport/binary-diffusivity factor (* 2.0 pressure)))))
    (testing "The product p*D_12 equals K_12"
      (doseq [p [pressure 100.0 250.0]]
        (is (relatively-close? factor (* p (transport/binary-diffusivity factor p)))
            (str "Pressure " p " Pa"))))))

(deftest sb4-viscosity-test
  (let [temperature 300.0
        base (transport/sb4-viscosity parameters temperature)]
    (testing "Hard-sphere viscosity in Pa*s"
      (is (relatively-close? 2.834534305222926e-5 base))
      (is (relatively-close? (first reference-viscosities)
                             (transport/sb4-viscosity parameters hot-temperature))))
    (testing "Square-root dependence on absolute temperature"
      (is (relatively-close? (* 2.0 base)
                             (transport/sb4-viscosity parameters (* 4.0 temperature)))))
    (testing "Inverse-square dependence on collision diameter"
      (let [larger-sb4 (update parameters :sb4-collision-diameter * 2.0)]
        (is (relatively-close? (/ base 4.0)
                               (transport/sb4-viscosity larger-sb4 temperature)))))))

(deftest h2-viscosity-test
  (testing "Sutherland viscosity in Pa*s"
    (is (relatively-close? 8.887048804384605e-6 (transport/h2-viscosity parameters 300.0)))
    (is (relatively-close? (second reference-viscosities)
                           (transport/h2-viscosity parameters hot-temperature))))
  (testing "Normalization at the reference temperature to within one ULP"
    (let [{:keys [h2-reference-temp h2-reference-viscosity]} parameters
          actual (transport/h2-viscosity parameters h2-reference-temp)]
      (is (<= (Math/abs (- actual h2-reference-viscosity))
              (Math/ulp (double h2-reference-viscosity)))))))

(deftest wilke-interaction-factor-test
  (let [[mu1 mu2] reference-viscosities
        [expected-phi12 expected-phi21] reference-wilke-factors
        {:keys [sb4-molar-mass h2-molar-mass]} parameters
        phi12 (transport/wilke-interaction-factor mu1 mu2 sb4-molar-mass h2-molar-mass)
        phi21 (transport/wilke-interaction-factor mu2 mu1 h2-molar-mass sb4-molar-mass)]
    (testing "Ordered Sb4/H2 factors are distinct and match the reference values"
      (is (relatively-close? expected-phi12 phi12))
      (is (relatively-close? expected-phi21 phi21))
      (is (< phi12 phi21)))
    (testing "Diagonal interaction factors equal one"
      (is (relatively-close? 1.0 (transport/wilke-interaction-factor mu1 mu1 sb4-molar-mass sb4-molar-mass)))
      (is (relatively-close? 1.0 (transport/wilke-interaction-factor mu2 mu2 h2-molar-mass h2-molar-mass))))))

(deftest mixture-viscosity-test
  (let [[mu1 mu2] reference-viscosities
        [phi12 phi21] reference-wilke-factors]
    (testing "Reference mixture viscosities and invariance under component exchange"
      (doseq [[x1 x2 expected] [[0.1 0.9 4.891634986306238e-5]
                                [0.5 0.5 5.249996699646302e-5]
                                [0.9 0.1 5.129093547843548e-5]]]
        (let [viscosity (transport/mixture-viscosity x1 x2 mu1 mu2 phi12 phi21)
              swapped (transport/mixture-viscosity x2 x1 mu2 mu1 phi21 phi12)]
          (is (relatively-close? expected viscosity) (str "Sb4 mole fraction " x1))
          (is (relatively-close? viscosity swapped) (str "Sb4 mole fraction " x1)))))
    (testing "Pure-component limits"
      (is (relatively-close? mu1 (transport/mixture-viscosity 1.0 0.0 mu1 mu2 phi12 phi21)))
      (is (relatively-close? mu2 (transport/mixture-viscosity 0.0 1.0 mu1 mu2 phi12 phi21))))
    (testing "Identical species retain their common viscosity at any composition"
      (is (relatively-close? mu1 (transport/mixture-viscosity 0.3 0.7 mu1 mu1 1.0 1.0))))))
