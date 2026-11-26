(ns sb-distillation.cpim.model-test
  (:require [clojure.test :refer [deftest is testing]]
            [sb-distillation.cpim.model :as model]
            [sb-distillation.cpim.transport :as transport]))

(def ^:private temperature 973.15)
(def ^:private diameter 0.012)

;; The reference fixtures use 50-digit Python Decimal arithmetic on the
;; pinned Typst equations and SI inputs, rounded to 16 significant digits.
;; The calculation uses D_12 and A_C directly in the harmonic products.
(def ^:private reference-constants
  {:rt 8091.21929670670
   :diameter 0.012
   :ak 14092.13962560272
   :sqrt-m1 0.6978825116020604
   :sqrt-m2 0.04489855231519163
   :dk1 0.8227250469449949
   :dk2 12.78806091762617
   :k12 20.42360712282102
   :mu1 5.105179111305377e-5
   :mu2 1.848032814750365e-5
   :phi12 0.04587204907763955
   :phi21 4.011869427733789})

(def ^:private alternate-parameters
  (assoc transport/default-parameters
         :gas-const 8.2
         :boltzmann-const 1.42e-23
         :avogadro-const 6.1e23
         :sb4-molar-mass 0.420
         :h2-molar-mass 0.003
         :sb4-collision-diameter 0.68e-9
         :h2-collision-diameter 0.31e-9
         :h2-reference-viscosity 9.1e-6
         :h2-reference-temp 310.0
         :h2-sutherland-constant 90.0))

(def ^:private alternate-constants
  {:rt 5330.0
   :diameter 0.008
   :ak 17156.34674339199
   :sqrt-m1 0.6480740698407860
   :sqrt-m2 0.05477225575051661
   :dk1 0.4793776131942094
   :dk2 5.672072411777620
   :k12 7.747887674036033
   :mu1 3.039609814548452e-5
   :mu2 1.493472252342394e-5
   :phi12 0.05959366942657983
   :phi21 4.099279051838347})

(def ^:private local-reference-cases
  [[[12.0 48.0] {:p 60.0 :x1 0.2 :x2 0.8
                 :mu 5.236318833170473e-5
                 :d-eff1 0.2407753114217077 :d-eff2 0.3315677594762190
                 :inverse-ac 1.118379786443496e-4 :a 5470.473774421573}]
   [[180.0 20.0] {:p 200.0 :x1 0.9 :x2 0.1
                  :mu 5.129093547843548e-5
                  :d-eff1 0.09084250856055518 :d-eff2 0.1013090403910355
                  :inverse-ac 1.371848648658792e-3 :a 693.0918594866533}]
   [[0.0 60.0] {:p 60.0 :x1 0.0 :x2 1.0
                :mu 1.848032814750365e-5
                :d-eff1 0.2407753114217077 :d-eff2 0.3315677594762190
                :inverse-ac 8.107229209266536e-5 :a 6577.483139351885}]
   [[60.0 0.0] {:p 60.0 :x1 1.0 :x2 0.0
                :mu 5.105179111305377e-5
                :d-eff1 0.2407753114217077 :d-eff2 0.3315677594762190
                :inverse-ac 4.561641215332542e-4 :a 1897.080835775028}]])

(defn- close?
  ([expected actual]
   ;; Relative tolerance covers rounding of fixtures and chained double
   ;; operations across different SI scales. Exact zeros use zero tolerance.
   (close? expected actual 0.0 1.0e-12))
  ([expected actual absolute-tolerance relative-tolerance]
   (and (number? actual)
        (<= (Math/abs (- expected actual))
            (+ absolute-tolerance (* relative-tolerance (Math/abs expected)))))))

(defn- check-coefficient-map
  [expected actual]
  (is (= (set (keys expected)) (set (keys actual))))
  (doseq [[coefficient value] expected]
    (is (close? value (get actual coefficient)) (name coefficient))))

(defn- prepare-reference-tube
  []
  (model/prepare-coefficients transport/default-parameters temperature diameter))

(deftest molecular-resistance-test
  (let [{:keys [rt ak]} reference-constants
        base (model/molecular-resistance rt diameter)]
    (testing "Independent resistance values for two isothermal tubes"
      (is (close? ak base))
      (is (close? (:ak alternate-constants)
                  (model/molecular-resistance (:rt alternate-constants) 0.008))))
    (testing "Inverse diameter and square-root temperature scaling at fixed R"
      (is (close? (/ base 2.0) (model/molecular-resistance rt (* 2.0 diameter))))
      (is (close? (* 2.0 base) (model/molecular-resistance (* 4.0 rt) diameter))))))

(deftest prepare-coefficients-test
  (testing "Complete map for the default parameters"
    (check-coefficient-map reference-constants (prepare-reference-tube)))
  (testing "Complete map with changed physical constants, species and Sutherland parameters"
    (check-coefficient-map alternate-constants
                           (model/prepare-coefficients alternate-parameters 650.0 0.008))))

(deftest effective-diffusivity-test
  (let [{:keys [k12 dk1 dk2]} reference-constants]
    (testing "Numerical references for both species at different pressures"
      (doseq [[_ {:keys [p d-eff1 d-eff2]}] (take 2 local-reference-cases)]
        (is (close? d-eff1 (model/effective-diffusivity p k12 dk1)))
        (is (close? d-eff2 (model/effective-diffusivity p k12 dk2)))))
    (testing "Both coefficients decrease strictly with increasing pressure"
      (doseq [dk [dk1 dk2]]
        (is (apply > (map #(model/effective-diffusivity % k12 dk) [1.0 60.0 200.0])))))
    (testing "Limits and equal diffusion resistances for both species"
      (doseq [dk [dk1 dk2]]
        (let [crossover (/ k12 dk)
              ratio 1.0e-6
              low-pressure (* crossover ratio)
              high-pressure (/ crossover ratio)]
          ;; The neglected resistance is ratio times the retained one.
          (is (close? dk (model/effective-diffusivity low-pressure k12 dk) 0.0 ratio))
          (is (close? (/ k12 high-pressure)
                      (model/effective-diffusivity high-pressure k12 dk) 0.0 ratio))
          (is (close? (/ dk 2.0) (model/effective-diffusivity crossover k12 dk))))))))

(deftest inverse-continuum-resistance-test
  (let [{:keys [rt sqrt-m1 sqrt-m2]} reference-constants]
    (testing "Numerical references at different compositions and pressures"
      (doseq [[_ {:keys [p x1 x2 mu inverse-ac]}] local-reference-cases]
        (is (close? inverse-ac
                    (model/inverse-continuum-resistance p mu rt diameter x1 x2 sqrt-m1 sqrt-m2)))))
    (testing "Pressure, diameter, and viscosity scaling at fixed composition"
      (let [[_ {:keys [p x1 x2 mu inverse-ac]}] (first local-reference-cases)
            resistance (fn [pressure viscosity tube-diameter]
                         (model/inverse-continuum-resistance pressure viscosity rt tube-diameter
                                                             x1 x2 sqrt-m1 sqrt-m2))]
        (is (close? (* 2.0 inverse-ac) (resistance (* 2.0 p) mu diameter)))
        (is (close? (* 4.0 inverse-ac) (resistance p mu (* 2.0 diameter))))
        (is (close? (/ inverse-ac 2.0) (resistance p (* 2.0 mu) diameter)))
        (is (zero? (resistance 0.0 mu diameter)))))))

(deftest pressure-resistance-test
  (let [ak (:ak reference-constants)]
    (testing "Numerical references"
      (doseq [[_ {:keys [inverse-ac a]}] local-reference-cases]
        (is (close? a (model/pressure-resistance inverse-ac ak)))))
    (testing "Molecular and continuum limits and equal resistances"
      (let [ratio 1.0e-6
            low-inverse-ac (/ ratio ak)
            high-inverse-ac (/ 1.0 (* ratio ak))]
        (is (close? ak (model/pressure-resistance low-inverse-ac ak) 0.0 ratio))
        (is (close? (/ 1.0 high-inverse-ac)
                    (model/pressure-resistance high-inverse-ac ak) 0.0 ratio))
        (is (close? (/ ak 2.0) (model/pressure-resistance (/ 1.0 ak) ak)))))))

(deftest local-coefficients-test
  (let [coefficients (prepare-reference-tube)]
    (testing "Mixed states and both pure limits, including the H2-only outlet"
      (doseq [[[p1 p2] expected] local-reference-cases]
        (testing (str "[p_1 p_2] = " [p1 p2] " Pa")
          (let [{:keys [x1 x2] :as actual} (model/local-coefficients coefficients p1 p2)]
            (check-coefficient-map expected actual)
            ;; Mole fractions are dimensionless and normalized to within a few ULPs.
            (is (close? 1.0 (+ x1 x2) 1.0e-15 0.0))))))))

(deftest pressure-derivatives-test
  (let [coefficients (prepare-reference-tube)
        n1 0.012
        n2 0.003]
    (testing "Reference derivatives and total pressure gradient along increasing z"
      (doseq [[[p1 p2] [expected-dp1 expected-dp2] expected-g]
              [[[12.0 48.0] [-317.2748197297876 270.7249949427800] -46.54982478700762]
               [[180.0 20.0] [103.4497839020895 -109.3475006174605] -5.897716715371040]
               [[0.0 60.0] [-403.2582534611419 347.2885696089053] -55.96968385223663]]]
        (testing (str "[p_1 p_2] = " [p1 p2] " Pa")
          (let [[dp1 dp2 :as actual] (model/pressure-derivatives coefficients p1 p2 n1 n2)]
            (is (= 2 (count actual)))
            (is (close? expected-dp1 dp1))
            (is (close? expected-dp2 dp2))
            (is (close? expected-g (+ dp1 dp2)))))))
    (testing "Zero flux, superposition and homogeneity at fixed pressures"
      (let [rhs (fn [flux1 flux2]
                  (model/pressure-derivatives coefficients 12.0 48.0 flux1 flux2))
            base (rhs n1 n2)
            summed (mapv + (rhs n1 0.0) (rhs 0.0 n2))
            scaled (rhs (* 2.5 n1) (* 2.5 n2))]
        (is (= [0.0 0.0] (rhs 0.0 0.0)))
        (doseq [[value superposed multiplied] (map vector base summed scaled)]
          (is (close? value superposed))
          (is (close? (* 2.5 value) multiplied)))))))

(deftest pressure-derivative-limits-test
  (let [coefficients (prepare-reference-tube)
        {:keys [rt ak sqrt-m1 sqrt-m2 dk1 dk2 k12]} reference-constants
        [_ {:keys [x1 x2 mu]}] (first local-reference-cases)
        n1 0.012
        n2 0.003
        flux-weight (+ (* sqrt-m1 n1) (* sqrt-m2 n2))
        mass-weight (+ (* x1 sqrt-m1) (* x2 sqrt-m2))
        ac-times-p (/ (* 32.0 mu rt) (* diameter diameter mass-weight))
        continuum-rt-h-per-p (/ (* rt (- (* x1 n2) (* x2 n1))) k12)
        knudsen-rt-h (* rt (- (/ (* x1 n2) dk2) (/ (* x2 n1) dk1)))]
    (testing "At low positive pressure each gradient approaches its own Knudsen law"
      (let [knudsen [(- (/ (* rt n1) dk1)) (- (/ (* rt n2) dk2))]
            gk (reduce + knudsen)]
        (doseq [p [1.0e-4 1.0e-5]]
          (let [ak-over-ac (/ (* ak p) ac-times-p)
                rhs (model/pressure-derivatives coefficients (* x1 p) (* x2 p) n1 n2)]
            ;; The continuum composition term and the correction to G_K are O(p).
            (doseq [[x expected actual] (map vector [x1 x2] knudsen rhs)]
              (let [error-bound (+ (Math/abs (* p continuum-rt-h-per-p))
                                   (* x (Math/abs gk) ak-over-ac))]
                (is (close? expected actual error-bound 1.0e-12) (str "Pressure " p " Pa"))))))))
    (testing "At high finite pressure the gradients approach the continuum equations"
      (doseq [p [1.0e4 1.0e5]]
        (let [ac (/ ac-times-p p)
              gc (- (* ac flux-weight))
              continuum-rt-h (* p continuum-rt-h-per-p)
              expected [(+ (* x1 gc) continuum-rt-h) (- (* x2 gc) continuum-rt-h)]
              [dp1 dp2 :as rhs] (model/pressure-derivatives coefficients (* x1 p) (* x2 p) n1 n2)
              g-error-bound (* (Math/abs gc) (/ ac ak))]
          ;; The Knudsen composition contribution is O(1) against an O(p)
          ;; continuum gradient. The relative correction to G_C is O(1/p).
          (doseq [[x reference actual] (map vector [x1 x2] expected rhs)]
            (let [error-bound (+ (Math/abs knudsen-rt-h) (* x g-error-bound))]
              (is (close? reference actual error-bound 1.0e-12) (str "Pressure " p " Pa"))))
          ;; Summing the large opposite gradients loses precision on the scale
          ;; of the individual derivatives, so allow a small absolute roundoff.
          (let [roundoff (* 16.0 (Math/ulp (double continuum-rt-h)))]
            (is (close? gc (+ dp1 dp2) (+ g-error-bound roundoff) 1.0e-12)
                (str "Total pressure gradient at " p " Pa"))))))))
