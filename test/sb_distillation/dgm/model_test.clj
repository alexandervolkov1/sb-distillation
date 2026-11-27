(ns sb-distillation.dgm.model-test
  (:require [clojure.test :refer [deftest is testing]]
            [sb-distillation.cpim.model :as cpim]
            [sb-distillation.dgm.model :as dgm]
            [sb-distillation.physics.common :as physics]
            [sb-distillation.physics.transport :as transport]
            [sb-distillation.solver :as solver]))

(def ^:private temperature 973.15)
(def ^:private diameter 0.012)

;; SI coefficients for the default tube, evaluated with 60-digit Python Decimal
;; arithmetic from the transport equations in theory revision e6a744f3a66f.
;; Rational arithmetic keeps the downstream reference independent of the
;; double-precision pressure resistance and gradient functions.
(def ^:private reference-transport
  (update-vals {:rt 8091.21929670670M
                :dk1 0.82272504694499487197M
                :dk2 12.788060917626169289M
                :k12 20.423607122821020799M
                :mu1 0.000051051791113053772884M
                :mu2 0.000018480328147503648549M
                :phi12 0.045872049077639552702M
                :phi21 4.0118694277337886753M
                :b0 0.0000045M}
               rationalize))

(defn- coefficients []
  (dgm/prepare-coefficients
   (physics/prepare-coefficients transport/default-parameters temperature diameter)))

(defn- reference-state [[p1 p2]]
  (let [{:keys [mu1 mu2 phi12 phi21 b0 dk1 dk2]} reference-transport
        pressure (+ (rationalize p1) (rationalize p2))
        x1 (/ (rationalize p1) pressure)
        x2 (/ (rationalize p2) pressure)
        mu (+ (/ (* x1 mu1) (+ x1 (* x2 phi12)))
              (/ (* x2 mu2) (+ x2 (* x1 phi21))))
        b-over-mu (/ b0 mu)]
    {:pressure pressure
     :fractions [x1 x2]
     :b-over-mu b-over-mu
     :beta-per-pa [(/ (* x1 b-over-mu) dk1) (/ (* x2 b-over-mu) dk2)]}))

(defn- reference-gradient
  "G - DGM total pressure gradient in Pa/m from the summed species equations.
  Pressures in Pa and molar fluxes in mol/(m^2*s)."
  [{:keys [pressure beta-per-pa]} [n1 n2]]
  (let [{:keys [rt dk1 dk2]} reference-transport
        knudsen (+ (/ (rationalize n1) dk1) (/ (rationalize n2) dk2))]
    (/ (* (- rt) knudsen) (+ 1 (* pressure (reduce + beta-per-pa))))))

(defn- close? [expected actual tolerance]
  (<= (abs (- (double expected) (double actual))) tolerance))

(defn- roundoff [values]
  (* 32.0 (apply max (map #(Math/ulp (double %)) values))))

(deftest derivative-references-test
  (let [tube (coefficients)]
    ;; The partial gradients are 60-digit Decimal evaluations of the DGM
    ;; species equations, rounded here. The cases span mixed and pure states,
    ;; both signs of total pressure gradient, and signed component fluxes.
    (doseq [[[p1 p2 :as state] [n1 n2 :as fluxes] [expected1 expected2]]
            [[[12.0 48.0] [0.012 0.003] [-273.598987059751087 227.049162272743468]]
             [[180.0 20.0] [0.012 0.003] [114.042169580117514 -119.939886295488554]]
             [[8.0 12.0] [-0.004 0.010] [75.7145096246314903 -58.4201829278261126]]
             [[0.0 60.0] [0.012 0.003] [-403.258253461141908 347.288569608905278]]
             [[60.0 0.0] [0.012 0.003] [57.0659369725898399 -73.2087399826371669]]]]
      (testing (pr-str [state fluxes])
        (let [[dp1 dp2 :as dgm-rhs] (dgm/pressure-derivatives tube p1 p2 n1 n2)
              cpim-rhs (cpim/pressure-derivatives tube p1 p2 n1 n2)
              {:keys [pressure fractions b-over-mu] :as local} (reference-state state)
              {:keys [dk1 dk2]} reference-transport
              g (reference-gradient local fluxes)
              [x1 x2] fractions
              delta (* x1 x2 pressure b-over-mu (- (/ dk1) (/ dk2)) g)
              [delta1 delta2] (map - cpim-rhs dgm-rhs)]
          (is (= 2 (count dgm-rhs)))
          (is (close? expected1 dp1 1.0e-10))
          (is (close? expected2 dp2 1.0e-10))
          (is (close? g (+ dp1 dp2) 1.0e-10))
          (is (close? g (reduce + cpim-rhs) 1.0e-10))
          (is (close? delta delta1 1.0e-10))
          (is (close? (- delta) delta2 1.0e-10)))))))

(deftest permeability-test
  (let [shared (physics/prepare-coefficients transport/default-parameters temperature diameter)
        prepared (dgm/prepare-coefficients shared)]
    (is (close? (:b0 reference-transport) (:b0 prepared) 1.0e-20))
    (is (close? (* 4.0 (:b0 reference-transport))
                (:b0 (dgm/prepare-coefficients (assoc shared :diameter 0.024))) 1.0e-20))
    (is (= shared (dissoc prepared :b0))))
  (is (close? 7.03125e306 (:b0 (dgm/prepare-coefficients {:diameter 1.5e154})) 1.0e292))
  (doseq [d [nil "0.012" -0.012 -1e-500M 0.0 1e-500M ##Inf ##-Inf ##NaN]]
    (let [{:keys [parameter value]}
          (try (dgm/prepare-coefficients {:diameter d})
               (catch clojure.lang.ExceptionInfo ex (ex-data ex)))]
      (is (= :diameter parameter) (str "Diameter " d " m"))
      (is (= (pr-str d) (pr-str value)))))
  (doseq [[d invalid?] [[1.0e-200 zero?] [1.0e200 #(Double/isInfinite %)]]]
    (let [{:keys [quantity value]}
          (try (dgm/prepare-coefficients {:diameter d})
               (catch clojure.lang.ExceptionInfo ex (ex-data ex)))]
      (is (= :b0 quantity) (str "Diameter " d " m"))
      (is (and (number? value) (invalid? value)) (str "Diameter " d " m")))))

(deftest flux-linearity-test
  (let [tube (coefficients)
        rhs #(dgm/pressure-derivatives tube 12.0 48.0 %1 %2)
        base (rhs 0.012 0.003)]
    (is (every? zero? (rhs 0.0 0.0)))
    (doseq [[value superposed scaled]
            (map vector base (map + (rhs 0.012 0.0) (rhs 0.0 0.003)) (rhs -0.03 -0.0075))]
      (is (close? value superposed 1.0e-10))
      (is (close? (* -2.5 value) scaled 1.0e-10)))))

(deftest diameter-and-pressure-scaling-test
  (let [rhs (dgm/pressure-derivatives (coefficients) 12.0 48.0 0.012 0.003)]
    (doseq [factor [0.5 2.0]
            :let [tube (dgm/prepare-coefficients
                        (physics/prepare-coefficients transport/default-parameters
                                                      temperature (* diameter factor)))
                  scaled (dgm/pressure-derivatives tube (/ 12.0 factor) (/ 48.0 factor)
                                                   0.012 0.003)]
            [original actual] (map vector rhs scaled)]
      ;; Scaling d and 1/p together scales both diffusivities and preserves beta_i.
      (is (close? (/ original factor) actual 1.0e-10)))))

(deftest pure-component-and-pressure-domain-test
  (let [tube (coefficients)
        {:keys [rt dk1 dk2 mu1 mu2 b0]} reference-transport]
    (doseq [[state fluxes active dk viscosity] [[[60.0 0.0] [0.012 0.0] 0 dk1 mu1]
                                               [[0.0 60.0] [0.0 0.012] 1 dk2 mu2]]]
      (let [gradient (apply dgm/pressure-derivatives tube (concat state fluxes))
            pure-limit (/ (* (- rt) 0.012) (+ dk (/ (* 60 b0) viscosity)))]
        (is (close? pure-limit (gradient active) 1.0e-10))
        (is (zero? (gradient (- 1 active))))))
    (is (solver/valid-pressure-state? [0.0 60.0]))
    (is (solver/valid-pressure-state? [20.0 40.0]))
    (doseq [state [[60.0 0.0] [0.0 0.0] [-1.0 60.0] [60.0 -1.0] [0.0 ##NaN] [0.0 ##Inf]
                   [Double/MAX_VALUE Double/MAX_VALUE]]]
      (is (not (solver/valid-pressure-state? state)) (pr-str state))
      (is (thrown? clojure.lang.ExceptionInfo
                   (solver/integrate-trial dgm/pressure-derivatives
                                           {:length 0.5 :coefficients tube :n2 0.003
                                            :outlet-state state}
                                           0.012))))))

(deftest equal-masses-and-zero-total-gradient-test
  (testing "Equal molecular masses give equal Knudsen diffusivities and both species gradients"
    (let [params (assoc transport/default-parameters
                        :sb4-molar-mass (:h2-molar-mass transport/default-parameters))
          tube (dgm/prepare-coefficients (physics/prepare-coefficients params temperature diameter))]
      (is (= (:dk1 tube) (:dk2 tube)))
      (doseq [[p1 p2] [[12.0 48.0] [180.0 20.0] [0.0 60.0]]
              [expected actual] (map vector
                                     (cpim/pressure-derivatives tube p1 p2 0.012 0.003)
                                     (dgm/pressure-derivatives tube p1 p2 0.012 0.003))]
        (is (close? expected actual 1.0e-10)))))
  (testing "Knudsen-weighted counterflow cancels total pressure without canceling species transport"
    (let [{:keys [dk1 dk2]} reference-transport
          tube (coefficients)
          n1 0.012
          n2 (double (* -0.012 (/ dk2 dk1)))
          [dp1 dp2] (dgm/pressure-derivatives tube 12.0 48.0 n1 n2)
          [cpim1 cpim2] (cpim/pressure-derivatives tube 12.0 48.0 n1 n2)]
      (is (close? 0 (reference-gradient (reference-state [12.0 48.0]) [n1 n2]) 1.0e-10))
      (is (close? 0 (+ dp1 dp2) 1.0e-10))
      (is (< dp1 0.0 dp2))
      (is (close? cpim1 dp1 1.0e-10))
      (is (close? cpim2 dp2 1.0e-10)))))

(deftest molecular-limit-test
  (let [{:keys [rt dk1 dk2 k12]} reference-transport
        tube (coefficients)
        n1 0.012
        n2 0.003
        knudsen [(/ (* (- rt) (rationalize n1)) dk1)
                 (/ (* (- rt) (rationalize n2)) dk2)]
        gk (reduce + knudsen)
        leading (/ (* rt (abs (- (* 3/4 (rationalize n1)) (* 1/4 (rationalize n2))))) k12)
        {:keys [beta-per-pa]} (reference-state [1.0 3.0])]
    (doseq [pressure [1.0e-4 1.0e-5 1.0e-6]
            :let [rhs (dgm/pressure-derivatives tube (* 0.25 pressure) (* 0.75 pressure) n1 n2)]
            [expected actual beta-slope] (map vector knudsen rhs beta-per-pa)]
      ;; The mutual and viscous corrections are O(p), with |G| <= |G_K|.
      (is (close? expected actual
                  (+ (* pressure (+ leading (* beta-slope (abs gk)))) (roundoff rhs)))
          (str "Pressure " pressure " Pa")))))

(deftest continuum-asymptotics-test
  (let [{:keys [rt dk1 dk2 k12]} reference-transport
        tube (coefficients)
        {:keys [b-over-mu beta-per-pa fractions]} (reference-state [1.0 3.0])
        [x1 x2] fractions
        alpha (reduce + beta-per-pa)]
    (doseq [[label [n1 n2] zero-leading?]
            [["Nonzero leading diffusion" [0.03125 0.015625] false]
             ["Fluxes proportional to mole fractions" [0.03125 0.09375] true]]]
      (testing label
        (let [k1 (/ (rationalize n1) dk1)
              k2 (/ (rationalize n2) dk2)
              g-times-p-limit (/ (* (- rt) (+ k1 k2)) alpha)
              mutual1 (/ (* rt x2 (rationalize n1)) k12)
              mutual2 (/ (* rt x1 (rationalize n2)) k12)
              leading1 (- mutual2 mutual1)
              cpim-offset1 (* (- rt) (- (* x2 k1) (* x1 k2)))
              dgm-offsets (mapv #(- (* (- rt) %1) (* %2 g-times-p-limit)) [k1 k2] beta-per-pa)
              delta-limit (* x1 x2 b-over-mu (- (/ dk1) (/ dk2)) g-times-p-limit)]
          (is (= zero-leading? (zero? leading1)))
          (is (> (abs delta-limit) 1))
          (when zero-leading?
            (is (every? zero? dgm-offsets))
            (is (> (abs cpim-offset1) 1)))
          (doseq [pressure [1.0e4 1.0e5 1.0e6]
                  :let [p1 (* 0.25 pressure)
                        p2 (* 0.75 pressure)
                        dgm-rhs (dgm/pressure-derivatives tube p1 p2 n1 n2)
                        cpim-rhs (cpim/pressure-derivatives tube p1 p2 n1 n2)
                        ;; CPIM subtracts the two diffusion terms before forming its
                        ;; derivatives, even when their leading parts cancel.
                        tolerance (roundoff (concat dgm-rhs cpim-rhs
                                                     [(* pressure mutual1) (* pressure mutual2)]))
                        correction (/ (+ 1 (* pressure alpha)))
                        p-g-error (* (abs g-times-p-limit) correction)]]
            (doseq [rhs [dgm-rhs cpim-rhs]]
              ;; DGM can cancel the finite Knudsen and viscous terms to produce
              ;; gradients much smaller than either term.
              (is (close? g-times-p-limit (* pressure (reduce + rhs))
                          (+ p-g-error (* pressure (roundoff (concat rhs [(* rt k1) (* rt k2)])))))))
            (doseq [[slope dgm-value cpim-value dgm-offset cpim-offset beta-slope x]
                    (map vector [leading1 (- leading1)] dgm-rhs cpim-rhs
                         dgm-offsets [cpim-offset1 (- cpim-offset1)] beta-per-pa fractions)]
              (is (close? dgm-offset (- dgm-value (* pressure slope))
                          (+ (* beta-slope p-g-error) tolerance)))
              (is (close? cpim-offset (- cpim-value (* pressure slope))
                          (+ (/ (* x (abs g-times-p-limit)) pressure) tolerance))))
            (when zero-leading?
              ;; N_i = x_i N gives dp_i/dz = -RT (N_i/D_Ki)/(1 + p*alpha).
              (doseq [[knudsen actual] (map vector [k1 k2] dgm-rhs)]
                (let [p-gradient-limit (/ (* (- rt) knudsen) alpha)]
                  (is (close? p-gradient-limit (* pressure actual)
                              (+ (* (abs p-gradient-limit) correction)
                                 (* pressure (roundoff [(* rt knudsen)]))))))))
            (let [delta1 (- (first cpim-rhs) (first dgm-rhs))]
              (is (close? delta-limit delta1 (+ (* (abs delta-limit) correction) tolerance)))
              (is (> (abs delta1) (* 0.9 (abs delta-limit)))))))))))
