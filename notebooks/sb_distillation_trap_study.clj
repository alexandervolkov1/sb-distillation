(ns sb-distillation-trap-study
  (:require [scicloj.kindly.v4.kind :as kind]
            [sb-distillation.plot :as plot]))

;; Preliminary model assumptions and units are encoded below.
;; SI units internally; charts use mm, Celsius, days, and Torr.

(def gas-constant 8.314462618)
(def boltzmann-constant 1.380649e-23)
(def avogadro-constant 6.02214076e23)
(def pa-per-torr 133.32236842105263)
(def h2-molar-mass 0.00201588)
(def sb4-molar-mass 0.48704)
(def h2-collision-diameter-m 0.289e-9)
(def sb4-collision-diameter-m 0.6e-9)

(def mass-kg 0.5)
(def tube-length-m 0.5)
(def h2-flow-ml-h 200.0)
(def gauge-h2-torr 0.03)
(def reference-temp-k 293.15)
(def reference-pressure-pa 101325.0)

(def diameter-grid-mm (mapv #(/ % 2.0) (range 16 33))) ; 8.0 to 16.0 mm
(def temperature-grid-c (mapv double (range 650 751 2)))
(def x-domain [8.0 16.0])
(def y-domain [650.0 750.0])
(def major-levels [0.25 0.5 0.75 1.0 2.0 3.0 5.0 7.0 10.0
                   14.0 20.0 30.0 40.0 60.0 80.0])
(def minor-levels [4.0 6.0 8.0 9.0 12.0 16.0 18.0 25.0 35.0 50.0 70.0])

(defn liquid-sb-vapor-pressure-pa [temp-c]
  (let [temp-k (+ temp-c 273.15)]
    (* 1.0e5 (Math/pow 10.0
                       (- 2.26041 (/ 4475.449 (- temp-k 152.352)))))))

(defn sb4-viscosity [temp-k]
  (let [molecule-mass (/ sb4-molar-mass avogadro-constant)
        d sb4-collision-diameter-m]
    (/ (* 5.0 (Math/sqrt (/ (* molecule-mass boltzmann-constant temp-k)
                             Math/PI)))
       (* 16.0 d d))))

(defn h2-viscosity [temp-k]
  (* 8.76e-6
     (Math/pow (/ temp-k 293.85) 1.5)
     (/ (+ 293.85 72.0) (+ temp-k 72.0))))

(defn wilke-phi [eta-i eta-j molar-mass-i molar-mass-j]
  (/ (Math/pow (+ 1.0
                  (* (Math/sqrt (/ eta-i eta-j))
                     (Math/pow (/ molar-mass-j molar-mass-i) 0.25)))
               2.0)
     (Math/sqrt (* 8.0 (+ 1.0 (/ molar-mass-i molar-mass-j))))))

(defn mixture-viscosity [temp-k h2-fraction]
  (let [sb4-fraction (- 1.0 h2-fraction)
        eta-h (h2-viscosity temp-k)
        eta-s (sb4-viscosity temp-k)
        phi-hs (wilke-phi eta-h eta-s h2-molar-mass sb4-molar-mass)
        phi-sh (wilke-phi eta-s eta-h sb4-molar-mass h2-molar-mass)]
    (+ (/ (* h2-fraction eta-h)
          (+ h2-fraction (* sb4-fraction phi-hs)))
       (/ (* sb4-fraction eta-s)
          (+ sb4-fraction (* h2-fraction phi-sh))))))

(defn knudsen-diffusivity [temp-k diameter-m molar-mass]
  (* (/ diameter-m 3.0)
     (Math/sqrt (/ (* 8.0 gas-constant temp-k)
                   (* Math/PI molar-mass)))))

(defn binary-diffusivity [temp-k pressure-pa]
  (let [m-h (/ h2-molar-mass avogadro-constant)
        m-s (/ sb4-molar-mass avogadro-constant)
        reduced-mass (/ (* m-h m-s) (+ m-h m-s))
        collision-diameter-m (/ (+ h2-collision-diameter-m
                                   sb4-collision-diameter-m) 2.0)]
    (/ (* 3.0 boltzmann-constant temp-k
          (Math/sqrt (/ (* 2.0 boltzmann-constant temp-k)
                        (* Math/PI reduced-mass))))
       (* 16.0 pressure-pa collision-diameter-m collision-diameter-m))))

;; Pure Sb4 calculated with the old transitional conductance formula.
(defn pure-old-mass-flow [temp-c diameter-m]
  (let [temp-k (+ temp-c 273.15)
        p (liquid-sb-vapor-pressure-pa temp-c)
        viscosity (sb4-viscosity temp-k)
        mean-free-path (/ (* boltzmann-constant temp-k)
                          (* (Math/sqrt 2.0) Math/PI
                             sb4-collision-diameter-m sb4-collision-diameter-m
                             (/ p 2.0)))
        ratio (/ diameter-m 2.0 mean-free-path)
        z (/ (+ 1.0 (* 2.507 ratio)) (+ 1.0 (* 3.095 ratio)))
        c-mol (* (/ (* Math/PI (Math/pow diameter-m 3.0))
                    (* 12.0 tube-length-m))
                 (Math/sqrt (/ (* 8.0 gas-constant temp-k)
                               (* Math/PI sb4-molar-mass))))
        c-vis (/ (* Math/PI (Math/pow diameter-m 4.0) (/ p 2.0))
                 (* 128.0 viscosity tube-length-m))]
    (* (/ sb4-molar-mass (* gas-constant temp-k))
       (+ c-vis (* z c-mol)) p)))

;; Pure Sb4 calculated independently from the DGM single-component limit.
(defn pure-dgm-mass-flow [temp-c diameter-m]
  (let [temp-k (+ temp-c 273.15)
        p (liquid-sb-vapor-pressure-pa temp-c)
        area (/ (* Math/PI diameter-m diameter-m) 4.0)
        dk (knudsen-diffusivity temp-k diameter-m sb4-molar-mass)
        b0 (/ (* diameter-m diameter-m) 32.0)
        viscosity (sb4-viscosity temp-k)
        flux (/ (+ (* dk p) (/ (* b0 p p) (* 2.0 viscosity)))
                (* gas-constant temp-k tube-length-m))]
    (* flux area sb4-molar-mass)))

(defn hydrogen-molar-flow [flow-ml-h]
  (/ (* reference-pressure-pa flow-ml-h 1.0e-6)
     (* 3600.0 gas-constant reference-temp-k)))

;; C_trap [L/s] is effective at the gauge-side/reference temperature.
;; nil represents zero pressure drop across an ideal trap.
;; Only H2 is assumed to pass through the trap.
(defn receiver-h2-torr [trap-conductance-l-s]
  (let [q-torr-l-s (/ (* (/ reference-pressure-pa pa-per-torr)
                         h2-flow-ml-h)
                      3.6e6)]
    (+ gauge-h2-torr
       (if (nil? trap-conductance-l-s)
         0.0
         (/ q-torr-l-s trap-conductance-l-s)))))

;; Binary DGM pressure gradients, integrated from receiver to source.
;; Pressure-times-binary-diffusivity is constant in the hard-sphere model.
(defn dgm-properties [temp-k diameter-m]
  (let [eta-h (h2-viscosity temp-k)
        eta-s (sb4-viscosity temp-k)]
    {:eta-h eta-h
     :eta-s eta-s
     :phi-hs (wilke-phi eta-h eta-s h2-molar-mass sb4-molar-mass)
     :phi-sh (wilke-phi eta-s eta-h sb4-molar-mass h2-molar-mass)
     :d-h (knudsen-diffusivity temp-k diameter-m h2-molar-mass)
     :d-s (knudsen-diffusivity temp-k diameter-m sb4-molar-mass)
     :d-hs-p (binary-diffusivity temp-k 1.0)
     :b0 (/ (* diameter-m diameter-m) 32.0)
     :rt (* gas-constant temp-k)}))

(defn dgm-gradient [properties h2-flux sb4-flux [p-h p-s]]
  (when-not (and (pos? p-h) (<= 0.0 p-s))
    (throw (ex-info "Nonphysical partial pressure." {:p-h p-h :p-s p-s})))
  (let [{:keys [eta-h eta-s phi-hs phi-sh
                d-h d-s d-hs-p b0 rt]} properties
        p (+ p-h p-s)
        y-h (/ p-h p)
        y-s (/ p-s p)
        eta (+ (/ (* y-h eta-h) (+ y-h (* y-s phi-hs)))
               (/ (* y-s eta-s) (+ y-s (* y-h phi-sh))))
        interaction (/ (- (* p-s h2-flux) (* p-h sb4-flux)) d-hs-p)
        f-h (* rt (+ (/ h2-flux d-h) interaction))
        f-s (* rt (- (/ sb4-flux d-s) interaction))
        a-h (/ (* b0 p-h) (* eta d-h))
        a-s (/ (* b0 p-s) (* eta d-s))
        total-gradient (/ (+ f-h f-s) (+ 1.0 a-h a-s))]
    [(- f-h (* a-h total-gradient))
     (- f-s (* a-s total-gradient))]))

(defn v-add [a b factor]
  [(+ (nth a 0) (* factor (nth b 0)))
   (+ (nth a 1) (* factor (nth b 1)))])

(defn rk4-step [gradient state ds]
  (let [k1 (gradient state)
        k2 (gradient (v-add state k1 (* ds 0.5)))
        k3 (gradient (v-add state k2 (* ds 0.5)))
        k4 (gradient (v-add state k3 ds))]
    [(+ (nth state 0)
        (* (/ ds 6.0) (+ (nth k1 0) (* 2.0 (nth k2 0))
                         (* 2.0 (nth k3 0)) (nth k4 0))))
     (+ (nth state 1)
        (* (/ ds 6.0) (+ (nth k1 1) (* 2.0 (nth k2 1))
                         (* 2.0 (nth k3 1)) (nth k4 1))))]))

(defn source-partial-pressures
  [temp-c diameter-m receiver-torr sb4-molar-flow]
  (let [temp-k (+ temp-c 273.15)
        area (/ (* Math/PI diameter-m diameter-m) 4.0)
        nh (/ (hydrogen-molar-flow h2-flow-ml-h) area)
        ns (/ sb4-molar-flow area)
        steps 120
        ds (/ tube-length-m steps)
        properties (dgm-properties temp-k diameter-m)
        gradient #(dgm-gradient properties nh ns %)]
    (loop [i 0
           state [(* receiver-torr pa-per-torr) 0.0]]
      (if (= i steps)
        state
        (recur (inc i) (rk4-step gradient state ds))))))

;; Find Sb4 flow by matching its source partial pressure to saturation.
;; Unphysical trial states can be used as an upper bracket only;
;; the final solution must be physical and satisfy the boundary condition.
(defn mixed-dgm-mass-flow [temp-c diameter-m receiver-torr]
  (let [target (liquid-sb-vapor-pressure-pa temp-c)
        initial (/ (pure-dgm-mass-flow temp-c diameter-m)
                   sb4-molar-mass)
        error
        (fn [sb4-molar-flow]
          (try
            (- (nth (source-partial-pressures
                     temp-c diameter-m receiver-torr sb4-molar-flow) 1)
               target)
            (catch clojure.lang.ExceptionInfo _ ##Inf)))
        high
        (loop [h initial
               i 0]
          (if (or (>= (error h) 0.0) (= i 20))
            h
            (recur (* h 1.5) (inc i))))]
    (when-not (>= (error high) 0.0)
      (throw (ex-info "Could not bracket Sb4 flow."
                      {:temp-c temp-c :diameter-m diameter-m
                       :receiver-torr receiver-torr})))
    (let [flow
          (loop [lo 0.0 hi high i 0]
            (if (= i 20)
              (* 0.5 (+ lo hi))
              (let [mid (* 0.5 (+ lo hi))]
                (if (neg? (error mid))
                  (recur mid hi (inc i))
                  (recur lo mid (inc i))))))
          p-source (source-partial-pressures
                    temp-c diameter-m receiver-torr flow)]
      (when (> (/ (Math/abs (- (nth p-source 1) target)) target) 1.0e-3)
        (throw (ex-info "Boundary condition not satisfied."
                        {:temp-c temp-c :diameter-m diameter-m
                         :receiver-torr receiver-torr
                         :p-source p-source :target target})))
      (* flow sb4-molar-mass))))

(def scenarios
  [{:id :vacuum-old :title "Vacuum: original Knudsen model"
    :mode :pure-old}
   {:id :vacuum-dgm :title "Vacuum: DGM limit"
    :mode :pure-dgm}
   {:id :h2-free :title "H2: ideal trap (0.03 Torr)"
    :mode :mixture :conductance nil}
   {:id :h2-light :title "H2: trap conductance 0.60 L/s"
    :mode :mixture :conductance 0.60}
   {:id :h2-medium :title "H2: trap conductance 0.16 L/s"
    :mode :mixture :conductance 0.16}
   {:id :h2-strong :title "H2: trap conductance 0.044 L/s"
    :mode :mixture :conductance 0.044}
   {:id :h2-severe :title "H2: trap conductance 0.014 L/s"
    :mode :mixture :conductance 0.014}])

(defn scenario-flow [scenario temp-c diameter-mm]
  (let [diameter-m (* diameter-mm 1.0e-3)]
    (case (:mode scenario)
      :pure-old (pure-old-mass-flow temp-c diameter-m)
      :pure-dgm (pure-dgm-mass-flow temp-c diameter-m)
      :mixture (mixed-dgm-mass-flow
                temp-c diameter-m
                (receiver-h2-torr (:conductance scenario))))))

(defn compute-surface [scenario]
  (let [rows
        (mapv
         (fn [temp-c]
           (mapv
            (fn [diameter-mm]
              (let [mass-flow (scenario-flow scenario temp-c diameter-mm)]
                {:temp-c temp-c
                 :diameter-mm diameter-mm
                 :flow-g-h (* 3.6e6 mass-flow)
                 :days (/ mass-kg (* 86400.0 mass-flow))}))
            diameter-grid-mm))
         temperature-grid-c)]
    (assoc scenario :rows rows)))

;; Surfaces are computed once and then reused by plots and tables.
(def surfaces (mapv compute-surface scenarios))

(defn level-label [days]
  (if (< days 1.0)
    (str (int (* 24.0 days)) " h")
    (str (int days) " d")))

;; Interpolate an isocontour between adjacent diameter nodes.
;; This avoids rerunning the DGM solver for every contour.
(defn contour-points [rows levels]
  (vec
   (for [level levels
         row rows
         :let [pair
               (first
                (filter
                 (fn [[a b]]
                   (<= (min (:days a) (:days b))
                       level
                       (max (:days a) (:days b))))
                 (partition 2 1 row)))]
         :when pair
         :let [[a b] pair
               log-a (Math/log (:days a))
               log-b (Math/log (:days b))
               ratio (if (= log-a log-b)
                       0.5
                       (/ (- (Math/log level) log-a)
                          (- log-b log-a)))]]
     {:temp-c (:temp-c a)
      :diameter-mm (+ (:diameter-mm a)
                      (* ratio (- (:diameter-mm b) (:diameter-mm a))))
      :days level
      :label (level-label level)})))

(defn contour-labels [points levels]
  (vec
   (keep
    (fn [level]
      (let [curve (filterv #(= level (:days %)) points)]
        (when (seq curve)
          (apply min-key
                 (fn [p]
                   (let [x (/ (- (:diameter-mm p) 8.0) 8.0)
                         y (/ (- (:temp-c p) 650.0) 100.0)]
                     (Math/abs (- x y))))
                 curve))))
    levels)))

(def fill-scale
  {:type "log" :domain [0.1 100.0]
   :range ["#1a9850" "#91cf60" "#ffffbf" "#fc8d59" "#d73027"]
   :clamp true})

(def line-scale
  {:type "log" :domain [0.1 100.0]
   :range ["#0b5d2a" "#4f7f2a" "#b39b00" "#a64b1f" "#7f1d1d"]
   :clamp true})

(defn plot-surface [surface]
  (let [rows (:rows surface)
        values (vec
                (for [row rows
                      p row]
                  (let [d (:diameter-mm p)
                        t (:temp-c p)]
                    (assoc p
                           :d0 (max 8.0 (- d 0.25))
                           :d1 (min 16.0 (+ d 0.26))
                           :t0 (max 650.0 (- t 1.0))
                           :t1 (min 750.0 (+ t 1.04))))))
        major (contour-points rows major-levels)
        minor (contour-points rows minor-levels)
        labels (contour-labels major major-levels)
        fill-layer
        {:data {:values values}
         :mark {:type "rect" :stroke nil :strokeWidth 0}
         :encoding
         {:x {:field "d0" :type "quantitative"
              :scale {:domain x-domain}
              :axis {:title "Tube diameter (mm)" :values [8 10 12 14 16]}}
          :x2 {:field "d1"}
          :y {:field "t0" :type "quantitative"
              :scale {:domain y-domain}
              :axis {:title "Source temperature (°C)"
                     :values [650 670 690 710 730 750]}}
          :y2 {:field "t1"}
          :color {:field "days" :type "quantitative"
                  :scale fill-scale
                  :legend {:title "Days" :orient "top-right"
                           :direction "vertical" :offset 12
                           :gradientLength 180 :gradientThickness 13
                           :values [0.1 0.25 1.0 5.0 20.0 100.0]
                           :labelFontSize 13 :titleFontSize 14
                           :fillColor "white" :strokeColor "#cccccc"
                           :padding 7}}
          :tooltip [{:field "diameter-mm" :type "quantitative"
                     :title "Diameter (mm)" :format ".1f"}
                    {:field "temp-c" :type "quantitative"
                     :title "Temperature (°C)" :format ".0f"}
                    {:field "days" :type "quantitative"
                     :title "Duration (days)" :format ".2f"}]}}
        line-layer
        (fn [points width opacity]
          {:data {:values points}
           :mark {:type "line" :strokeWidth width :opacity opacity}
           :encoding {:x {:field "diameter-mm" :type "quantitative"
                          :scale {:domain x-domain}}
                      :y {:field "temp-c" :type "quantitative"
                          :scale {:domain y-domain}}
                      :detail {:field "days" :type "nominal"}
                      :order {:field "temp-c" :type "quantitative"
                              :sort "ascending"}
                      :color {:field "days" :type "quantitative"
                              :scale line-scale :legend nil}}})
        text-layer
        {:data {:values labels}
         :mark {:type "text" :color "black" :fontSize 15
                :fontWeight 900 :align "center" :baseline "middle"
                :stroke "white" :strokeWidth 0.4}
         :encoding {:x {:field "diameter-mm" :type "quantitative"
                        :scale {:domain x-domain}}
                    :y {:field "temp-c" :type "quantitative"
                        :scale {:domain y-domain}}
                    :text {:field "label" :type "nominal"}}}]
    (kind/vega-lite
     (merge plot/theme
            {:width 1040 :height 710
             :autosize {:type "fit" :contains "padding"}
             :title (:title surface)
             :resolve {:scale {:color "independent"}}
             :layer [fill-layer
                     (plot/grid-layer :x [8 10 12 14 16]
                                      {:domain x-domain} true)
                     (plot/grid-layer :y [650 670 690 710 730 750]
                                      {:domain y-domain} true)
                     (line-layer minor 1.2 0.60)
                     (line-layer major 2.4 0.95)
                     text-layer]}))))

;; Numerical comparison at a representative temperature and diameter.
(kind/table
 (mapv
  (fn [surface]
    (let [row (first (filter #(= 690.0 (:temp-c (first %)))
                             (:rows surface)))
          point (first (filter #(= 8.0 (:diameter-mm %)) row))]
      {:scenario (:title surface)
       :receiver-h2-torr (when (= (:mode surface) :mixture)
                           (receiver-h2-torr (:conductance surface)))
       :sb4-g-h (:flow-g-h point)
       :days (:days point)}))
  surfaces))

;; Each chart uses the same diameter, temperature, and logarithmic color scales.
(plot-surface (nth surfaces 0))
(plot-surface (nth surfaces 1))
(plot-surface (nth surfaces 2))
(plot-surface (nth surfaces 3))
(plot-surface (nth surfaces 4))
(plot-surface (nth surfaces 5))
(plot-surface (nth surfaces 6))
