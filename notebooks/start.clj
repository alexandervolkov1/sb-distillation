(ns start
  (:require [scicloj.kindly.v4.kind :as kind]))

;; # Модель вакуумной дистилляции сурьмы

(def setup
  {:amp-diameter 0.06
   :amp-length 0.180
   :tube-diameter 0.008
   :tube-length 0.500
   :tube-overlap 0.100
   :heater-pitch 0.010
   :sb-mass 0.500
   :sb-density 6400.0
   :source-temp 700.0
   :receiver-temp 500.0
   :tube-temp 700.0})

(defn ampoule-volume
  [{:keys [amp-diameter amp-length]}]
  (let [r (/ amp-diameter 2.0)
        l amp-length]
    (+ (* Math/PI r r l)
       (* 4/3 Math/PI r r r))))

(* 1000 (ampoule-volume setup))

(defn sb-volume
  [{:keys [sb-mass sb-density]}]
  (/ sb-mass sb-density))

{:sb-volume-ml (* 1.0e6 (sb-volume setup))
 :filling-percent (* 100.0
                     (/ (sb-volume setup)
                        (ampoule-volume setup)))}

(defn filled-volume
  [{:keys [amp-diameter amp-length]} h]
  (let [r (/ amp-diameter 2.0)
        a (Math/sqrt (* h (- (* 2.0 r) h)))

        segment (- (* r r (Math/acos (/ (- r h) r)))
                   (* (- r h) a))

        cylinder (* amp-length segment)
        caps (* Math/PI h h (- r (/ h 3.0)))]

    (+ cylinder caps)))

(* 1000 (filled-volume setup 0.03))

(defn melt-height [params]
  (let [target (sb-volume params)]
    (loop [low 0.0
           high (:amp-diameter params)
           n 0]
      (let [mid (/ (+ low high) 2.0)]
        (if (>= n 60)
          mid
          (if (< (filled-volume params mid) target)
            (recur mid high (inc n))
            (recur low mid (inc n))))))))

(* 1000 (melt-height setup))

(defn evaporation-area
  [{:keys [amp-diameter amp-length]} h]
  (let [r (/ amp-diameter 2.0)
        a (Math/sqrt (* h (- (* 2.0 r) h)))]
    (+ (* 2.0 amp-length a)
       (* Math/PI a a))))

(* 1.0e4 (evaporation-area setup (melt-height setup)))

(defn sb-vapor-pressure
  [temp-c]
  (let [temp-k (+ temp-c 273.15)
        pressure-bar
        (Math/pow 10.0
                  (- 2.26041
                     (/ 4475.449 (- temp-k 152.352))))]
    (* pressure-bar 750.061683)))

(sb-vapor-pressure (:source-temp setup))

(let [points
      (mapv (fn [t]
              {:temp-c t
               :pressure-torr (sb-vapor-pressure t)})
            (range 650 751 5))]

  (kind/vega-lite
   {:data {:values points}
    :mark {:type "line" :point true}
    :width 650
    :height 350
    :encoding
    {:x {:field "temp-c"
         :type "quantitative"
         :title "Температура, °C"}
     :y {:field "pressure-torr"
         :type "quantitative"
         :title "Давление, Torr"}}}))

(defn sb-solid-vapor-pressure
  [temp-c]
  (let [temp-k (+ temp-c 273.15)]
    (Math/pow 10.0
              (- 10.571 (/ 10300.0 temp-k)))))

{:source (sb-vapor-pressure (:source-temp setup))
 :receiver (sb-solid-vapor-pressure (:receiver-temp setup))}

(defn mean-free-path
  [temp-c pressure-torr molecule-diameter]
  (let [k 1.380649e-23
        temp-k (+ temp-c 273.15)
        pressure-pa (* pressure-torr 133.322368)]
    (/ (* k temp-k)
       (* (Math/sqrt 2.0)
          Math/PI
          molecule-diameter
          molecule-diameter
          pressure-pa))))

(let [lambda (mean-free-path
              (:source-temp setup)
              (sb-vapor-pressure (:source-temp setup))
              0.6e-9)]
  {:mean-free-path-mm (* 1000 lambda)
   :knudsen (/ lambda (:tube-diameter setup))})

(defn molecular-conductance
  [{:keys [tube-diameter tube-length tube-temp]}]
  (let [gas-constant 8.314462618
        sb4-molar-mass (* 4.0 0.12176)
        temp-k (+ tube-temp 273.15)

        mean-speed
        (Math/sqrt
         (/ (* 8.0 gas-constant temp-k)
            (* Math/PI sb4-molar-mass)))]

    (/ (* Math/PI
          (Math/pow tube-diameter 3.0)
          mean-speed)
       (* 12.0 tube-length))))

(* 1000 (molecular-conductance setup))

(defn sb4-viscosity
  [temp-c molecule-diameter]
  (let [k 1.380649e-23
        avogadro 6.02214076e23
        molar-mass (* 4.0 0.12176)
        molecule-mass (/ molar-mass avogadro)
        temp-k (+ temp-c 273.15)]

    (* (/ 5.0
          (* 16.0 molecule-diameter molecule-diameter))
       (Math/sqrt
        (/ (* molecule-mass k temp-k)
           Math/PI)))))

(sb4-viscosity (:tube-temp setup) 0.6e-9)

(defn viscous-conductance
  [{:keys [tube-diameter tube-length tube-temp
           source-temp receiver-temp]}]

  (let [p1 (* 133.322368
              (sb-vapor-pressure source-temp))

        p2 (* 133.322368
              (sb-solid-vapor-pressure receiver-temp))

        p-mean (/ (+ p1 p2) 2.0)

        eta (sb4-viscosity tube-temp 0.6e-9)]

    (/ (* Math/PI
          (Math/pow tube-diameter 4.0)
          p-mean)
       (* 128.0 eta tube-length))))

{:molecular-l-s (* 1000 (molecular-conductance setup))
 :viscous-l-s   (* 1000 (viscous-conductance setup))}

(defn transition-conductance
  [params]
  (let [cm (molecular-conductance params)
        cv (viscous-conductance params)

        ratio (/ cv cm)
        x (* ratio (/ 12.1 135.0))

        correction
        (/ (+ 1.0 (* 192.0 x))
           (+ 1.0 (* 237.0 x)))]

    (+ cv (* cm correction))))

(* 1000 (transition-conductance setup))

(defn sb-transfer-rate
  [{:keys [source-temp receiver-temp tube-temp] :as params}]
  (let [p1 (sb-vapor-pressure source-temp)
        p2 (sb-solid-vapor-pressure receiver-temp)

        pressure-drop-pa (* 133.322368 (- p1 p2))

        throughput (* (transition-conductance params)
                      pressure-drop-pa)

        temp-k (+ tube-temp 273.15)
        molar-mass (* 4.0 0.12176)

        mol-per-sec (/ throughput (* 8.314462618 temp-k))
        kg-per-sec (* mol-per-sec molar-mass)]

    (* kg-per-sec 3.6e6)))

(let [rate (sb-transfer-rate setup)
      hours (/ (* 1000.0 (:sb-mass setup)) rate)]

  {:rate-g-per-hour rate
   :time-hours hours
   :time-days (/ hours 24.0)})
