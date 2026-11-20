(ns sb-distillation
  (:require [scicloj.kindly.v4.kind :as kind]
            [sb-distillation.plot :as plot]))

;; # Vacuum Distillation of Antimony
;;
;; A preliminary model of vacuum distillation
;; of pure antimony without a carrier gas.
;;
;; ## 1. Ampoule Volume
;;
;; The ampoule consists of a cylindrical section and two
;; hemispherical end caps. The internal volume of the ampoule
;; [$m^3$] is:
;;
;; $$V = \pi r^2 L + \frac{4}{3}\pi r^3$$
;;
;; where:
;; - $r$ is the internal radius [$m$]
;; - $L$ is the straight cylindrical length [$m$]

(defn amp-volume
  [radius straight-length]
  (+ (* Math/PI radius radius straight-length)
     (* 4/3 Math/PI radius radius radius)))

;; ## 2. Liquid Fill Volume
;;
;; The ampoule is horizontal, and the liquid surface is assumed
;; to be flat. The fill height $h$ [$m$] is measured from the bottom.
;; The liquid volume consists of two parts:
;; - the volume within the cylindrical section;
;; - the combined volume within the two hemispherical end caps.
;;
;; The cross-sectional area of the liquid-filled segment in the cylindrical
;; section of the ampoule [$m^2$] is:
;;
;; $$A(h) = r^2\arccos\left(\frac{r-h}{r}\right)
;;    - (r-h)\sqrt{2rh-h^2}$$
;;
;; The total liquid volume [$m^3$] is:
;;
;; $$V(h) = L A(h) + \pi h^2\left(r-\frac{h}{3}\right)$$
;;
;; Valid for $0 \leqslant h \leqslant 2r$.

(defn filled-volume
  [radius straight-length liquid-height]
  (let [cylinder-volume
        (* straight-length
           (- (* radius radius (Math/acos (/ (- radius liquid-height)
                                             radius)))
              (* (- radius liquid-height)
                 (Math/sqrt (* liquid-height
                               (- (* 2.0 radius) liquid-height))))))
        end-caps-volume
        (* Math/PI liquid-height liquid-height
           (- radius (/ liquid-height 3.0)))]
    (+ cylinder-volume end-caps-volume)))

;; ## 3. Liquid Surface Area
;;
;; The free surface of the liquid consists of two parts:
;; - a rectangular area within the cylindrical section;
;; - two semicircular areas within the hemispherical end caps.
;;
;; The total liquid surface area [$m^2$] is:
;;
;; $$S(h) = 2L\sqrt{2rh-h^2} + \pi(2rh-h^2)$$
;;
;; Valid for $0 \leqslant h \leqslant 2r$.

(defn liquid-surface-area
  [radius straight-length liquid-height]
  (let [a-squared (* liquid-height (- (* 2.0 radius) liquid-height))
        cylinder-area (* 2.0 straight-length (Math/sqrt a-squared))
        end-caps-area (* Math/PI a-squared)]
    (+ cylinder-area end-caps-area)))

;; ## 4. Liquid Height from Mass
;;
;; For a given liquid mass and density, the liquid volume [$m^3$] is:
;;
;; $$V = \frac{m}{\rho}$$
;;
;; where:
;; - $m$ is the liquid mass [$kg$]
;; - $\rho$ is the liquid density [$kg/m^3$]
;;
;; The liquid height is calculated by solving:
;;
;; $$V(h) = \frac{m}{\rho}$$
;;
;; The bisection method is used over the interval
;; $0 \leqslant h \leqslant 2r$.

(defn liquid-height-from-mass
  [radius straight-length mass density]
  (when-not (pos? density)
    (throw (ex-info "Density must be positive."
                    {:density density})))
  (let [target-volume (/ mass density)
        total-volume (amp-volume radius straight-length)]
    (when-not (<= 0.0 target-volume total-volume)
      (throw (ex-info "Liquid volume is outside ampoule capacity."
                      {:volume target-volume
                       :capacity total-volume})))
    (cond
      (zero? target-volume) 0.0

      (= target-volume total-volume)
      (* 2.0 radius)

      :else

      (loop [low 0.0
             high (* 2.0 radius)
             steps 0]
        (let [mid (/ (+ low high) 2.0)]
          (if (= steps 60)
            mid
            (if (< (filled-volume radius straight-length mid)
                   target-volume)
              (recur mid high (inc steps))
              (recur low mid (inc steps)))))))))

;; ## 5. Saturated Vapor Pressure of Antimony
;;
;; The saturated vapor pressure depends on temperature and on the
;; phase of antimony.
;;
;; ### 5.1. Liquid Antimony
;;
;; $$\log_{10}P = 2.26041 - \frac{4475.449}{T-152.352}$$
;;
;; $T$ is in $K$, $P$ is in $bar$.
;; Experimental range: 785-1100 °C.
;;
;; Source: Jellinek and Wannow (1935).
;; [NIST Chemistry WebBook](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440360&Mask=4)

(defn sb-liquid-vapor-pressure
  [temp-c]
  (let [temp-k (+ temp-c 273.15)
        pressure-bar
        (Math/pow 10.0
                  (- 2.26041
                     (/ 4475.449 (- temp-k 152.352))))]
    (* 1.0e5 pressure-bar)))

;; ### 5.2. Solid Antimony
;;
;; $$\log_{10}P = 10.571 - \frac{10300}{T}$$
;;
;; $T$ is in $K$, $P$ is in $Torr$.
;; Experimental range: 420-550 °C.
;;
;; Source: Rosenblatt and Birchenall (1961).
;; [Journal of Chemical Physics](https://doi.org/10.1063/1.1701217)

(defn sb-solid-vapor-pressure
  [temp-c]
  (let [temp-k (+ temp-c 273.15)
        pressure-torr
        (Math/pow 10.0
                  (- 10.571 (/ 10300.0 temp-k)))]
    (* 133.32236842105263 pressure-torr)))

;; ### 5.3. Saturated Vapor Pressure
;;
;; Temperature range: 450-800 °C.
;; Logarithmic pressure scale.
;; Values outside the experimental ranges are extrapolated.

(defn pa-to-torr
  [pressure-pa]
  (/ pressure-pa 133.32236842105263))

(let [;; Physical constants and plotting range.
      melting-point 630.63
      x-domain [450 800]
      y-domain [1.0e-4 5.0]

      ;; Major and minor tick positions.
      x-major-values (vec (range 450 801 50))
      x-minor-values (vec (remove #(zero? (mod % 50))
                                  (range 450 801 10)))

      y-major-values (vec (filter #(<= % (second y-domain))
                                  [1.0e-4 1.0e-3 1.0e-2 1.0e-1 1.0 10.0]))

      y-minor-values
      (vec
       (for [exp (range -4 1)
             m (range 2 10)
             :let [y (* m (Math/pow 10.0 exp))]
             :when (<= y (second y-domain))]
         y))

      ;; Main data: solid and liquid vapor pressure curves.
      points
      (vec
       (concat
        (for [temp-c (concat (range 450 631 5) [melting-point])]
          {:temp-c temp-c
           :pressure-torr
           (pa-to-torr (sb-solid-vapor-pressure temp-c))
           :phase "Solid"})

        (for [temp-c (concat [melting-point] (range 635 801 5))]
          {:temp-c temp-c
           :pressure-torr
           (pa-to-torr (sb-liquid-vapor-pressure temp-c))
           :phase "Liquid"})))

      ;; Special annotation: melting point.
      melting-data
      [{:temp-c melting-point
        :pressure-torr
        (pa-to-torr (sb-liquid-vapor-pressure melting-point))
        :label "melting point"}]

      ;; Plot appearance.
      curve-width 3

      axis-label-font-size 15
      axis-title-font-size 17
      legend-font-size 15
      melting-label-font-size 14

      ;; Reusable axis definitions.
      x-axis
      {:title "Temperature (°C)"
       :values x-major-values
       :grid false
       :labelFontSize axis-label-font-size
       :titleFontSize axis-title-font-size}

      y-axis
      {:title "Vapor pressure (Torr)"
       :values y-major-values
       :grid false
       :labelFontSize axis-label-font-size
       :titleFontSize axis-title-font-size}

      ;; Reusable legend definition.
      phase-legend
      {:title nil
       :orient "top-left"
       :offset 8
       :fillColor "white"
       :strokeColor "lightgray"
       :padding 6
       :labelFontSize legend-font-size
       :symbolSize 180}

      ;; Helper layers.
      data-lines-layer
      {:data {:values points}
       :mark {:type "line"
              :strokeWidth curve-width}
       :encoding
       {:x {:field "temp-c"
            :type "quantitative"
            :scale {:domain x-domain}
            :axis x-axis}

        :y {:field "pressure-torr"
            :type "quantitative"
            :scale {:type "log"
                    :domain y-domain}
            :axis y-axis}

        :color {:field "phase"
                :type "nominal"
                :legend phase-legend}}}

      melting-point-layer
      {:data {:values melting-data}
       :mark {:type "circle"
              :size 100
              :color "black"}
       :encoding
       {:x {:field "temp-c"
            :type "quantitative"}
        :y {:field "pressure-torr"
            :type "quantitative"
            :scale {:type "log"}}}}

      melting-label-layer
      {:data {:values melting-data}
       :mark {:type "text"
              :color "black"
              :fontSize melting-label-font-size
              :align "left"
              :dx 10
              :dy 5}
       :encoding
       {:x {:field "temp-c"
            :type "quantitative"}
        :y {:field "pressure-torr"
            :type "quantitative"
            :scale {:type "log"}}
        :text {:field "label"
               :type "nominal"}}}]

  (kind/vega-lite
   (merge plot/theme
          {:title "Saturated vapor pressure of antimony"

           :layer
           [(plot/grid-layer :x x-minor-values
                             {:domain x-domain} false)
            (plot/grid-layer :x x-major-values
                             {:domain x-domain} true)
            (plot/grid-layer :y y-minor-values
                             {:type "log" :domain y-domain} false)
            (plot/grid-layer :y y-major-values
                             {:type "log" :domain y-domain} true)

            data-lines-layer
            melting-point-layer
            melting-label-layer]})))

;; ## 6. Mean Free Path and Knunsen Number
;;
;; The vapor is assumed to consist entirely os $Sb_4$.
;; Molecules are approximated as hard spheres.
;;
;; $$\lambda = \frac{k_B T}{\sqrt{2}\pi d^2 P}$$
;;
;; $T$ is in $K$, $d$ is in $m$, $P$ is in $Pa$.
;; The mean free path is returned in $m$.
;;
;; $$\mathrm{Kn} = \frac{\lambda}{D}$$
;;
;; $D$ is the internal tube diameter in m.

(def boltzmann-const 1.380649e-23)

(defn sb4-mean-free-path
  [temp-c pressure-pa collision-diameter-m]
  (let [temp-k (+ temp-c 273.15)]
    (when-not (and (pos? temp-k)
                   (pos? pressure-pa)
                   (pos? collision-diameter-m))
      (throw (ex-info "Invalid gas parameters."
                      {:temp-c temp-c
                       :pressure-pa pressure-pa
                       :collision-diameter-m collision-diameter-m})))
    (/ (* boltzmann-const temp-k)
       (* (Math/sqrt 2.0)
          Math/PI
          collision-diameter-m
          collision-diameter-m
          pressure-pa))))

(defn knudsen-number
  [mean-free-path-m tube-diameter-m]
  (when-not (pos? tube-diameter-m)
    (throw (ex-info "Tube diameter must be positive."
                    {:tube-diameter-m tube-diameter-m})))
  (/ mean-free-path-m tube-diameter-m))

;; ### 6.1. Knudsen Number at the Source
;;
;; The local vapor pressure is approximated by the saturation
;; pressure of liquid antimony.
;; The collision diameter is assumed to be 0.6 $nm$.
;; The pressure correlation is extrapolated below its
;; experimental temperature range.

(let [collision-diameter 0.6e-9
      tube-diameters-mm [8 10 16 20]

      x-domain [650 750]
      y-domain [0.001 0.1]

      x-major-values [650 700 750]
      x-minor-values
      (vec (remove #(zero? (mod % 50))
                   (range 650 751 10)))

      y-major-values [0.001 0.01 0.1]
      y-minor-values
      (vec
       (for [exp (range -3 -1)
             m (range 2 10)]
         (* m (Math/pow 10.0 exp))))

      ;; Calculate $Kn$ for each temperature and tube diameter.
      points
      (vec
       (for [diameter-mm tube-diameters-mm
             temp-c (range 650 751 2)
             :let [pressure
                   (sb-liquid-vapor-pressure temp-c)

                   mean-free-path
                   (sb4-mean-free-path
                    temp-c pressure collision-diameter)]]
         {:temp-c temp-c
          :kn (knudsen-number mean-free-path
                              (* diameter-mm 1.0e-3))
          :diameter (str diameter-mm " mm")}))

      ;; Main curves
      curves-layer
      {:data {:values points}
       :mark {:type "line"
              :strokeWidth 3}
       :encoding
       {:x {:field "temp-c"
            :type "quantitative"
            :scale {:domain x-domain}
            :axis {:title "Temperature (°C)"
                   :values x-major-values}}
        :y {:field "kn"
            :type "quantitative"
            :scale {:type "log"
                    :domain y-domain}
            :axis {:title "Knudsen number"
                   :values y-major-values}}

        :color {:field "diameter"
                :type "nominal"
                :sort ["8 mm" "10 mm" "16 mm" "20 mm"]
                :legend {:title nil}}}}]

  (kind/vega-lite
   (merge plot/theme
          {:title "Knudsen number at the source"

           :layer
           [(plot/grid-layer :x x-minor-values
                             {:domain x-domain} false)
            (plot/grid-layer :x x-major-values
                             {:domain x-domain} true)
            (plot/grid-layer :y y-minor-values
                             {:type "log" :domain y-domain} false)
            (plot/grid-layer :y y-major-values
                             {:type "log" :domain y-domain} true)
            curves-layer]})))

;; ### 6.2. Knudsen Number at the Receiver
;;
;; The local vapor pressure is approximated by
;; the saturation pressure of solid antimony.
;;
;; The collision diameter is assumed to be 0.6 nm.
;; The temperature range is 450-550 °C.

(let [collision-diameter-m 0.6e-9
      tube-diameters-mm [8 10 16 20]

      x-domain [450 550]
      y-domain [0.1 100.0]

      ;; Major and minor tick positions.
      x-major-values [450 500 550]
      x-minor-values
      (vec (remove #(zero? (mod % 50))
                   (range 450 551 10)))

      y-major-values [0.1 1.0 10.0 100.0]
      y-minor-values
      (vec
       (for [exp (range -1 2)
             m (range 2 10)]
         (* m (Math/pow 10.0 exp))))

      ;; Calculate Kn for each temperature and tube diameter.
      points
      (vec
       (for [diameter-mm tube-diameters-mm
             temp-c (range 450 551 2)
             :let [pressure-pa
                   (sb-solid-vapor-pressure temp-c)

                   mean-free-path-m
                   (sb4-mean-free-path
                    temp-c pressure-pa collision-diameter-m)]]

         {:temp-c temp-c
          :kn (knudsen-number
               mean-free-path-m
               (* diameter-mm 1.0e-3))
          :diameter (str diameter-mm " mm")}))

      ;; Main curves.
      curves-layer
      {:data {:values points}
       :mark {:type "line"
              :strokeWidth 3}
       :encoding
       {:x {:field "temp-c"
            :type "quantitative"
            :scale {:domain x-domain}
            :axis {:title "Temperature (°C)"
                   :values x-major-values}}

        :y {:field "kn"
            :type "quantitative"
            :scale {:type "log"
                    :domain y-domain}
            :axis {:title "Knudsen number"
                   :values y-major-values}}

        :color {:field "diameter"
                :type "nominal"
                :sort ["8 mm" "10 mm" "16 mm" "20 mm"]
                :legend {:title nil}}}}]

  (kind/vega-lite
   (merge plot/theme
          {:title "Knudsen number at the receiver"

           :layer
           [(plot/grid-layer :x x-minor-values
                             {:domain x-domain} false)
            (plot/grid-layer :x x-major-values
                             {:domain x-domain} true)
            (plot/grid-layer :y y-minor-values
                             {:type "log" :domain y-domain} false)
            (plot/grid-layer :y y-major-values
                             {:type "log" :domain y-domain} true)

            curves-layer]})))

^:kindly/hide-code
(kind/hiccup
 [:script
  "
  document.querySelectorAll('pre > code.source-clojure')
    .forEach(function (code) {
      const pre = code.parentElement;

      if (pre.closest('details')) return;

      const details = document.createElement('details');
      //details.open = true;
      const summary = document.createElement('summary');

      summary.textContent = 'Code';
      summary.style.cursor = 'pointer';

      pre.replaceWith(details);
      details.append(summary, pre);
    });
  "])
