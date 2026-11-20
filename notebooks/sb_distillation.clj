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
                :legend phase-legend}
        :tooltip
        [{:field "temp-c"
          :type "quantitative"
          :title "Temperature (°C)"
          :format ".1f"}

         {:field "pressure-torr"
          :type "quantitative"
          :title "Vapor pressure (Torr)"
          :format ".4g"}

         {:field "phase"
          :type "nominal"
          :title "Phase"}]}}

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
;; $D$ is the internal tube diameter in $m$.

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

;; ## 7. Sb4 Mass Flow Through the Tube
;;
;; The vapor is assumed to consist entirely of Sb4.
;; The connecting tube is straight, cylindrical, and isothermal.
;; Its wall temperature equals the source temperature.
;; No condensation occurs inside the tube.
;;
;; The inlet pressure is approximated by the saturated vapor
;; pressure of liquid antimony at the source temperature.
;; The outlet is treated as an ideal vapor sink ($p_2 = 0$).
;;
;; ### 7.1. Tube Conductance
;;
;; The transitional-flow conductance is approximated
;; using the Knudsen formula:
;;
;; $$C = C_{\mathrm{vis}} + Z C_{\mathrm{mol}}$$
;;
;; where:
;; - $C$ is the total conductance [$m^3/s$]
;; - $C_{\mathrm{vis}}$ is the viscous conductance [$m^3/s$]
;; - $C_{\mathrm{mol}}$ is the molecular conductance [$m^3/s$]
;; - $Z$ is the dimensionless Knudsen correction factor
;;
;; The viscous conductance is:
;;
;; $$C_{\mathrm{vis}} =
;; \frac{\pi D^4 \bar{p}}{128\eta L}$$
;;
;; The molecular conductance is:
;;
;; $$C_{\mathrm{mol}} =
;; \frac{\pi D^3}{12L}\sqrt{\frac{8RT}{\pi M}}$$
;;
;; where:
;; - $D$ is the internal tube diameter [$m$]
;; - $L$ is the tube length [$m$]
;; - $\bar{p}=(p_1+p_2)/2$ is the mean pressure [$Pa$]
;; - $\eta$ is the dynamic viscosity [$Pa \cdot s$]
;; - $T$ is the absolute temperature [$K$]
;; - $R$ is the universal gas constant [$J/(mol \cdot K)$]
;; - $M$ is the molar mass of Sb4 [$kg/mol$]
;;
;; The Knudsen correction factor is:
;;
;; $$Z = \frac{1+2.507r/\bar{\lambda}}
;; {1+3.095r/\bar{\lambda}}$$
;;
;; where:
;; - $r=D/2$ is the tube radius [$m$]
;; - $\bar{\lambda}$ is the mean free path at $\bar{p}$ [$m$]
;;
;; ### 7.2. Gas Viscosity
;;
;; The dynamic viscosity is estimated using
;; the hard-sphere gas model:
;;
;; $$\eta = \frac{5}{16d^2}
;; \sqrt{\frac{mk_BT}{\pi}}$$
;;
;; where:
;; - $d$ is the molecular collision diameter [$m$]
;; - $m=M/N_A$ is the mass of one Sb4 molecule [$kg$]
;; - $k_B$ is the Boltzmann constant [$J/K$]
;; - $N_A$ is the Avogadro constant [$mol^{-1}$]
;;
;; The collision diameter is assumed to be 0.6 nm.

(def gas-constant 8.314462618)
(def avogadro-constant 6.02214076e23)
(def sb4-molar-mass 0.48704)

(defn sb4-hard-sphere-viscosity
  [temp-c collision-diameter-m]
  (let [temp-k (+ temp-c 273.15)
        molecule-mass (/ sb4-molar-mass avogadro-constant)]
    (when-not (and (pos? temp-k)
                   (pos? collision-diameter-m))
      (throw (ex-info "Invalid gas parameters."
                      {:temp-c temp-c
                       :collision-diameter-m collision-diameter-m})))
    (/ (* 5.0
          (Math/sqrt
           (/ (* molecule-mass boltzmann-const temp-k)
              Math/PI)))
       (* 16.0 collision-diameter-m collision-diameter-m))))

;; ### 7.3. Mass Flow Rate
;;
;; The mass flow rate is:
;;
;; $$\dot{m} = \frac{M}{RT} C (p_1-p_2)$$
;;
;; where:
;; - $\dot{m}$ is the mass flow rate [$kg/s$]
;; - $p_1$ is the inlet pressure [$Pa$]
;; - $p_2$ is the outlet pressure [$Pa$]
;;
;; The function returns the mass flow rate in kg/s.

(defn sb4-tube-mass-flow
  [temp-c tube-diameter-m tube-length-m collision-diameter-m]
  (when-not (and (pos? tube-diameter-m)
                 (pos? tube-length-m))
    (throw (ex-info "Invalid tube geometry."
                    {:tube-diameter-m tube-diameter-m
                     :tube-length-m tube-length-m})))

  (let [temp-k (+ temp-c 273.15)

        ;; Boundary pressures.
        p1 (sb-liquid-vapor-pressure temp-c)
        p2 0.0
        p-mean (/ (+ p1 p2) 2.0)

        ;; Gas properties.
        viscosity
        (sb4-hard-sphere-viscosity
         temp-c collision-diameter-m)

        mean-free-path
        (sb4-mean-free-path
         temp-c p-mean collision-diameter-m)

        ;; Tube geometry.
        radius (/ tube-diameter-m 2.0)

        ;; Molecular conductance.
        c-molecular
        (* (/ (* Math/PI
                 (Math/pow tube-diameter-m 3.0))
              (* 12.0 tube-length-m))
           (Math/sqrt
            (/ (* 8.0 gas-constant temp-k)
               (* Math/PI sb4-molar-mass))))

        ;; Viscous conductance.
        c-viscous
        (/ (* Math/PI
              (Math/pow tube-diameter-m 4.0)
              p-mean)
           (* 128.0 viscosity tube-length-m))

        ;; Knudsen correction factor.
        ratio (/ radius mean-free-path)

        z (/ (+ 1.0 (* 2.507 ratio))
             (+ 1.0 (* 3.095 ratio)))

        ;; Total conductance.
        conductance
        (+ c-viscous (* z c-molecular))]

    (* (/ sb4-molar-mass (* gas-constant temp-k))
       conductance
       (- p1 p2))))

(defn sb4-distillation-days
  [mass-kg temp-c tube-diameter-m tube-length-m collision-diameter-m]
  (when-not (pos? mass-kg)
    (throw (ex-info "Mass must be positive."
                    {:mass-kg mass-kg})))
  (/ mass-kg
     (* 86400.0
        (sb4-tube-mass-flow
         temp-c
         tube-diameter-m
         tube-length-m
         collision-diameter-m))))

;; ### 7.4. Mass Flow versus Source Temperature
;;
;; Tube length: 0.5 m.
;; Tube diameters: 8, 10, 16, and 20 mm.
;; Collision diameter: 0.6 nm.
;;
;; The calculated mass flow is converted to g/h.
;;
;; These are preliminary estimates based on a mean-pressure
;; Knudsen conductance approximation and ideal boundary conditions.
;; The liquid vapor pressure correlation is extrapolated
;; below its experimental range of 785-1100 °C.

(let [;; Model parameters.
      tube-length-m 0.5
      collision-diameter-m 0.6e-9
      tube-diameters-mm [8 10 16 20]

      ;; Plot ranges.
      x-domain [650 750]
      y-domain [0.1 1000.0]

      ;; Major and minor tick positions.
      x-major-values [650 700 750]

      x-minor-values
      (vec
       (remove #(zero? (mod % 50))
               (range 650 751 10)))

      y-major-values
      [0.1 1.0 10.0 100.0 1000.0]

      y-minor-values
      (vec
       (for [exp (range -1 3)
             m (range 2 10)]
         (* m (Math/pow 10.0 exp))))

      ;; Mass flow values in g/h.
      points
      (vec
       (for [diameter-mm tube-diameters-mm
             temp-c (range 650 751 2)]

         {:temp-c temp-c
          :mass-flow-g-h
          (* 3.6e6
             (sb4-tube-mass-flow
              temp-c
              (* diameter-mm 1.0e-3)
              tube-length-m
              collision-diameter-m))
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
            :axis {:title "Source temperature (°C)"
                   :values x-major-values}}

        :y {:field "mass-flow-g-h"
            :type "quantitative"
            :scale {:type "log"
                    :domain y-domain}
            :axis {:title "Mass flow (g/h)"
                   :values y-major-values}}

        :color {:field "diameter"
                :type "nominal"
                :sort ["8 mm" "10 mm" "16 mm" "20 mm"]
                :legend {:title nil}}

        :tooltip
        [{:field "temp-c"
          :type "quantitative"
          :title "Temperature (°C)"
          :format ".0f"}

         {:field "diameter"
          :type "nominal"
          :title "Tube diameter"}

         {:field "mass-flow-g-h"
          :type "quantitative"
          :title "Mass flow (g/h)"
          :format ".3f"}]}}]

  (kind/vega-lite
   (merge plot/theme
          {:title "Sb4 mass flow through the tube"

           :layer
           [(plot/grid-layer :x x-minor-values
                             {:domain x-domain} false)
            (plot/grid-layer :x x-major-values
                             {:domain x-domain} true)
            (plot/grid-layer :y y-minor-values
                             {:type "log"
                              :domain y-domain} false)
            (plot/grid-layer :y y-major-values
                             {:type "log"
                              :domain y-domain} true)

            curves-layer]})))

;; ### 7.5. Mass Flow versus Tube Diameter
;;
;; Source temperatures: 650-750 °C, in steps of 20 °C.
;; Tube diameter: 8-20 mm.
;; Tube length: 0.5 m.
;; Collision diameter: 0.6 nm.
;;
;; The calculated mass flow is converted to g/h.

(let [;; Model parameters.
      tube-length-m 0.5
      collision-diameter-m 0.6e-9
      temperatures-c (vec (range 650 751 20))

      ;; Plot ranges.
      x-domain [8 20]
      y-domain [0.1 1000.0]

      ;; Major and minor tick positions.
      x-major-values (vec (range 8 21 2))
      x-minor-values (vec (range 9 20 2))

      y-major-values
      [0.1 1.0 10.0 100.0 1000.0]

      y-minor-values
      (vec
       (for [exp (range -1 3)
             m (range 2 10)]
         (* m (Math/pow 10.0 exp))))

      ;; Mass flow values in g/h.
      points
      (vec
       (for [temp-c temperatures-c
             diameter-step (range 80 201)
             :let [diameter-mm (/ diameter-step 10.0)]]

         {:diameter-mm diameter-mm
          :mass-flow-g-h
          (* 3.6e6
             (sb4-tube-mass-flow
              temp-c
              (* diameter-mm 1.0e-3)
              tube-length-m
              collision-diameter-m))
          :temperature (str temp-c " °C")}))

      ;; Main curves.
      curves-layer
      {:data {:values points}
       :mark {:type "line"
              :strokeWidth 3}
       :encoding
       {:x {:field "diameter-mm"
            :type "quantitative"
            :scale {:domain x-domain}
            :axis {:title "Tube diameter (mm)"
                   :values x-major-values}}

        :y {:field "mass-flow-g-h"
            :type "quantitative"
            :scale {:type "log"
                    :domain y-domain}
            :axis {:title "Mass flow (g/h)"
                   :values y-major-values}}

        :color {:field "temperature"
                :type "nominal"
                :sort (mapv #(str % " °C") temperatures-c)
                :legend {:title nil}}

        :tooltip
        [{:field "diameter-mm"
          :type "quantitative"
          :title "Tube diameter (mm)"
          :format ".1f"}

         {:field "temperature"
          :type "nominal"
          :title "Source temperature"}

         {:field "mass-flow-g-h"
          :type "quantitative"
          :title "Mass flow (g/h)"
          :format ".3f"}]}}]

  (kind/vega-lite
   (merge plot/theme
          {:title "Sb4 mass flow versus tube diameter"

           :layer
           [(plot/grid-layer :x x-minor-values
                             {:domain x-domain} false)
            (plot/grid-layer :x x-major-values
                             {:domain x-domain} true)
            (plot/grid-layer :y y-minor-values
                             {:type "log"
                              :domain y-domain} false)
            (plot/grid-layer :y y-major-values
                             {:type "log"
                              :domain y-domain} true)

            curves-layer]})))

;; ### 7.6. Distillation Time Contours
;;
;; Estimated time to distill 0.5 kg of antimony.
;; The tube limits the flow; its wall is at source temperature.
;; No carrier gas, no tube condensation, ideal vapor sink.
;; The flow is assumed constant during distillation.

(let [;; Physical parameters.
      mass-kg 0.5
      tube-length-m 0.5
      collision-diameter-m 0.6e-9

      ;; Plot ranges.
      x-domain [8.0 20.0]
      y-domain [650.0 750.0]

      ;; Major contours are labeled; minor contours are not.
      major-levels
      [0.25 0.5 0.75 1.0 2.0 3.0 5.0 7.0 10.0 14.0 20.0]

      minor-levels
      [4.0 6.0 8.0 9.0 12.0 16.0 18.0]

      ;; Only major grid lines are shown.
      x-major-values (vec (range 8 21 2))
      y-major-values (vec (range 650 751 20))

      ;; Duration in days.
      duration
      (fn [temp-c diameter-mm]
        (sb4-distillation-days
         mass-kg temp-c (* diameter-mm 1.0e-3)
         tube-length-m collision-diameter-m))

      ;; Fine rectangular mesh for the background color field.
      ;; Cell size: 0.2 mm by 1 °C.
      heatmap-data
      (vec
       (for [di (range 40 100)
             ti (range 650 750)
             :let [d0 (/ di 5.0)
                   d1 (/ (inc di) 5.0)
                   t0 (double ti)
                   t1 (inc t0)]]
         {:d0 d0
          :d1 (min (second x-domain) (+ d1 0.012))
          :t0 t0
          :t1 (min (second y-domain) (+ t1 0.12))
          :days (duration (/ (+ t0 t1) 2.0)
                          (/ (+ d0 d1) 2.0))}))

      ;; Logarithmic color scale for the filled field.
      fill-color-scale
      {:type "log"
       :domain [(duration 750.0 20.0)
                (duration 650.0 8.0)]
       :range ["#1a9850"
               "#91cf60"
               "#ffffbf"
               "#fc8d59"
               "#d73027"]}

      ;; Darker version of the same scale for contour lines and labels.
      line-color-scale
      {:type "log"
       :domain [(duration 750.0 20.0)
                (duration 650.0 8.0)]
       :range ["#0b5d2a"
               "#4f7f2a"
               "#b39b00"
               "#a64b1f"
               "#7f1d1d"]}

      ;; Format duration labels.
      format-days
      (fn [days]
        (if (< days 1.0)
          (str (int (* 24.0 days)) " h")
          (str (int days) " d")))

      ;; Solve t(T, D) = target-days for D by bisection.
      diameter-at-days
      (fn [temp-c target-days]
        (let [d-min (first x-domain)
              d-max (second x-domain)]
          (when (<= (duration temp-c d-max)
                    target-days
                    (duration temp-c d-min))
            (loop [low d-min
                   high d-max
                   steps 0]
              (let [mid (/ (+ low high) 2.0)]
                (if (= steps 45)
                  mid
                  (if (> (duration temp-c mid) target-days)
                    (recur mid high (inc steps))
                    (recur low mid (inc steps)))))))))

      temperatures-c
      (mapv #(/ % 2.0) (range 1300 1501))

      ;; Generate points for each contour.
      contour-points
      (fn [levels]
        (vec
         (for [days levels
               temp-c temperatures-c
               :let [diameter-mm
                     (diameter-at-days temp-c days)]
               :when diameter-mm]
           {:diameter-mm diameter-mm
            :temp-c temp-c
            :days days
            :label (format-days days)})))

      major-points (contour-points major-levels)
      minor-points (contour-points minor-levels)

;; Place horizontal labels along the plot diagonal.
      labels
      (vec
       (keep
        (fn [days]
          (let [curve
                (filterv #(= days (:days %)) major-points)]
            (when (seq curve)
              (let [best
                    (apply min-key
                           (fn [{:keys [diameter-mm temp-c]}]
                             (let [x (/ (- diameter-mm (first x-domain))
                                        (- (second x-domain) (first x-domain)))
                                   y (/ (- temp-c (first y-domain))
                                        (- (second y-domain) (first y-domain)))]
                               (Math/abs (- x y))))
                           curve)

                    half-width-mm
                    (+ 0.14 (* 0.06 (count (:label best))))

                    half-height-c 2.0]

                (assoc best
                       :d0 (- (:diameter-mm best) half-width-mm)
                       :d1 (+ (:diameter-mm best) half-width-mm)
                       :t0 (- (:temp-c best) half-height-c)
                       :t1 (+ (:temp-c best) half-height-c))))))
        major-levels))

;; Continuous background color field.
      fill-layer
      {:data {:values heatmap-data}
       :mark {:type "rect"
              :stroke nil
              :strokeWidth 0}
       :encoding
       {:x {:field "d0"
            :type "quantitative"
            :scale {:domain x-domain}
            :axis {:title "Tube diameter (mm)"
                   :values x-major-values}}
        :x2 {:field "d1"}

        :y {:field "t0"
            :type "quantitative"
            :scale {:domain y-domain}
            :axis {:title "Source temperature (°C)"
                   :values y-major-values}}
        :y2 {:field "t1"}

        :color {:field "days"
                :type "quantitative"
                :scale fill-color-scale
                :legend {:title "Days"
                         :orient "top-right"
                         :direction "vertical"
                         :offset 12
                         :gradientLength 200
                         :gradientThickness 14
                         :values [0.1 1.0 10.0]
                         :labelFontSize 13
                         :titleFontSize 14
                         :fillColor "white"
                         :strokeColor "#cccccc"
                         :cornerRadius 4
                         :padding 8
                         :zindex 1}}

        :tooltip
        [{:field "days"
          :type "quantitative"
          :title "Duration (days)"
          :format ".2f"}]}}

      ;; Reusable contour layer.
      contour-layer
      (fn [points stroke-width opacity]
        {:data {:values points}
         :mark {:type "line"
                :strokeWidth stroke-width
                :opacity opacity}
         :encoding
         {:x {:field "diameter-mm"
              :type "quantitative"
              :scale {:domain x-domain}}

          :y {:field "temp-c"
              :type "quantitative"
              :scale {:domain y-domain}}

          :detail {:field "days"
                   :type "nominal"}

          :order {:field "temp-c"
                  :type "quantitative"
                  :sort "ascending"}

          :color {:field "days"
                  :type "quantitative"
                  :scale line-color-scale
                  :legend nil}

          :tooltip
          [{:field "days"
            :type "quantitative"
            :title "Duration (days)"}
           {:field "diameter-mm"
            :type "quantitative"
            :title "Tube diameter (mm)"
            :format ".2f"}
           {:field "temp-c"
            :type "quantitative"
            :title "Source temperature (°C)"
            :format ".1f"}]}})

      ;; Horizontal contour labels.
      label-text-layer
      {:data {:values labels}
       :mark {:type "text"
              :color "black"
              :fontSize 16
              :fontWeight 900
              :align "center"
              :baseline "middle"
              :stroke "white"
              :strokeWidth 0.4}
       :encoding
       {:x {:field "diameter-mm"
            :type "quantitative"
            :scale {:domain x-domain}}

        :y {:field "temp-c"
            :type "quantitative"
            :scale {:domain y-domain}}

        :text {:field "label"
               :type "nominal"}}}]

  (kind/vega-lite
   (merge plot/theme
          {:width 1040
           :height 710

           :autosize {:type "fit"
                      :contains "padding"}

           :title "Estimated distillation time for 0.5 kg Sb"

           ;; Independent color scales for the filled field and contour lines.
           :resolve {:scale {:color "independent"}}

           :layer
           [fill-layer

            ;; Major grid lines only.
            (plot/grid-layer :x x-major-values
                             {:domain x-domain} true)
            (plot/grid-layer :y y-major-values
                             {:domain y-domain} true)

            (contour-layer minor-points 1.2 0.55)
            (contour-layer major-points 2.4 0.95)

            label-text-layer]})))

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
