(ns sb-distillation
  (:require [scicloj.kindly.v4.kind :as kind]))

;; # Vacuum Distillation of Antimony
;;
;; A preliminary model of pure antimony distillation
;; under vacuum, without carrier gas.
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
;; $0 \leqslant h \leqslant 2r$

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
    (loop [low 0.0
           high (* 2.0 radius)
           steps 0]
      (let [mid (/ (+ low high) 2.0)]
        (if (= steps 60)
          mid
          (if (< (filled-volume radius straight-length mid)
                 target-volume)
            (recur mid high (inc steps))
            (recur low mid (inc steps))))))))
