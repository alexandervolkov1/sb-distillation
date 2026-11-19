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

;; Basic consistency check:

(let [radius 0.03
      straight-length 0.18
      total-volume (amp-volume radius straight-length)
      volume (fn [liquid-height]
               (filled-volume radius straight-length liquid-height))
      eps 1.0e-12]
  (assert (< (Math/abs (volume 0.0)) eps))
  (assert (< (Math/abs (- (volume radius) (/ total-volume 2.0))) eps))
  (assert (< (Math/abs (- (volume (* 2.0 radius)) total-volume)) eps))
  :ok)

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
;; Valid for $0 \leqslant h \leqslant 2r$

(defn liquid-surface-area
  [radius straight-length liquid-height]
  (let [a (* liquid-height (- (* 2.0 radius) liquid-height))
        cylinder-area (* 2.0 straight-length (Math/sqrt a))
        end-caps-area (* Math/PI a)]
    (+ cylinder-area end-caps-area)))

;; Basic consistency check:

(let [radius 0.03
      straight-length 0.18
      liquid-area (fn [liquid-height]
                    (liquid-surface-area radius
                                         straight-length
                                         liquid-height))
      max-area (+ (* 2.0 radius straight-length)
                  (* Math/PI radius radius))
      eps 1.0e-12]
  (assert (< (Math/abs (liquid-area 0.0)) eps))
  (assert (< (Math/abs (- (liquid-area radius) max-area)) eps))
  (assert (< (Math/abs (liquid-area (* 2.0 radius))) eps))
  :ok)
