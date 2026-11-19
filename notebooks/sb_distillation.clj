(ns sb-distillation
  (:require [scicloj.kindly.v4.kind :as kind]))

;; # Vacuum Distillation of Antimony
;;
;; A first-approximation model of pure antimony distillation
;; under vacuum, without carrier gas.
;;
;; ## 1.Ampoule Volume
;;
;; Ampoule consist of a cylindrical section and two hemispherical caps.
;; The interal volume of ampoule [$m^3$] is:
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

;; ## 2.Liquid Fill Volume
;;
;; The ampoule is horisontal, and the liquid surface is assumed
;;  to be flat. The fill height $h$ [$m$] is measured from the bottom.
;;  The liquid volume consist of two parts:
;; - a circular segment along the cylindrical session;
;; - a spherical cap formed by the two hemispherical ends.
;;
;; Cross-sectional area of the liquid-filled segment in the cylindrical
;; section of the ampoule [$m^2$] is:
;;
;; $$A(h) = r^2\arccos\left(\frac{r-h}{r}\right)
;;    - (r-h)\sqrt{2rh-h^2}$$
;;
;; The total liquid volume [$m^3$] is:
;;
;; $$V(h) = L A(h) + \pi h^2\left(r-\frac{h}{3}\right)$$

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

;; The small test:

(let [radius 0.03
      straight-length 0.18]
  {:empty (filled-volume radius straight-length 0.0)
   :half-filled (filled-volume radius straight-length radius)
   :full (filled-volume radius straight-length (* 2.0 radius))
   :total-volume (amp-volume radius straight-length)})
