(ns sb-distillation.cpim.thermodynamics)

(defn sb-liquid-vapor-pressure
  "p_{1,sat} - Saturation vapor pressure
  over liquid antimony in Pa.
  Temperature in K.
  Antoine correlation with a fit range of 1058-1373 K.
  Values outside this range are extrapolated."
  [temp]
  (* 1.0e5
     (Math/pow
      10.0
      (- 2.26041
         (/ 4475.449 (- temp 152.352))))))
