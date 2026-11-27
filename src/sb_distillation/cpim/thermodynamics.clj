(ns sb-distillation.cpim.thermodynamics
  "CPIM compatibility names for the shared thermodynamic correlations."
  (:require [sb-distillation.physics.thermodynamics :as thermo]))

(def sb-liquid-vapor-pressure thermo/sb-liquid-vapor-pressure)
(def sb-solid-vapor-pressure thermo/sb-solid-vapor-pressure)
