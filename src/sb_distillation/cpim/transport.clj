(ns sb-distillation.cpim.transport
  "CPIM compatibility names for the shared transport coefficients."
  (:require [sb-distillation.physics.transport :as transport]))

(def default-parameters transport/default-parameters)
(def knudsen-diffusivity transport/knudsen-diffusivity)
(def binary-diffusivity-factor transport/binary-diffusivity-factor)
(def binary-diffusivity transport/binary-diffusivity)
(def sb4-viscosity transport/sb4-viscosity)
(def h2-viscosity transport/h2-viscosity)
(def wilke-interaction-factor transport/wilke-interaction-factor)
(def mixture-viscosity transport/mixture-viscosity)
