^:kindly/hide-code
(ns sb-distillation.analysis
  (:require
   [scicloj.kindly.v4.kind :as kind]
   [sb-distillation.diagrams.apparatus :as apparatus]
   [sb-distillation.diagrams.boundary :as boundary]))

;; # Численное моделирование переноса сурьмы

;; ## Постановка задачи

^:kindly/hide-code
(kind/hiccup
 (apparatus/boundary-conditions-svg))

;; ## Граничные условия

^:kindly/hide-code
(kind/hiccup
 (boundary/boundary-summary-svg))
