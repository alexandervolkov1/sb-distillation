^:kindly/hide-code
(ns sb-distillation.analysis
  (:require
   [scicloj.kindly.v4.kind :as kind]
   [sb-distillation.diagrams.apparatus :as apparatus]
   [sb-distillation.diagrams.boundary :as boundary]
   [sb-distillation.experiments.time-grid :as time-grid]
   [sb-distillation.visualization.time-maps :as time-maps]))

;; # Численное моделирование переноса сурьмы

;; ## Постановка задачи

^:kindly/hide-code
(kind/hiccup
 (apparatus/boundary-conditions-svg))

;; ## Граничные условия

^:kindly/hide-code
(kind/hiccup
 (boundary/boundary-summary-svg))

;; ## Карты времени перегонки сурьмы
;;
;; Цвет показывает, сколько суток требуется для перегонки 500 г сурьмы.
;; Цветовая шкала общая для всех карт. Толстая изолиния соответствует одним суткам.
;; При нулевой подаче H₂ обе панели показывают одно аналитическое решение
;; для чистого пара Sb₄ с выходом в вакуум.

^:kindly/hide-code
(kind/hiccup
 [:table.table
  [:thead [:tr [:th "Параметр"] [:th "Значение"]]]
  [:tbody
   [:tr [:td "Температура источника"] [:td "680–720 °C, шаг 1 °C"]]
   [:tr [:td "Внутренний диаметр трубки"] [:td "6–12 мм, шаг 0,1 мм"]]
   [:tr [:td "Длина трубки"] [:td "0,500 м"]]
   [:tr [:td "Масса сурьмы"] [:td "500 г"]]
   [:tr [:td "Подача H₂ на картах"] [:td "0, 100, 200 и 300 мл/ч"]]
   [:tr [:td "Опорное состояние для объёмного расхода"] [:td "293,15 K, 101 325 Па"]]
   [:tr [:td "Перепад на ловушке"] [:td "0 Па"]]
   [:tr [:td "Измерения вакуумметра"]
    [:td "Для 100, 150, 200, 250 и 300 мл/ч соответственно 0,02010; 0,02870; 0,03680; 0,04465; 0,05110 Торр. Использованы без газовой поправки."]]
   [:tr [:td "Допуски"] [:td "RK4: 10⁻⁸ Па и 10⁻⁶ относительно. Пристрелка: 10⁻³ Па."]]]])

^:kindly/hide-code
(let [study time-grid/default-study
      {:keys [rows]} (time-grid/cached-grid study "target/distillation-maps")]
  (time-maps/maps rows study))

;; В пределах рассчитанной сетки время уменьшается при увеличении диаметра
;; или температуры при фиксированной подаче водорода.
;; Сравнение моделей при ненулевой подаче приведено в
;; [отчёте о сравнении CPIM и DGM](https://github.com/alexandervolkov1/sb-distillation/blob/2b89fd01e6a67254aa34d70a01836732e084ad05/scripts/cpim-verification/COMPARISON.md).
