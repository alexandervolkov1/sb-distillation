(ns sb-distillation.diagrams.apparatus)

(defn boundary-conditions-svg
  "Black-and-white diagram of the CPIM apparatus
  and pressure boundary conditions."
  []
  [:svg {:xmlns "http://www.w3.org/2000/svg"
         :viewBox "0 0 720 465"
         :width "100%"
         :font-family "Old Standard, serif"
         :font-size 20
         :fill "black"}

   [:defs
    [:marker {:id "arrow"
              :markerWidth 8
              :markerHeight 8
              :refX 7
              :refY 3
              :orient "auto-start-reverse"}
     [:path {:d "M0 0 L7 3 L0 6 Z"}]]

    [:pattern {:id "hatch"
               :width 6
               :height 6
               :patternUnits "userSpaceOnUse"}
     [:path {:d "M0 6 L6 0"
             :stroke "black"
             :stroke-width 0.7}]]]

   [:rect {:width 720 :height 465 :fill "white"}]

   ;; Ampoules and connecting tube
   [:g {:fill "none" :stroke "black" :stroke-width 2}
    [:path {:d "M210 155 V125 Q210 90 175 90
                H60 Q25 90 25 125 V265
                Q25 300 60 300 H175 Q210 300 210 265 V205"}]

    [:path {:d "M510 155 V125 Q510 90 545 90
                H660 Q695 90 695 125 V265
                Q695 300 660 300 H545 Q510 300 510 265 V205"}]

    [:path {:d "M210 155 H510 M210 205 H510"}]]

   ;; Liquid antimony
   [:path {:d "M26 250 Q117 240 209 250 V265
               Q209 299 175 299 H60 Q26 299 26 265 Z"
           :fill "url(#hatch)"
           :stroke "black"}]

   ;; Deposited crystals
   (into
    [:g {:fill "white"
         :stroke "black"
         :stroke-width 1.5}]
    (for [[x h] [[545 37]
                 [580 50]
                 [615 34]
                 [650 46]]]
      [:polygon
       {:points
        (str x ",298 "
             (- x 10) "," (- 298 (/ h 2)) " "
             x "," (- 298 h) " "
             (+ x 10) "," (- 298 (/ h 2)))}]))

   ;; Gas inlet and outlet
   [:g {:stroke "black"
        :stroke-width 1.5
        :fill "none"
        :marker-end "url(#arrow)"}
    [:path {:d "M117 25 V90"}]
    [:path {:d "M602 90 V45 H690"}]]

   [:text {:x 132 :y 42} "H₂ от РРГ"]
   [:text {:x 650 :y 30 :text-anchor "middle"} "К ловушке"]

   ;; Species fluxes
   (into
    [:g {:stroke "black"
         :stroke-width 1.5
         :marker-end "url(#arrow)"}]
    (for [y [170 190]]
      [:line {:x1 235 :y1 y :x2 485 :y2 y}]))

   ;; Apparatus labels
   [:g {:text-anchor "middle"}
    [:text {:x 117 :y 120} "Источник"]
    [:text {:x 602 :y 120} "Приёмник"]
    [:text {:x 117 :y 225} "Sb (жидк.)"]
    [:text {:x 602 :y 225} "Sb (тв.)"]
    [:text {:x 360 :y 140} "N₁, N₂"]]

   ;; Coordinate references
   (into
    [:g {:stroke "black"
         :stroke-width 1
         :stroke-dasharray "4 4"}]
    (for [x [210 510]]
      [:line {:x1 x :y1 215 :x2 x :y2 350}]))

   ;; Tube length
   [:line {:x1 210 :y1 335
           :x2 510 :y2 335
           :stroke "black"
           :marker-start "url(#arrow)"
           :marker-end "url(#arrow)"}]

   [:g {:text-anchor "middle"}
    [:text {:x 360 :y 325} "L"]
    [:text {:x 210 :y 365} "z = 0"]
    [:text {:x 510 :y 365} "z = L"]

    ;; Boundary conditions
    [:text {:x 117 :y 392} "p₁(0) = p₁,sat(Tₛ)"]
    [:text {:x 117 :y 420} "p₂(0) = ?"]

    [:text {:x 602 :y 392} "p₁(L) = 0"]
    [:text {:x 602 :y 420}
     "p₂(L) = p"
     [:tspan {:baseline-shift "sub" :font-size 14} "g"]
     " + Q₂/C"
     [:tspan {:baseline-shift "sub" :font-size 14} "trap"]]]

   ;; Backward integration
   [:line {:x1 495 :y1 392
           :x2 225 :y2 392
           :stroke "black"
           :stroke-width 1.5
           :marker-end "url(#arrow)"}]

   [:text {:x 360 :y 420
           :text-anchor "middle"
           :font-size 16}
    "Интегрирование"]])
