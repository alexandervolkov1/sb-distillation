(ns sb-distillation.diagrams.boundary)

(defn boundary-summary-svg
  "Black-and-white diagram of the CPIM
  boundary-value problem."
  []
  [:svg {:xmlns "http://www.w3.org/2000/svg"
         :viewBox "0 0 720 300"
         :width "100%"
         :font-family "Old Standard, serif"
         :font-size 18
         :fill "black"}

   [:defs
    [:marker {:id "boundary-arrow"
              :viewBox "0 0 10 8"
              :markerWidth 10
              :markerHeight 8
              :markerUnits "userSpaceOnUse"
              :refX 10
              :refY 4
              :orient "auto"}
     [:path {:d "M0 0 L10 4 L0 8 Z"
             :fill "black"}]]]

   [:rect {:width 720
           :height 300
           :fill "white"}]

   ;; Tube walls
   [:g {:stroke "black"
        :stroke-width 1.6}
    (for [y [50 135]]
      [:line {:x1 60 :y1 y
              :x2 660 :y2 y}])]

   ;; Temperature and boundary labels
   [:text {:x 60 :y 24} "Вход"]

   [:text {:x 660 :y 24
           :text-anchor "end"}
    "Выход"]

   [:text {:x 360 :y 39
           :text-anchor "middle"}
    "T = Tₛ = const"]

   ;; Species fluxes
   [:g {:stroke "black"
        :stroke-width 1.6
        :marker-end "url(#boundary-arrow)"}

    [:line {:x1 85 :y1 84
            :x2 635 :y2 84}]

    [:line {:x1 85 :y1 121
            :x2 635 :y2 121
            :stroke-dasharray "7 4"}]]

   [:g {:text-anchor "middle"}
    [:text {:x 360 :y 76} "N₁ — искомый поток"]
    [:text {:x 360 :y 113} "N₂ — заданный поток"]]

   ;; Coordinate references
   (into
    [:g {:stroke "black"
         :stroke-width 1
         :stroke-dasharray "4 4"}]
    (for [x [60 660]]
      [:line {:x1 x :y1 50
              :x2 x :y2 155}]))

   [:g {:text-anchor "middle"}
    [:text {:x 60 :y 175} "z = 0"]
    [:text {:x 660 :y 175} "z = L"]]

   ;; Inlet boundary
   [:text {:x 60 :y 202}
    "p₁(0) = p₁,sat(Tₛ)"]

   [:text {:x 60 :y 227}
    "p₂(0) — неизвестно"]

   ;; Outlet boundary
   [:g {:text-anchor "end"}
    [:text {:x 660 :y 202}
     "p₁(L) = 0"]

    [:text {:x 660 :y 227}
     "p₂(L) = p"
     [:tspan {:baseline-shift "sub"
              :font-size 13} "g"]
     " + Δp"
     [:tspan {:baseline-shift "sub"
              :font-size 13} "trap"]]]

   ;; Backward integration
   [:line {:x1 660 :y1 257
           :x2 60 :y2 257
           :stroke "black"
           :stroke-width 1.6
           :marker-end "url(#boundary-arrow)"}]

   [:text {:x 360 :y 281
           :text-anchor "middle"}
    "Направление интегрирования"]])
