(ns sb-distillation.report.diagram)

(defn- subscript [s]
  [:tspan {:baseline-shift "sub"
           :font-size "12"}
   s])

(defn apparatus-sketch []
  [:div {:style {:margin "1.5rem 0 2rem 0"}}
   [:svg {:viewBox "0 0 980 360"
          :width "100%"
          :style {:display "block"
                  :background "white"}}
    [:defs
     [:marker {:id "arrow"
               :markerWidth "10"
               :markerHeight "10"
               :refX "8"
               :refY "3"
               :orient "auto"
               :markerUnits "strokeWidth"}
      [:path {:d "M0,0 L0,6 L9,3 z"
              :fill "black"}]]]

    ;; Заголовок рисунка
    [:text {:x 490 :y 28
            :text-anchor "middle"
            :font-size 20
            :font-family "Old Standard"}
     "Схема установки и расчётной области"]

    ;; Основная трубка
    [:line {:x1 290 :y1 170 :x2 690 :y2 170
            :stroke "black" :stroke-width 4}]
    [:text {:x 490 :y 156
            :text-anchor "middle"
            :font-size 18
            :font-family "Old Standard"}
     "соединительная трубка"]

    ;; Размер L
    [:line {:x1 290 :y1 120 :x2 690 :y2 120
            :stroke "black" :stroke-width 1.5
            :marker-start "url(#arrow)"
            :marker-end "url(#arrow)"}]
    [:text {:x 490 :y 110
            :text-anchor "middle"
            :font-size 18
            :font-family "Old Standard"}
     "L"]

    ;; Размер d
    [:line {:x1 470 :y1 160 :x2 470 :y2 180
            :stroke "black" :stroke-width 1.2
            :marker-start "url(#arrow)"
            :marker-end "url(#arrow)"}]
    [:text {:x 455 :y 174
            :text-anchor "end"
            :font-size 18
            :font-family "Old Standard"}
     "d"]

    ;; Координата z
    [:line {:x1 320 :y1 205 :x2 660 :y2 205
            :stroke "black" :stroke-width 1.5
            :marker-end "url(#arrow)"}]
    [:text {:x 670 :y 210
            :font-size 18
            :font-family "Old Standard"}
     "z"]

    ;; Левая ампула
    [:path {:d "M110,170
                 C110,120 145,95 190,95
                 L235,95
                 C260,95 278,110 286,130
                 L286,210
                 C278,230 260,245 235,245
                 L190,245
                 C145,245 110,220 110,170 Z"
            :fill "none"
            :stroke "black"
            :stroke-width 2.2}]
    [:line {:x1 286 :y1 170 :x2 290 :y2 170
            :stroke "black" :stroke-width 4}]

    [:text {:x 198 :y 78
            :text-anchor "middle"
            :font-size 18
            :font-family "Old Standard"}
     "исходная ампула"]
    [:text {:x 198 :y 273
            :text-anchor "middle"
            :font-size 18
            :font-family "Old Standard"}
     "Tₛ"]
    [:text {:x 198 :y 296
            :text-anchor "middle"
            :font-size 17
            :font-family "Old Standard"}
     "жидкая Sb"]

    ;; Подача водорода
    [:line {:x1 60 :y1 120 :x2 110 :y2 120
            :stroke "black" :stroke-width 1.6
            :marker-end "url(#arrow)"}]
    [:text {:x 58 :y 108
            :text-anchor "start"
            :font-size 18
            :font-family "Old Standard"}
     "H₂"]

    ;; Левое граничное условие
    [:text {:x 300 :y 80
            :font-size 18
            :font-family "Old Standard"}
     "z = 0"]
    [:text {:x 300 :y 100
            :font-size 17
            :font-family "Old Standard"}
     "p" (subscript "1")
     "(0) = p" (subscript "1,sat")
     "(T" (subscript "s") ")"]

    ;; Правая ампула
    [:path {:d "M694,130
                 C702,110 720,95 745,95
                 L790,95
                 C835,95 870,120 870,170
                 C870,220 835,245 790,245
                 L745,245
                 C720,245 702,230 694,210 Z"
            :fill "none"
            :stroke "black"
            :stroke-width 2.2}]
    [:line {:x1 690 :y1 170 :x2 694 :y2 170
            :stroke "black" :stroke-width 4}]

    [:text {:x 782 :y 78
            :text-anchor "middle"
            :font-size 18
            :font-family "Old Standard"}
     "приёмная ампула"]
    [:text {:x 782 :y 273
            :text-anchor "middle"
            :font-size 18
            :font-family "Old Standard"}
     "Tᵣ"]
    [:text {:x 782 :y 296
            :text-anchor "middle"
            :font-size 17
            :font-family "Old Standard"}
     "осаждение Sb"]

    ;; Правое граничное условие
    [:text {:x 700 :y 55
            :font-size 18
            :font-family "Old Standard"}
     "z = L"]
    [:text {:x 700 :y 100
            :font-size 17
            :font-family "Old Standard"}
     "p" (subscript "1") "(L) = 0"]
    [:text {:x 700 :y 120
            :font-size 17
            :font-family "Old Standard"}
     "p" (subscript "2")
     "(L) = p" (subscript "g")
     " + Q" (subscript "2")
     " / C" (subscript "trap")]

    ;; Ловушка
    [:line {:x1 870 :y1 170 :x2 900 :y2 170
            :stroke "black" :stroke-width 3}]
    [:path {:d "M900,170
                 L900,205
                 Q900,225 920,225
                 Q940,225 940,205
                 L940,135
                 Q940,115 960,115
                 Q980,115 980,135
                 L980,210"
            :fill "none"
            :stroke "black"
            :stroke-width 2.2}]
    [:text {:x 938 :y 100
            :text-anchor "middle"
            :font-size 18
            :font-family "Old Standard"}
     "ловушка"]

    ;; Откачка
    [:line {:x1 980 :y1 210 :x2 980 :y2 255
            :stroke "black" :stroke-width 1.6
            :marker-end "url(#arrow)"}]
    [:text {:x 945 :y 278
            :font-size 18
            :font-family "Old Standard"}
     "к насосу"]]])
