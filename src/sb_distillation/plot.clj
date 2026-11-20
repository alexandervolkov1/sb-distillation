
(ns sb-distillation.plot)

(def theme
  {:width 800
   :height 600

   :config
   {:axis {:grid false
           :labelFontSize 15
           :titleFontSize 17}

    :legend {:orient "top-left"
             :offset 8
             :fillColor "white"
             :strokeColor "lightgray"
             :padding 6
             :labelFontSize 15
             :symbolSize 180}

    :title {:orient "bottom"
            :anchor "middle"
            :fontSize 18
            :offset 18}}})

(defn grid-layer
  [axis values scale major?]
  (let [color (if major? "#bbbbbb" "#dddddd")
        width (if major? 1.2 0.6)]
    {:data {:values
            (mapv (fn [value] {:position value}) values)}
     :mark {:type "rule"
            :color color
            :strokeWidth width}
     :encoding
     {axis {:field "position"
            :type "quantitative"
            :scale scale}}}))
