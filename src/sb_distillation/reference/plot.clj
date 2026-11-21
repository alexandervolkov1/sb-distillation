(ns sb-distillation.reference.plot)

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
  "Слой линий сетки Vega-Lite для заданной оси и значений."
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

(defn- level-label [days]
  (if (< days 1.0)
    (str (int (* 24.0 days)) " h")
    (str (int days) " d")))

(defn contour-points
  "Изолинии на сетке из исследования ловушки. Каждая строка содержит точки одной
  температуры с ключами :temp-c, :diameter-mm и :days. Уровни и времена положительны.
  Берётся первое пересечение уровня, с интерполяцией по логарифму времени."
  [rows levels]
  (vec
   (for [level levels
         row rows
         :let [pair (first
                     (filter (fn [[a b]]
                               (<= (min (:days a) (:days b))
                                   level
                                   (max (:days a) (:days b))))
                             (partition 2 1 row)))]
         :when pair
         :let [[a b] pair
               log-a (Math/log (:days a))
               log-b (Math/log (:days b))
               ratio (if (= log-a log-b)
                       0.5
                       (/ (- (Math/log level) log-a)
                          (- log-b log-a)))]]
     {:temp-c (:temp-c a)
      :diameter-mm (+ (:diameter-mm a)
                      (* ratio (- (:diameter-mm b) (:diameter-mm a))))
      :days level
      :label (level-label level)})))
