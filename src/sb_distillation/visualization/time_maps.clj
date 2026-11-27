(ns sb-distillation.visualization.time-maps
  (:require [clojure.java.io :as io]
            [scicloj.kindly.v4.kind :as kind])
  (:import [java.util Base64]))

(defn- adjacent-pairs [points]
  (partition 2 1 points))

(defn- regular-rows [points]
  (let [rows (->> points
                  (group-by :temp-c)
                  (sort-by key)
                  (mapv (fn [[_ row]] (vec (sort-by :diameter-mm row)))))
        diameters (mapv :diameter-mm (first rows))
        decreasing? #(every? (fn [[a b]] (> (:days a) (:days b))) (adjacent-pairs %))]
    ;; A level then has one connected branch. The first crossing in each
    ;; horizontal row and at each vertical boundary is unique.
    (when-not (and (seq diameters)
                   (apply < diameters)
                   (every? #(= diameters (mapv :diameter-mm %)) rows)
                   (every? decreasing? rows)
                   (every? decreasing? (apply map vector rows)))
      (throw (ex-info "Expected a rectangular grid with time decreasing along both axes."
                      (select-keys (first points) [:model :h2-flow]))))
    rows))

(defn- crossing [points coordinate level]
  (some (fn [[a b]]
          (when (>= (:days a) level (:days b))
            (let [la (Math/log (:days a))
                  lb (Math/log (:days b))
                  fraction (/ (- (Math/log level) la) (- lb la))]
              (assoc {:diameter-mm (double (:diameter-mm a)) :temp-c (double (:temp-c a))}
                     coordinate (+ (coordinate a) (* fraction (- (coordinate b) (coordinate a))))))))
        (adjacent-pairs points)))

(defn contour-lines
  "Find visible level paths in a monotone time sheet using log-day interpolation.
  Coordinates are mm and °C. Each returned path is ordered by temperature."
  [points levels]
  (let [rows (regular-rows points)
        left (mapv first rows)
        right (mapv peek rows)]
    (into {}
          (keep (fn [level]
                  (let [path (->> (concat (map #(crossing % :diameter-mm level) rows)
                                          [(crossing left :temp-c level)
                                           (crossing right :temp-c level)])
                                  (remove nil?)
                                  (sort-by (juxt :temp-c :diameter-mm))
                                  distinct
                                  vec)]
                    (when (> (count path) 1)
                      [level path]))))
          levels)))

(defn- level-text [days]
  (if (< days 1.0)
    (str (long (* 24.0 days)) " ч")
    (str (long days) " сут")))

(defn- labels [lines [xmin xmax] [ymin ymax]]
  (reduce (fn [chosen [level path]]
            (let [index #(int (* % (dec (count path))))
                  primary? (== 1.0 level)
                  free? (fn [{:keys [diameter-mm temp-c]}]
                          (and (<= (+ xmin 0.34) diameter-mm (- xmax 0.34))
                               (<= (+ ymin 1.7) temp-c (- ymax 1.7))
                               (every? (fn [label]
                                         (or (> (abs (- diameter-mm (:diameter-mm label))) 0.95)
                                             (> (abs (- temp-c (:temp-c label))) 3.2)))
                                       chosen)))
                  point (some #(let [point (path (index %))
                                     candidate (if primary?
                                                 (-> point
                                                     (update :diameter-mm (fn [d] (min d (- xmax 0.38))))
                                                     (update :temp-c (fn [t] (min t (- ymax 1.8)))))
                                                 point)]
                                 (when (free? candidate) candidate))
                              [0.55 0.31 0.77])]
              (cond-> chosen
                point (conj (assoc point :label (level-text level))))))
          []
          (sort-by (fn [[level _]] [(not (== 1.0 level)) level]) lines)))

(defn- data [values]
  {:data {:values values}})

(defn- encodings [x-domain y-domain x-field y-field]
  {:x {:field x-field :type "quantitative" :scale {:domain x-domain :nice false :zero false}
       :axis {:title "Внутренний диаметр, мм" :format ".1~f" :tickMinStep 1}}
   :y {:field y-field :type "quantitative" :scale {:domain y-domain :nice false :zero false}
       :axis {:title "Температура источника, °C" :format "d" :tickMinStep 5}}})

(defn- chart [points domain x-domain y-domain x-step y-step]
  (let [levels (concat [0.25 0.5 0.75] (range 1 (inc (long (second domain)))))
        lines (contour-lines points levels)
        line-points (fn [primary?]
                      (for [[level path] (sort-by key lines)
                            :when (= primary? (== 1.0 level))
                            [order point] (map-indexed vector path)]
                        (assoc point :level level :order order)))
        coords (encodings x-domain y-domain "diameter-mm" "temp-c")
        label-data (labels lines x-domain y-domain)
        [xmin xmax] x-domain
        [ymin ymax] y-domain
        cells (mapv (fn [{:keys [model diameter-mm temp-c mass-flow seconds] :as point}]
                      (assoc point :model (case model :cpim "CPIM" :dgm "DGM")
                             :grams-hour (* 3.6e6 mass-flow) :hours (/ seconds 3600.0)
                             :x0 (max xmin (- diameter-mm (/ x-step 2.0)))
                             :x1 (min xmax (+ diameter-mm (/ x-step 2.0)))
                             :y0 (max ymin (- temp-c (/ y-step 2.0)))
                             :y1 (min ymax (+ temp-c (/ y-step 2.0)))))
                    points)
        line-layer (fn [primary?]
                     (assoc (data (line-points primary?))
                            :mark {:type "line" :clip true :color "#253836"
                                   :strokeWidth (if primary? 2.25 0.8)
                                   :opacity (if primary? 1.0 0.72)}
                            :encoding (assoc coords :detail {:field "level" :type "nominal"}
                                             :order {:field "order" :type "quantitative"})))]
    {:width "container" :height 320
     :usermeta {:embedOptions {:actions false}}
     :autosize {:type "fit" :contains "padding"}
     :config {:font "Old Standard, Georgia, serif"
              :axis {:grid false :labelFontSize 12 :titleFontSize 13 :titlePadding 9}
              :legend {:labelFontSize 11 :titleFontSize 12 :gradientLength 142}
              :view {:stroke nil}}
     :layer [(assoc (data cells)
                    :mark {:type "rect" :clip true}
                    :encoding (assoc (encodings x-domain y-domain "x0" "y0")
                                     :x2 {:field "x1"} :y2 {:field "y1"}
                                     :color {:field "days" :type "quantitative"
                                             :scale {:type "log" :domain domain
                                                     :scheme "redyellowgreen" :reverse true}
                                             :legend {:title "Время, сут" :format ".2~g"
                                                      :values (filterv #(<= (first domain) % (second domain))
                                                                       [0.25 0.5 1 2 4 8 16 32 64])
                                                      :orient "bottom" :gradientLength 185}}
                                     :tooltip [{:field "model" :title "Модель" :type "nominal"}
                                               {:field "h2-flow" :title "H₂, мл/ч" :format "d"}
                                               {:field "temp-c" :title "Температура, °C" :format "d"}
                                               {:field "diameter-mm" :title "Диаметр, мм" :format ".1f"}
                                               {:field "grams-hour" :title "Расход Sb, г/ч" :format ".3f"}
                                               {:field "hours" :title "Время, ч" :format ".2f"}
                                               {:field "days" :title "Время, сут" :format ".3f"}]))
             (line-layer false)
             (line-layer true)
             (assoc (data label-data)
                    :mark {:type "text" :fontSize 11 :fontWeight "bold"
                           :fill "#fffdf3" :stroke "#fffdf3" :strokeWidth 2.5}
                    :encoding (assoc coords :text {:field "label"}))
             (assoc (data label-data)
                    :mark {:type "text" :fontSize 11 :fontWeight "bold"
                           :fill "#17332e"}
                    :encoding (assoc coords :text {:field "label"}))]}))

(defn- font-face [file weight]
  (let [bytes (with-open [stream (io/input-stream (io/resource (str "sb_distillation/visualization/fonts/" file)))]
                (.readAllBytes stream))]
    (str "@font-face{font-family:'Old Standard';font-weight:" weight
         ";font-style:normal;font-display:swap;src:url(data:font/otf;base64,"
         (.encodeToString (Base64/getEncoder) bytes) ") format('opentype');}")))

(defn maps [points {:keys [temperature-c diameter-tenths-mm feeds-ml-h]}]
  (let [[t0 t1 dt] temperature-c
        [d0 d1 dd] diameter-tenths-mm
        days (map :days points)
        domain [(apply min days) (apply max days)]
        grouped (group-by (juxt :h2-flow :model) points)]
    (kind/hiccup
     (into [:section {:class "distillation-maps"}
            [:style
             (str (font-face "OldStandard-Regular.otf" 400)
                  (font-face "OldStandard-Bold.otf" 700)
                  ".distillation-maps{font-family:'Old Standard',Georgia,serif;color:#18332e}"
                  ".time-map-row{margin:1.6rem 0 2rem}"
                  ".time-map-row h3{margin-bottom:.6rem}"
                  ".time-map-pair{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,360px),1fr));gap:1rem}"
                  ".time-map{min-width:0;overflow:hidden;padding:.65rem .7rem .35rem;border:1px solid #d7e0d9;border-radius:6px;background:#fff}"
                  ".time-map h4{font-size:1.13rem;margin:0 .25rem .2rem}"
                  ".time-map .vega-embed{display:block;width:100%;max-width:100%}")]
            ;; A chart can first measure its card while the parser has only
            ;; inserted one of the two grid children.
            [:script "window.addEventListener('load', () => requestAnimationFrame(() => window.dispatchEvent(new Event('resize'))));"]]
           (for [feed feeds-ml-h]
             [:section {:class "time-map-row" :data-feed feed}
              [:h3 (str "Подача H₂ — " feed " мл/ч")]
              (into [:div {:class "time-map-pair"}]
                    (for [[model title] [[:cpim "CPIM"] [:dgm "DGM"]]]
                      [:div {:class "time-map" :data-model (name model)}
                       [:h4 title]
                       (kind/vega-lite
                        (chart (grouped [feed model]) domain
                               [(/ d0 10.0) (/ d1 10.0)] [t0 t1] (/ dd 10.0) dt))]))])))))
