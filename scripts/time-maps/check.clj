(require '[clojure.string :as str]
         '[sb-distillation.experiments.time-grid :as grid]
         '[sb-distillation.visualization.time-maps :as viz])

(def study grid/default-study)
(def result (grid/cached-grid study "target/distillation-maps"))
(def groups (group-by (juxt :h2-flow :model) (:rows result)))
(defn degrees-and-mm [row] ((juxt :temp-c :diameter-mm) row))

(assert (= 8 (count groups)))
(assert (every? #(= 2501 (count %) (count (distinct (map degrees-and-mm %)))) (vals groups)))
(assert (= (map #(dissoc % :model) (groups [0 :cpim]))
           (map #(dissoc % :model) (groups [0 :dgm]))))
(assert (not= (grid/cache-key study) (grid/cache-key (assoc study :length 0.6))))
(assert (not= (grid/cache-key study)
              (grid/cache-key (assoc-in study [:solver-options :residual-tol] 1e-4))))

(doseq [model grid/models
        :let [name (name model)
              archived (rest (str/split-lines (slurp (str "target/cpim-verification/" name "-dense.csv"))))]]
  ;; Walk the map rows grouped by model, using the same integer-index key as
  ;; the archived study. The two scalar values should match exactly.
  (let [current (into {} (for [feed [100 200 300]
                               {:keys [temp-c diameter-mm mass-flow seconds]} (groups [feed model])]
                           [[temp-c diameter-mm feed] [mass-flow seconds]]))]
    (assert (= 7503 (count archived) (count current)))
    (doseq [line archived
            :let [[t d feed _ _ flow time] (str/split line #",")
                  key (mapv parse-double [t d feed])
                  key [(long (key 0)) (key 1) (long (key 2))]]]
      (assert (= [(parse-double flow) (parse-double time)] (current key)) (pr-str [model key])))
    (println "Archived rows match exactly" model (count archived))))

(defn diameter-at [path temperature]
  (some (fn [[a b]]
          (when (<= (:temp-c a) temperature (:temp-c b))
            (let [fraction (/ (- temperature (:temp-c a)) (- (:temp-c b) (:temp-c a)))]
              (+ (:diameter-mm a) (* fraction (- (:diameter-mm b) (:diameter-mm a)))))))
        (partition 2 1 path)))

(defn refined-diameter [model feed temperature level coarse]
  (let [centre (Math/round (* coarse 100.0))
        points (mapv #(grid/calculate-point study model feed temperature (/ % 100.0))
                     (range (max 600 (- centre 7)) (min 1201 (+ centre 8))))]
    (some (fn [[a b]]
            (when (>= (:days a) level (:days b))
              (let [la (Math/log (:days a))
                    lb (Math/log (:days b))]
                (+ (:diameter-mm a)
                   (* 0.01 (/ (- (Math/log level) la) (- lb la)))))))
          (partition 2 1 points))))

(let [max-days (apply max (map :days (:rows result)))
      levels (concat [0.25 0.5 0.75] (range 1 (inc (long max-days))))
      samples
      (vec (for [[[feed model] points] (sort-by key groups)
                 :let [paths (viz/contour-lines points levels)]
                 [temp preferred fixed-level]
                 (cond-> [[685.5 6.6 nil] [700.5 9.0 nil] [715.5 11.2 nil]]
                   (paths 1) (conj [(/ (+ (:temp-c (first (paths 1)))
                                          (:temp-c (peek (paths 1)))) 2.0) 12.0 1]))
                 :let [candidates (for [[level path] paths
                                        :let [d (diameter-at path temp)]
                                        :when (and d (if fixed-level (= level fixed-level)
                                                        (<= 6.15 d 11.85)))]
                                    [level d])
                       [level coarse] (apply min-key #(abs (- preferred (second %))) candidates)
                       fine (refined-diameter model feed temp level coarse)
                       error (abs (- coarse fine))]]
             {:feed feed :model model :temperature temp :level level
              :coarse-mm coarse :refined-mm fine :error-mm error}))]
  (assert (= 28 (count samples)))
  (doseq [sample samples]
    (assert (<= (:error-mm sample) 0.015) (pr-str sample))
    (prn sample))
  (prn {:maximum-contour-error-mm (apply max (map :error-mm samples))
        :calculation-seconds (:calculation-seconds result)})
  (doseq [[[feed model] points] (sort-by key groups)]
    (prn {:feed feed :model model :minimum (apply min-key :days points)
          :maximum (apply max-key :days points)})))
