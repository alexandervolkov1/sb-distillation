(require '[clojure.edn :as edn]
         '[clojure.java.io :as io]
         '[clojure.string :as str]
         '[sb-distillation.cpim.model :as model]
         '[sb-distillation.cpim.solver :as solver]
         '[sb-distillation.cpim.transport :as transport]
         '[sb-distillation.solver :as shared-solver])

(def gauge-torr {100 0.02010, 200 0.03680, 300 0.05110})
(def stations (mapv #(/ % 20.0) (range 11)))

(defn problem [temp diameter feed outlet]
  (solver/prepare-problem
   transport/default-parameters
   {:temp temp :diameter diameter :length 0.5
    :ref-temp 293.15 :ref-pressure 101325.0
    :ref-volumetric-flow (/ (* feed 1.0e-6) 3600.0)
    :gauge-pressure outlet :trap-pressure-drop 0.0}))

(defn operating-problem [[temp diameter feed]]
  (problem (+ temp 273.15) (* diameter 1.0e-3) feed
           (* (/ 101325.0 760.0) (gauge-torr feed))))

(defn sampled-pressures [{:keys [coefficients n2]} {:keys [profile n1]} z]
  (if-let [node (some #(when (== (:z %) z) %) profile)]
    (:state node)
    (let [[{za :z pa :state} {zb :z pb :state}]
          (some (fn [[a b]] (when (> (:z a) z (:z b)) [a b]))
                (partition 2 1 profile))
          h (- zb za)
          t (/ (- z za) h)
          s (- 1.0 t)
          [p1a p2a] pa
          [p1b p2b] pb
          da (model/pressure-derivatives coefficients p1a p2a n1 n2)
          db (model/pressure-derivatives coefficients p1b p2b n1 n2)]
      ;; Cubic Hermite interpolation at common stations, using signed dz.
      (mapv (fn [a b a-slope b-slope]
              (+ (* s s (+ 1.0 (* 2.0 t)) a)
                 (* t t (+ 1.0 (* 2.0 s)) b)
                 (* h t s (- (* s a-slope) (* t b-slope)))))
            pa pb da db))))

(defn profile-errors [problem solution expected]
  (apply mapv max
         (map (fn [z expected-state]
                (mapv #(abs (- %1 %2)) (sampled-pressures problem solution z)
                      expected-state))
              stations expected)))

(defn measure [problem]
  (let [evaluate shared-solver/evaluate-trial
        trials (volatile! [])
        start (System/nanoTime)
        result (with-redefs [shared-solver/evaluate-trial
                            (fn [pressure-derivatives problem flux opts]
                              (let [trial (evaluate pressure-derivatives problem flux opts)]
                                (vswap! trials conj [(:attempts trial) (count (:profile trial))])
                                trial))]
                 (solver/solve problem))
        ms (/ (- (System/nanoTime) start) 1.0e6)
        ;; Previous lower endpoints are replaced during bracketing. The two
        ;; final endpoints and all subsequent midpoint profiles are cached.
        cached (when (= :success (:status result))
                 (drop (- (:bracket-trials result) 2) @trials))]
    {:result result
     :ms ms
     :evaluations (count @trials)
     :attempts (reduce + (map first @trials))
     :peak-cached-states (reduce + (map second cached))}))

(defn verify-solution [problem result]
  (let [{:keys [profile mass-flow residual z state]} result]
    (assert (and (= 0.0 z) (= state (:state (peek profile)))))
    (assert (= {:z (:length problem) :state (:outlet-state problem)} (first profile)))
    (assert (apply > (map :z profile)))
    (doseq [{:keys [state]} profile
            :let [[p1 p2] state
                  pressure (+ p1 p2)]]
      (assert (and (Double/isFinite pressure) (>= p1 0.0) (> p2 0.0)))
      (let [x1 (/ p1 pressure) x2 (/ p2 pressure)]
        (assert (and (<= 0.0 x1 1.0) (<= 0.0 x2 1.0)
                     (<= (abs (- (+ x1 x2) 1.0)) 1e-15)))))
    (assert (and (<= (abs residual) (:residual-tol solver/default-solver-options))
                 (= residual (- (first state) (:p1-in problem)))))
    (assert (and (Double/isFinite mass-flow) (pos? mass-flow)))))

(defn grid-case [[temp diameter feed ref-flux ref-p2 ref-mass ref-profile]]
  (let [point [temp diameter feed]
        problem (operating-problem point)
        {:keys [result] :as metrics} (measure problem)]
    (if (= :success (:status result))
      (do
        (verify-solution problem result)
        (merge (dissoc metrics :result)
               {:point point :success true
                :flow (:mass-flow result)
                :hours (/ 0.5 (:mass-flow result) 3600)
                :final-attempts (:attempts result)
                :flux-error (abs (- (/ (:n1 result) ref-flux) 1.0))
                :mass-error (abs (- (/ (:mass-flow result) ref-mass) 1.0))
                :p2-error (abs (- (second (:state result)) ref-p2))
                :residual (abs (:residual result))
                :profile-error (profile-errors problem result ref-profile)}))
      (let [{:keys [trial lower-trial]} result
            comparison (solver/evaluate-trial problem ref-flux)
            above (solver/evaluate-trial problem (* 1.01 ref-flux))]
        (assert (not (contains? result :mass-flow)))
        {:point point :success false :stage (:stage result) :reason (:reason result)
         :trial (dissoc trial :profile) :lower (select-keys lower-trial [:n1 :residual])
         :reference (select-keys comparison [:n1 :status :residual])
         :above-reference (select-keys above [:n1 :status :residual])}))))

(defn extrema [rows key]
  (mapv #(select-keys (apply % key rows) [:point key]) [min-key max-key]))

(defn performance [points]
  (let [rows (mapv (fn [point]
                     (let [{:keys [result] :as metrics} (measure point)]
                       (when (= :success (:status result))
                         (verify-solution point result))
                       (assoc (dissoc metrics :result)
                              :status [(:status result) (:stage result) (:reason result)])))
                   points)
        column (fn [key] (mapv key rows))]
    (println "Performance, measured serial solve calls:" (count rows))
    (println "Result frequencies:" (frequencies (column :status)))
    (doseq [key [:ms :evaluations :attempts :peak-cached-states]
            :let [values (column key)]]
      (println key {:min (apply min values) :max (apply max values)
                    :mean (/ (reduce + values) (double (count rows)))}))
    (println "Serial sum (s):" (/ (reduce + (column :ms)) 1000.0))))

(def variants
  [["default"     1e-8  1e-6  1e-3]
   ["abs-10"      1e-9  1e-6  1e-3]
   ["abs-100"     1e-10 1e-6  1e-3]
   ["rel-10"      1e-8  1e-7  1e-3]
   ["rel-100"     1e-8  1e-8  1e-3]
   ["joint-10"    1e-9  1e-7  1e-3]
   ["joint-100"   1e-10 1e-8  1e-3]
   ["tight-ode"   1e-12 1e-10 1e-3]
   ["root-10"     1e-12 1e-10 1e-4]
   ["root-100"    1e-12 1e-10 1e-5]
   ["root-1000"   1e-12 1e-10 1e-6]
   ["comparison"  1e-12 1e-10 1e-9]
   ["map-proposal" 1e-9 1e-7  1e-4]])

(defn convergence [grid path]
  (with-open [writer (io/writer path)]
    (binding [*out* writer]
      (println "temp,diameter,feed,variant,mode,kind,n1,z,p1,p2,residual")
      (doseq [[temp diameter feed ref-flux] grid
              :let [point [temp diameter feed]]
              :when (#{[680 6 300] [700 10 200] [720 12 100]} point)
              :let [problem (operating-problem point)]
              [variant abs-tol rel-tol residual-tol] variants
              :let [ode-options {:abs-tol abs-tol :rel-tol rel-tol}]
              [mode result] [["fixed" (solver/evaluate-trial problem ref-flux ode-options)]
                             ["shot" (solver/solve problem {:ode-options ode-options
                                                           :residual-tol residual-tol})]]]
        (when (= :success (:status result))
          (doseq [[kind samples] [["node" (:profile result)]
                                  ["station" (for [z stations]
                                               {:z z :state (sampled-pressures problem result z)})]]
                  {:keys [z state]} samples]
            (println (str/join "," (concat point [variant mode kind (:n1 result) z]
                                            state [(:residual result)])))))))))

(let [[ref-path pressure-path full-grid?] *command-line-args*
      {:keys [grid limits]} (edn/read-string (slurp ref-path))
      points (mapv #(operating-problem (take 3 %)) grid)]
  (dotimes [_ 5]
    (doseq [point points] (solver/solve point)))
  (let [rows (mapv grid-case grid)
        good (filter :success rows)
        bad (remove :success rows)]
    (println "Grid successes / failures:" (count good) (count bad))
    (println "Failure frequencies:" (frequencies (map (juxt :stage :reason) bad)))
    (doseq [row (take 3 bad)] (prn row))
    (doseq [key [:flow :hours :final-attempts :attempts :evaluations
                 :flux-error :mass-error :p2-error :residual]]
      (println key (extrema good key)))
    (println "Max profile station errors (Pa):" (apply mapv max (map :profile-error good)))
    (doseq [[diameter feed outlet flux limit-flux pressures limit-pressures] limits
            :let [problem (problem 973.15 diameter feed outlet)
                  result (solver/solve problem {:residual-tol 1.0e-8
                                               :ode-options {:abs-tol 1e-11 :rel-tol 1e-9}})]]
      (assert (= :success (:status result)) (pr-str [diameter outlet result]))
      (println "Limit" [diameter outlet] "relative flux errors [CPIM ref, limit]"
               [(abs (- (/ (:n1 result) flux) 1.0)) (abs (- (/ (:n1 result) limit-flux) 1.0))]
               "pressure errors [CPIM ref, limit]"
               [(profile-errors problem result pressures) (profile-errors problem result limit-pressures)]))
    (convergence grid pressure-path)
    (performance points)
    (when (= "--full-grid" full-grid?)
      (performance
       (for [temp (range 680 721), diameter (range 60 121), feed [100 200 300]]
         (operating-problem [temp (/ diameter 10.0) feed]))))))
