(require '[sb-distillation.cpim.solver :as cpim]
         '[sb-distillation.cpim.transport :as transport])

(def gauge-torr {100 0.02010, 200 0.03680, 300 0.05110})

(defn operating-grid []
  (vec
   (for [temp (range 680 721 5), diameter (range 6 13), feed [100 200 300]]
     (cpim/prepare-problem
      transport/default-parameters
      {:temp (+ temp 273.15) :diameter (* diameter 1.0e-3) :length 0.5
       :ref-temp 293.15 :ref-pressure 101325.0
       :ref-volumetric-flow (/ (* feed 1.0e-6) 3600.0)
       :gauge-pressure (* (/ 101325.0 760.0) (gauge-torr feed))
       :trap-pressure-drop 0.0}))))

(defn checked-time [solve problem]
  (let [start (System/nanoTime)
        result (solve problem)
        elapsed (- (System/nanoTime) start)]
    (assert (= :success (:status result)) (pr-str result))
    elapsed))

(let [[mode warmups passes] *command-line-args*
      solve (if (= "before" mode)
              #(cpim/solve % {})
              (let [options {:model (case mode "cpim" :cpim "dgm" :dgm)}
                    solve-model (requiring-resolve 'sb-distillation.solver/solve)]
                #(solve-model % options)))
      warmups (parse-long warmups)
      passes (parse-long passes)
      problems (operating-grid)]
  (prn {:mode mode :cases (count problems) :warmup-passes warmups :measured-passes passes
        :jvm (System/getProperty "java.runtime.version")
        :jvm-args (vec (.getInputArguments (java.lang.management.ManagementFactory/getRuntimeMXBean)))})
  (dotimes [_ warmups]
    (doseq [problem problems]
      (assert (= :success (:status (solve problem))))))
  (dotimes [pass passes]
    (let [times (mapv #(checked-time solve %) problems)
          total (reduce + times)]
      (prn {:pass (inc pass) :mean-ms (/ total (* 1.0e6 (count times)))
            :max-ms (/ (apply max times) 1.0e6)
            :solve-seconds (/ total 1.0e9)}))))
