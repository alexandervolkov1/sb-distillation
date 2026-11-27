(ns sb-distillation.solver-test
  (:require [clojure.test :refer [deftest is]]
            [sb-distillation.cpim.solver :as cpim]
            [sb-distillation.cpim.solver-test :refer [conditions]]
            [sb-distillation.dgm.model :as dgm]
            [sb-distillation.physics.transport :as transport]
            [sb-distillation.solver :as solver]))

(defn- baseline []
  (solver/prepare-problem transport/default-parameters conditions))

(defn- exception-data [f]
  (try (f) nil (catch clojure.lang.ExceptionInfo ex (ex-data ex))))

(deftest model-selection-test
  (let [problem (baseline)
        cpim-result (cpim/solve problem)
        default-result (solver/solve problem)]
    (is (= :cpim (:model default-result)))
    (is (= default-result (solver/solve problem {:model :cpim})))
    (is (= cpim-result (dissoc default-result :model)))
    (is (= cpim-result (cpim/solve problem {:model :cpim})))
    (doseq [model [:unknown nil]]
      (is (= {:model model :supported-models [:cpim :dgm]}
             (exception-data #(solver/solve problem {:model model})))))
    (is (= {:model :dgm :supported-models [:cpim]}
           (exception-data #(cpim/solve problem {:model :dgm}))))))

(deftest dgm-baseline-test
  (let [problem (baseline)
        result (solver/solve problem {:model :dgm})
        {:keys [status model n1 residual z state profile mass-flow]} result
        dgm-problem (update problem :coefficients dgm/prepare-coefficients)
        trial (solver/evaluate-trial dgm/pressure-derivatives dgm-problem n1)]
    (is (= [:success :dgm 0.0] [status model z]))
    (is (= {:z (:length problem) :state (:outlet-state problem)} (first profile)))
    (is (= {:z z :state state} (peek profile)))
    (is (apply > (map :z profile)))
    (is (every? (comp solver/valid-pressure-state? :state) profile))
    (is (= (- (first state) (:p1-in problem)) residual))
    (is (<= (abs residual) (:residual-tol solver/default-solver-options)))
    (is (and (Double/isFinite mass-flow) (pos? mass-flow)))
    (is (= trial (select-keys result (keys trial))))))

(deftest model-identifies-failed-solutions-test
  (let [problem (baseline)
        options {:ode-options {:max-attempts 1}}]
    (doseq [model [:cpim :dgm]]
      (let [result (solver/solve problem (assoc options :model model))]
        (is (= [:failed :bracketing :integration-failed model]
               ((juxt :status :stage :reason :model) result)))))
    (is (= (cpim/solve problem options)
           (dissoc (solver/solve problem options) :model)))))
