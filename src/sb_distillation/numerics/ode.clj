(ns sb-distillation.numerics.ode)

(defn rk4-step
  "y_{n+1} - Advance the state by one classical fourth-order
  Runge-Kutta step.
  f takes the independent coordinate and state vector
  and returns a derivative vector of the same dimension.
  h is a signed integration step."
  [f z y h]
  (let [half-h (* 0.5 h)
        k1 (f z y)
        k2 (f (+ z half-h)
              (mapv #(+ %1 (* half-h %2)) y k1))
        k3 (f (+ z half-h)
              (mapv #(+ %1 (* half-h %2)) y k2))
        k4 (f (+ z h)
              (mapv #(+ %1 (* h %2)) y k3))]
    (mapv (fn [yi a b c d]
            (+ yi (* (/ h 6.0)
                     (+ a (* 2.0 b) (* 2.0 c) d))))
          y k1 k2 k3 k4)))

(defn rk4-step-doubling
  "y_f, e - RK4 step with local error estimation
  using step doubling.
  Returns the fine state and the estimated signed
  correction (y_f - y_c)/15.
  h is a signed integration step."
  [f z y h]
  (let [half-h (* 0.5 h)
        coarse (rk4-step f z y h)
        middle (rk4-step f z y half-h)
        fine (rk4-step f (+ z half-h) middle half-h)
        error (mapv #(/ (- %1 %2) 15.0) fine coarse)]
    {:state fine
     :error error}))
