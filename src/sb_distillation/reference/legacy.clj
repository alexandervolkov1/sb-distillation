(ns sb-distillation.reference.legacy
  "Справочные расчёты из трёх исходных ноутбуков. Единицы и допущения описаны в README.")

(def gas-constant 8.314462618)
(def boltzmann-constant 1.380649e-23)
(def avogadro-constant 6.02214076e23)
(def pa-per-torr 133.32236842105263)
(def h2-molar-mass 0.00201588)
(def sb4-molar-mass 0.48704)
(def h2-collision-diameter-m 0.289e-9)
(def sb4-collision-diameter-m 0.6e-9)

(defn amp-volume
  "Внутренний объём ампулы, м³. Радиус и длина цилиндрической части в метрах."
  [radius straight-length]
  (+ (* Math/PI radius radius straight-length)
     (* 4/3 Math/PI radius radius radius)))

(defn filled-volume
  "Объём расплава, м³. Размеры в метрах, высота от дна в пределах [0, 2r]."
  [radius straight-length liquid-height]
  (let [cylinder-volume
        (* straight-length
           (- (* radius radius (Math/acos (/ (- radius liquid-height)
                                             radius)))
              (* (- radius liquid-height)
                 (Math/sqrt (* liquid-height
                               (- (* 2.0 radius) liquid-height))))))
        end-caps-volume
        (* Math/PI liquid-height liquid-height
           (- radius (/ liquid-height 3.0)))]
    (+ cylinder-volume end-caps-volume)))

(defn liquid-surface-area
  "Площадь плоской свободной поверхности, м². Размеры в метрах, высота в [0, 2r]."
  [radius straight-length liquid-height]
  (let [a-squared (* liquid-height (- (* 2.0 radius) liquid-height))
        cylinder-area (* 2.0 straight-length (Math/sqrt a-squared))
        end-caps-area (* Math/PI a-squared)]
    (+ cylinder-area end-caps-area)))

(defn liquid-height-from-mass
  "Высота расплава, м. Масса в кг, плотность в кг/м³."
  [radius straight-length mass density]
  (when-not (pos? density)
    (throw (ex-info "Density must be positive."
                    {:density density})))
  (let [target-volume (/ mass density)
        total-volume (amp-volume radius straight-length)]
    (when-not (<= 0.0 target-volume total-volume)
      (throw (ex-info "Liquid volume is outside ampoule capacity."
                      {:volume target-volume
                       :capacity total-volume})))
    (cond
      (zero? target-volume) 0.0

      (= target-volume total-volume)
      (* 2.0 radius)

      :else
      (loop [low 0.0
             high (* 2.0 radius)
             steps 0]
        (let [mid (/ (+ low high) 2.0)]
          (if (= steps 60)
            mid
            (if (< (filled-volume radius straight-length mid)
                   target-volume)
              (recur mid high (inc steps))
              (recur low mid (inc steps)))))))))

(defn- liquid-vapor-pressure-bar
  [temp-c]
  (let [temp-k (+ temp-c 273.15)]
    (Math/pow 10.0 (- 2.26041 (/ 4475.449 (- temp-k 152.352))))))

(defn sb-liquid-vapor-pressure
  "Давление насыщенного пара над жидкой сурьмой, Па. Температура в °C."
  [temp-c]
  (* 1.0e5 (liquid-vapor-pressure-bar temp-c)))

(defn- solid-vapor-pressure-torr
  [temp-c]
  (let [temp-k (+ temp-c 273.15)]
    (Math/pow 10.0 (- 10.571 (/ 10300.0 temp-k)))))

(defn sb-solid-vapor-pressure
  "Давление насыщенного пара над твёрдой сурьмой, Па. Температура в °C."
  [temp-c]
  (* pa-per-torr (solid-vapor-pressure-torr temp-c)))

(defn pa-to-torr
  [pressure-pa]
  (/ pressure-pa pa-per-torr))

(defn sb4-mean-free-path
  "Длина свободного пробега в модели твёрдых сфер, м."
  [temp-c pressure-pa collision-diameter-m]
  (let [temp-k (+ temp-c 273.15)]
    (when-not (and (pos? temp-k)
                   (pos? pressure-pa)
                   (pos? collision-diameter-m))
      (throw (ex-info "Invalid gas parameters."
                      {:temp-c temp-c
                       :pressure-pa pressure-pa
                       :collision-diameter-m collision-diameter-m})))
    (/ (* boltzmann-constant temp-k)
       (* (Math/sqrt 2.0)
          Math/PI
          collision-diameter-m
          collision-diameter-m
          pressure-pa))))

(defn knudsen-number
  [mean-free-path-m tube-diameter-m]
  (when-not (pos? tube-diameter-m)
    (throw (ex-info "Tube diameter must be positive."
                    {:tube-diameter-m tube-diameter-m})))
  (/ mean-free-path-m tube-diameter-m))

(defn- hard-sphere-viscosity
  [temp-k collision-diameter-m]
  (let [molecule-mass (/ sb4-molar-mass avogadro-constant)]
    (/ (* 5.0 (Math/sqrt (/ (* molecule-mass boltzmann-constant temp-k)
                             Math/PI)))
       (* 16.0 collision-diameter-m collision-diameter-m))))

(defn sb4-hard-sphere-viscosity
  "Динамическая вязкость Sb₄, Па·с. Температура в °C, диаметр молекулы в метрах."
  [temp-c collision-diameter-m]
  (let [temp-k (+ temp-c 273.15)]
    (when-not (and (pos? temp-k)
                   (pos? collision-diameter-m))
      (throw (ex-info "Invalid gas parameters."
                      {:temp-c temp-c
                       :collision-diameter-m collision-diameter-m})))
    (hard-sphere-viscosity temp-k collision-diameter-m)))

(defn molecular-conductance
  "Молекулярная проводимость длинной круглой трубки, м³/с. Аргументы в СИ."
  [temp-k tube-diameter-m tube-length-m molar-mass-kg-per-mol]
  (* (/ (* Math/PI (Math/pow tube-diameter-m 3.0)) (* 12.0 tube-length-m))
     (Math/sqrt (/ (* 8.0 gas-constant temp-k)
                   (* Math/PI molar-mass-kg-per-mol)))))

(defn viscous-conductance
  "Вязкая проводимость круглой трубки при среднем давлении, м³/с."
  [tube-diameter-m tube-length-m mean-pressure-pa viscosity-pa-s]
  (/ (* Math/PI (Math/pow tube-diameter-m 4.0) mean-pressure-pa)
     (* 128.0 viscosity-pa-s tube-length-m)))

(defn knudsen-correction
  "Поправка к молекулярной проводимости из sb_distillation.clj."
  [tube-radius-m mean-free-path-m]
  (let [ratio (/ tube-radius-m mean-free-path-m)]
    (/ (+ 1.0 (* 2.507 ratio))
       (+ 1.0 (* 3.095 ratio)))))

(defn empirical-knudsen-correction
  "Поправка из start.clj. Обе проводимости в одинаковых единицах."
  [molecular viscous]
  (let [x (* (/ viscous molecular) (/ 12.1 135.0))]
    (/ (+ 1.0 (* 192.0 x))
       (+ 1.0 (* 237.0 x)))))

(defn empirical-sb-transfer-rate
  "Расход Sb₄ по start.clj, г/ч. Размеры в м, температуры в °C."
  [{:keys [tube-diameter tube-length tube-temp source-temp receiver-temp]}]
  (let [temp-k (+ tube-temp 273.15)
        p1-torr (* (liquid-vapor-pressure-bar source-temp) 750.061683)
        p2-torr (solid-vapor-pressure-torr receiver-temp)
        p-mean (/ (+ (* 133.322368 p1-torr) (* 133.322368 p2-torr)) 2.0)
        pressure-drop-pa (* 133.322368 (- p1-torr p2-torr))
        molecular (molecular-conductance temp-k tube-diameter tube-length sb4-molar-mass)
        viscosity (hard-sphere-viscosity temp-k sb4-collision-diameter-m)
        viscous (viscous-conductance tube-diameter tube-length p-mean viscosity)
        conductance (+ viscous (* molecular (empirical-knudsen-correction molecular viscous)))
        mol-per-sec (/ (* conductance pressure-drop-pa) (* gas-constant temp-k))]
    (* mol-per-sec sb4-molar-mass 3.6e6)))

(defn sb4-tube-mass-flow
  "Расход чистого Sb₄, кг/с, с поправкой Кнудсена и нулевым давлением на выходе."
  [temp-c tube-diameter-m tube-length-m collision-diameter-m]
  (when-not (and (pos? tube-diameter-m)
                 (pos? tube-length-m))
    (throw (ex-info "Invalid tube geometry."
                    {:tube-diameter-m tube-diameter-m
                     :tube-length-m tube-length-m})))
  (let [temp-k (+ temp-c 273.15)
        p1 (sb-liquid-vapor-pressure temp-c)
        p-mean (/ p1 2.0)
        viscosity (sb4-hard-sphere-viscosity temp-c collision-diameter-m)
        mean-free-path (sb4-mean-free-path temp-c p-mean collision-diameter-m)
        z (knudsen-correction (/ tube-diameter-m 2.0) mean-free-path)
        molecular (molecular-conductance temp-k tube-diameter-m tube-length-m sb4-molar-mass)
        viscous (viscous-conductance tube-diameter-m tube-length-m p-mean viscosity)
        conductance (+ viscous (* z molecular))]
    (* (/ sb4-molar-mass (* gas-constant temp-k)) conductance p1)))

(defn sb4-distillation-days
  "Оценка длительности в сутках при постоянном массовом расходе."
  [mass-kg temp-c tube-diameter-m tube-length-m collision-diameter-m]
  (when-not (pos? mass-kg)
    (throw (ex-info "Mass must be positive."
                    {:mass-kg mass-kg})))
  (/ mass-kg
     (* 86400.0
        (sb4-tube-mass-flow
         temp-c
         tube-diameter-m
         tube-length-m
         collision-diameter-m))))

;; Параметры расчётов из sb_distillation_trap_study.clj.
(def tube-length-m 0.5)
(def h2-flow-ml-h 200.0)
(def gauge-h2-torr 0.03)
(def reference-temp-k 293.15)
(def reference-pressure-pa 101325.0)

(defn- sb4-viscosity [temp-k]
  (hard-sphere-viscosity temp-k sb4-collision-diameter-m))

(defn h2-viscosity
  "Динамическая вязкость H₂ по формуле Сазерленда, Па·с."
  [temp-k]
  (* 8.76e-6
     (Math/pow (/ temp-k 293.85) 1.5)
     (/ (+ 293.85 72.0) (+ temp-k 72.0))))

(defn- wilke-phi [eta-i eta-j molar-mass-i molar-mass-j]
  (/ (Math/pow (+ 1.0
                  (* (Math/sqrt (/ eta-i eta-j))
                     (Math/pow (/ molar-mass-j molar-mass-i) 0.25)))
               2.0)
     (Math/sqrt (* 8.0 (+ 1.0 (/ molar-mass-i molar-mass-j))))))

(defn- wilke-viscosity [eta-h eta-s phi-hs phi-sh y-h y-s]
  (+ (/ (* y-h eta-h) (+ y-h (* y-s phi-hs)))
     (/ (* y-s eta-s) (+ y-s (* y-h phi-sh)))))

(defn mixture-viscosity
  "Вязкость смеси H₂–Sb₄ по Уилки, Па·с. Доля H₂ молярная."
  [temp-k h2-fraction]
  (let [sb4-fraction (- 1.0 h2-fraction)
        eta-h (h2-viscosity temp-k)
        eta-s (sb4-viscosity temp-k)
        phi-hs (wilke-phi eta-h eta-s h2-molar-mass sb4-molar-mass)
        phi-sh (wilke-phi eta-s eta-h sb4-molar-mass h2-molar-mass)]
    (wilke-viscosity eta-h eta-s phi-hs phi-sh h2-fraction sb4-fraction)))

(defn knudsen-diffusivity
  "Коэффициент диффузии Кнудсена, м²/с. Молярная масса в кг/моль."
  [temp-k diameter-m molar-mass]
  (* (/ diameter-m 3.0)
     (Math/sqrt (/ (* 8.0 gas-constant temp-k)
                   (* Math/PI molar-mass)))))

(defn binary-diffusivity
  "Коэффициент взаимной диффузии H₂ и Sb₄ в модели твёрдых сфер, м²/с."
  [temp-k pressure-pa]
  (let [m-h (/ h2-molar-mass avogadro-constant)
        m-s (/ sb4-molar-mass avogadro-constant)
        reduced-mass (/ (* m-h m-s) (+ m-h m-s))
        collision-diameter-m (/ (+ h2-collision-diameter-m
                                   sb4-collision-diameter-m) 2.0)]
    (/ (* 3.0 boltzmann-constant temp-k
          (Math/sqrt (/ (* 2.0 boltzmann-constant temp-k)
                        (* Math/PI reduced-mass))))
       (* 16.0 pressure-pa collision-diameter-m collision-diameter-m))))

(defn pure-dgm-mass-flow
  "Расход чистого Sb₄ по однокомпонентной DGM, кг/с."
  [temp-c diameter-m]
  (let [temp-k (+ temp-c 273.15)
        p (sb-liquid-vapor-pressure temp-c)
        area (/ (* Math/PI diameter-m diameter-m) 4.0)
        dk (knudsen-diffusivity temp-k diameter-m sb4-molar-mass)
        b0 (/ (* diameter-m diameter-m) 32.0)
        viscosity (sb4-viscosity temp-k)
        flux (/ (+ (* dk p) (/ (* b0 p p) (* 2.0 viscosity)))
                (* gas-constant temp-k tube-length-m))]
    (* flux area sb4-molar-mass)))

(defn hydrogen-molar-flow
  "Расход H₂, моль/с. Объёмный расход в мл/ч при опорных температуре и давлении."
  [flow-ml-h]
  (/ (* reference-pressure-pa flow-ml-h 1.0e-6)
     (* 3600.0 gas-constant reference-temp-k)))

(defn receiver-h2-torr
  "Давление H₂ перед ловушкой, Торр. Проводимость в л/с, nil означает нулевой перепад."
  [trap-conductance-l-s]
  (let [q-torr-l-s (/ (* (/ reference-pressure-pa pa-per-torr)
                         h2-flow-ml-h)
                      3.6e6)]
    (+ gauge-h2-torr
       (if (nil? trap-conductance-l-s)
         0.0
         (/ q-torr-l-s trap-conductance-l-s)))))

(defn- dgm-properties [temp-k diameter-m]
  (let [eta-h (h2-viscosity temp-k)
        eta-s (sb4-viscosity temp-k)]
    {:eta-h eta-h
     :eta-s eta-s
     :phi-hs (wilke-phi eta-h eta-s h2-molar-mass sb4-molar-mass)
     :phi-sh (wilke-phi eta-s eta-h sb4-molar-mass h2-molar-mass)
     :d-h (knudsen-diffusivity temp-k diameter-m h2-molar-mass)
     :d-s (knudsen-diffusivity temp-k diameter-m sb4-molar-mass)
     ;; В модели твёрдых сфер произведение давления на Dₕₛ постоянно.
     :d-hs-p (binary-diffusivity temp-k 1.0)
     :b0 (/ (* diameter-m diameter-m) 32.0)
     :rt (* gas-constant temp-k)}))

(defn- dgm-gradient [properties h2-flux sb4-flux [p-h p-s]]
  (when-not (and (pos? p-h) (<= 0.0 p-s))
    (throw (ex-info "Nonphysical partial pressure." {:p-h p-h :p-s p-s})))
  (let [{:keys [eta-h eta-s phi-hs phi-sh
                d-h d-s d-hs-p b0 rt]} properties
        p (+ p-h p-s)
        y-h (/ p-h p)
        y-s (/ p-s p)
        eta (wilke-viscosity eta-h eta-s phi-hs phi-sh y-h y-s)
        interaction (/ (- (* p-s h2-flux) (* p-h sb4-flux)) d-hs-p)
        f-h (* rt (+ (/ h2-flux d-h) interaction))
        f-s (* rt (- (/ sb4-flux d-s) interaction))
        a-h (/ (* b0 p-h) (* eta d-h))
        a-s (/ (* b0 p-s) (* eta d-s))
        total-gradient (/ (+ f-h f-s) (+ 1.0 a-h a-s))]
    [(- f-h (* a-h total-gradient))
     (- f-s (* a-s total-gradient))]))

(defn- v-add [a b factor]
  [(+ (nth a 0) (* factor (nth b 0)))
   (+ (nth a 1) (* factor (nth b 1)))])

(defn rk4-step
  "Один шаг Рунге — Кутты четвёртого порядка для двух компонент."
  [gradient state ds]
  (let [k1 (gradient state)
        k2 (gradient (v-add state k1 (* ds 0.5)))
        k3 (gradient (v-add state k2 (* ds 0.5)))
        k4 (gradient (v-add state k3 ds))]
    [(+ (nth state 0)
        (* (/ ds 6.0) (+ (nth k1 0) (* 2.0 (nth k2 0))
                         (* 2.0 (nth k3 0)) (nth k4 0))))
     (+ (nth state 1)
        (* (/ ds 6.0) (+ (nth k1 1) (* 2.0 (nth k2 1))
                         (* 2.0 (nth k3 1)) (nth k4 1))))]))

(defn source-partial-pressures
  "Парциальные давления [H₂ Sb₄] на входе, Па. Расход Sb₄ в моль/с."
  [temp-c diameter-m receiver-torr sb4-molar-flow]
  (let [temp-k (+ temp-c 273.15)
        area (/ (* Math/PI diameter-m diameter-m) 4.0)
        nh (/ (hydrogen-molar-flow h2-flow-ml-h) area)
        ns (/ sb4-molar-flow area)
        steps 120
        ds (/ tube-length-m steps)
        properties (dgm-properties temp-k diameter-m)
        gradient #(dgm-gradient properties nh ns %)]
    (loop [i 0
           state [(* receiver-torr pa-per-torr) 0.0]]
      (if (= i steps)
        state
        (recur (inc i) (rk4-step gradient state ds))))))

(defn mixed-dgm-mass-flow
  "Расход Sb₄ в смеси с H₂ по бинарной DGM, кг/с."
  [temp-c diameter-m receiver-torr]
  (let [target (sb-liquid-vapor-pressure temp-c)
        initial (/ (pure-dgm-mass-flow temp-c diameter-m)
                   sb4-molar-mass)
        ;; Нефизичная проба задаёт верхнюю границу, конечное решение проверяется заново.
        error
        (fn [sb4-molar-flow]
          (try
            (- (nth (source-partial-pressures
                     temp-c diameter-m receiver-torr sb4-molar-flow) 1)
               target)
            (catch clojure.lang.ExceptionInfo _ ##Inf)))
        high
        (loop [h initial
               i 0]
          (if (or (>= (error h) 0.0) (= i 20))
            h
            (recur (* h 1.5) (inc i))))]
    (when-not (>= (error high) 0.0)
      (throw (ex-info "Could not bracket Sb4 flow."
                      {:temp-c temp-c :diameter-m diameter-m
                       :receiver-torr receiver-torr})))
    (let [flow
          (loop [lo 0.0 hi high i 0]
            (if (= i 20)
              (* 0.5 (+ lo hi))
              (let [mid (* 0.5 (+ lo hi))]
                (if (neg? (error mid))
                  (recur mid hi (inc i))
                  (recur lo mid (inc i))))))
          p-source (source-partial-pressures
                    temp-c diameter-m receiver-torr flow)]
      (when (> (/ (Math/abs (- (nth p-source 1) target)) target) 1.0e-3)
        (throw (ex-info "Boundary condition not satisfied."
                        {:temp-c temp-c :diameter-m diameter-m
                         :receiver-torr receiver-torr
                         :p-source p-source :target target})))
      (* flow sb4-molar-mass))))
