# sb-distillation

Точка входа будущего документа — `notebooks/sb_distillation/report.clj` (`sb-distillation.report`).

## Справочные расчёты

Функции загружаются в REPL из `sb-distillation.reference.legacy`. Исходные ноутбуки доступны в коммите `e08774d`.

| Функции | Происхождение и допущения |
| --- | --- |
| `amp-volume`, `filled-volume`, `liquid-surface-area`, `liquid-height-from-mass` | `sb_distillation.clj`. Горизонтальный цилиндр с двумя полусферическими торцами и плоской поверхностью. Высота отсчитывается от дна и лежит в пределах 0–2 радиуса. Длины в м, масса в кг, плотность в кг/м³. Обратный расчёт высоты использует 60 шагов бисекции. |
| `sb-liquid-vapor-pressure`, `sb-solid-vapor-pressure`, `pa-to-torr` | `sb_distillation.clj`. На входе °C, давление насыщения возвращается в Па. Старый ноутбук указывает диапазоны 785–1100 °C для жидкости и 420–550 °C для твёрдой фазы. Расчёт за их пределами является экстраполяцией. [Данные NIST для жидкости](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440360&Mask=4), [источник для твёрдой фазы, указанный в ноутбуке](https://doi.org/10.1063/1.1701217). |
| `sb4-mean-free-path`, `knudsen-number`, `molecular-conductance`, `viscous-conductance`, `knudsen-correction`, `sb4-tube-mass-flow`, `sb4-distillation-days` | `sb_distillation.clj` и `start.clj`. Модель твёрдых сфер. Число Кнудсена использует диаметр трубки, поправка позднего расчёта — радиус. Трубка прямая, круглая, изотермическая, конденсация в ней отсутствует. Для чистого Sb₄ вход задан давлением насыщения, выход — нулевым давлением. Проводимость в м³/с, расход в кг/с. Оценка времени предполагает постоянный расход. |
| `empirical-knudsen-correction`, `empirical-sb-transfer-rate` | `start.clj`. Ранняя аппроксимация со своими коэффициентами и пересчётом Торр. Давление на выходе равно давлению насыщения твёрдой сурьмы при температуре приёмника. Расход в г/ч. |
| `sb4-hard-sphere-viscosity`, `h2-viscosity`, `mixture-viscosity`, `knudsen-diffusivity`, `binary-diffusivity` | `sb_distillation.clj` и `sb_distillation_trap_study.clj`. Твёрдые сферы для Sb₄ и взаимной диффузии, формула Сазерленда для H₂, правило Уилки для смеси. Вязкость в Па·с, диффузия в м²/с. |
| `pure-dgm-mass-flow`, `mixed-dgm-mass-flow`, `source-partial-pressures`, `rk4-step`, `hydrogen-molar-flow`, `receiver-h2-torr` | `sb_distillation_trap_study.clj`. Изотермическая трубка длиной 0,5 м. Через ловушку проходит только H₂ с расходом 200 мл/ч при 293,15 К и 101325 Па. Давление H₂ за ловушкой 0,03 Торр. Эффективная проводимость ловушки в л/с относится к опорной температуре. Sb₄ на входе насыщен, на выходе его парциальное давление равно нулю. Для смеси используются 120 шагов Рунге — Кутты, до 20 расширений интервала и 20 шагов бисекции, проверка относительной невязки граничного условия с допуском 0,001. |

Графические заготовки лежат в `sb-distillation.reference.plot`. Функция `grid-layer` и карта `theme` перенесены из `plot.clj`. Функция `contour-points` интерполирует уровень по логарифму времени на табличной сетке из исследования ловушки.

## Проба типографики

```powershell
clojure -M:test
clojure -M:clay -r notebooks/typography.clj
```

Откройте `target/clay/typography.html`. Clay передаёт Quarto параметры из `clay.edn` и копирует каталог `styles/` рядом с HTML. Файл `styles/monograph.css` задаёт книжный набор, экранные поля и печать A4. Quarto преобразует формулы TeX в MathML.

STIX Two Text и STIX Two Math подключаются локальными WOFF2 из выпуска 2.13 b171. Лицензия SIL OFL 1.1 приложена в `styles/fonts/OFL.txt`. Для сравнения в образце используется системная Georgia. Ссылка на коммерческую «Обыкновенную новую» приведена в каталоге правообладателя. [Выпуск STIX Two](https://github.com/stipub/stixfonts/releases/tag/v2.13b171), [New Standard в каталоге Paratype](https://www.paratype.ru/catalog/font/pt/new-standard).
