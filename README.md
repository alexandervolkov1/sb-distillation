# sb-distillation

Численное моделирование переноса паров сурьмы в водороде на Clojure. Математическая модель описана в [отдельном проекте Typst](https://github.com/alexandervolkov1/sb-distillation-theory).

## Структура

- `src/sb_distillation/reference/legacy.clj` — справочные расчёты давления насыщенного пара, проводимости, диффузии, вязкости и переноса газовой смеси. [Единицы и допущения](https://github.com/alexandervolkov1/sb-distillation/blob/4784a763ff7d0c385fdeb27f44311637af7b2950/README.md#справочные-расчёты)
- `src/sb_distillation/reference/plot.clj` — независимые функции подготовки графиков и изолиний.
- `notebooks/sb_distillation/analysis.clj` — вычислительный ноутбук Clay.
- `test/sb_distillation_test.clj` — тесты справочных расчётов.

## Запуск

Из корня проекта выполните команды для запуска тестов и сборки ноутбука.

```powershell
clojure -M:test
clojure -M:clay -r notebooks/sb_distillation/analysis.clj
```

Clay создаёт HTML в `target/clay/sb_distillation.analysis.html`.
