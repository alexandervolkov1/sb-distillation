"""Compare the two verified scalar grids produced by study.clj."""

import argparse
import csv
from decimal import Decimal
from itertools import product
from math import isclose, isfinite
from pathlib import Path
from statistics import mean


TEMPERATURES = range(680, 721)
DIAMETERS = range(60, 121)
FEEDS = (100, 200, 300)


def read_grid(path):
    grid = {}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            tenths = Decimal(row["diameter_mm"]) * 10
            if tenths != tenths.to_integral_value():
                raise ValueError(f"Invalid diameter in {path}: {row}")
            point = int(row["temp_c"]), int(tenths), int(row["feed_ml_h"])
            flow, seconds = (float(row[key]) for key in ("mass_flow_kg_s", "time_s"))
            if (point in grid or row["status"] != "success"
                    or not all(isfinite(value) and value > 0 for value in (flow, seconds))
                    or not isclose(flow * seconds, 0.5, rel_tol=1e-13)
                    or not abs(float(row["residual_pa"])) <= 1e-3):
                raise ValueError(f"Invalid result in {path}: {row}")
            grid[point] = flow, seconds
    if set(grid) != set(product(TEMPERATURES, DIAMETERS, FEEDS)):
        raise ValueError(f"Unexpected grid points in {path}")
    return grid


def location(point):
    return f"{point[0]} °C, {point[1] / 10:.1f} мм"


def summary(values, feed):
    selected = {point: value for point, value in values.items() if point[2] == feed}
    lo = min(selected, key=selected.get)
    hi = max(selected, key=selected.get)
    return (f"{selected[lo]:+.3f} ({location(lo)}) | "
            f"{selected[hi]:+.3f} ({location(hi)}) | {mean(selected.values()):+.3f}")


def trend(values, feed, axis):
    changes = []
    for t, d in product(TEMPERATURES, DIAMETERS):
        next_t, next_d = (t + 1, d) if axis == "temperature" else (t, d + 1)
        if (next_t, next_d, feed) in values:
            changes.append(values[next_t, next_d, feed] - values[t, d, feed])
    return f"{sum(v < 0 for v in changes)} / {sum(v > 0 for v in changes)} / {sum(v == 0 for v in changes)}"


def report(cpim, dgm):
    flow = {point: 100 * (dgm[point][0] / cpim[point][0] - 1) for point in cpim}
    time = {point: 100 * (dgm[point][1] / cpim[point][1] - 1) for point in cpim}
    print("# Сравнение CPIM и DGM в рабочем диапазоне")
    print()
    print("Сетки этапа 27.5 содержат по 7503 решения при 680–720 °C с шагом 1 °C, "
          "диаметре 6–12 мм с шагом 0,1 мм и трёх расходах H₂. "
          "Относительное различие рассчитано как `100 × (DGM / CPIM − 1)`. "
          "Средние взяты по всем 2501 точкам для каждого расхода. Время относится к 500 г сурьмы.")
    for label, values in (("Массовый расход", flow), ("Время перегонки", time)):
        print(f"\n## {label}\n")
        print("| H₂, мл/ч | Минимум, % | Максимум, % | Среднее, % |")
        print("| ---: | --- | --- | ---: |")
        for feed in FEEDS:
            print(f"| {feed} | {summary(values, feed)} |")
    print("\n## Зависимость от диаметра и температуры\n")
    print("Числа в таблице показывают, сколько соседних пар дали уменьшение, увеличение "
          "или равенство относительной разницы массового расхода при увеличении одной координаты.")
    print("\n| H₂, мл/ч | При увеличении диаметра | При увеличении температуры |")
    print("| ---: | ---: | ---: |")
    for feed in FEEDS:
        print(f"| {feed} | {trend(flow, feed, 'diameter')} | {trend(flow, feed, 'temperature')} |")
    print("\nНа этой сетке разница убывает с температурой. По диаметру она убывает "
          "на всей сетке для 100 мл/ч и на большей части сетки для 200 и 300 мл/ч. "
          "Локальные максимумы при 200 и 300 мл/ч смещены внутрь диапазона диаметров. "
          "Разница времени меняется в противоположном направлении.")
    print(f"\nНа сетке разница массовых расходов составляет {min(flow.values()):+.3f}…"
          f"{max(flow.values()):+.3f} %. Предварительная оценка порядка 1–11 % "
          "согласуется с рассчитанным диапазоном.")
    if min(flow.values()) > 0:
        print("Знак различия сохраняется по всей сетке. DGM даёт больший массовый расход "
              "и меньшее время, поскольку время обратно пропорционально расходу.")
    print("\n## Базовый режим и угловые точки\n")
    print("| °C | мм | H₂, мл/ч | CPIM, г/ч | DGM, г/ч | Δ расхода, % | CPIM, ч | DGM, ч | Δ времени, % |")
    print("| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    corners = product((680, 720), (60, 120), FEEDS)
    for point in ((700, 100, 200), *corners):
        t, d, feed = point
        cpim_flow, cpim_time = cpim[point]
        dgm_flow, dgm_time = dgm[point]
        print(f"| {t} | {d / 10:.1f} | {feed} | {cpim_flow * 3.6e6:.4f} | "
              f"{dgm_flow * 3.6e6:.4f} | {flow[point]:+.3f} | "
              f"{cpim_time / 3600:.3f} | {dgm_time / 3600:.3f} | {time[point]:+.3f} |")
    print("\nНа независимых сетках максимальные относительные ошибки массового расхода "
          "составили 0,00289 % для CPIM и 0,00251 % для DGM. Даже их сумма на два порядка "
          "меньше наименьшего различия моделей. Экспериментальных данных для выбора "
          "модели по точности массового расхода пока нет. [Численная верификация](DGM.md)")
    print("\n## Воспроизведение\n")
    print("Входные CSV создаёт `study.clj` в полном прогоне с `--dense-csv`. "
          "Для опубликованного сравнения использованы сетки этапа 27.5. Версию решателя "
          "проверили по исходному коду, базовому режиму и угловым точкам. [Команды расчёта сеток](DGM.md#reference-and-reproduction)")
    print("\n```powershell\npython -X utf8 scripts/cpim-verification/compare.py target/cpim-verification/cpim-dense.csv target/cpim-verification/dgm-dense.csv\n```")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cpim", type=Path)
    parser.add_argument("dgm", type=Path)
    inputs = parser.parse_args()
    report(read_grid(inputs.cpim), read_grid(inputs.dgm))
