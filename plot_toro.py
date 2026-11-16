#!/usr/bin/env python3
"""Графики численного расчёта из solution.csv и run.txt (solver-cli)."""
import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
import numpy as np

BLUE, GRAY = "#2563a6", "#8593a3"
TITLES = {
    1: "Задача Сода", 2: "Две волны разрежения",
    3: "Сильный разрыв давления слева", 4: "Сильный разрыв давления справа",
    5: "Столкновение сильных ударных волн",
}
LABELS = {
    "density": ("Плотность", "ρ"), "velocity": ("Скорость", "u"),
    "pressure": ("Давление", "p"), "momentum": ("Плотность импульса", "ρu"),
    "energy": ("Плотность полной энергии", "E"),
    "internal_energy": ("Удельная внутренняя энергия", "e"),
    "sound_speed": ("Скорость звука", "a"),
    "entropy": ("Энтропийная функция", "p / ρᵞ"),
}
CSV_FIELDS = {"density", "momentum", "energy", "velocity", "pressure"}


@dataclass
class Case:
    directory: Path
    meta: dict
    times: np.ndarray
    x: np.ndarray
    values: dict

    @property
    def name(self):
        return self.directory.name

    @property
    def title(self):
        parts = self.name.rsplit("_", 1)
        if len(parts) == 2 and parts[0] == "toro" and parts[1].isdigit():
            number = int(parts[1])
            if number in TITLES:
                return f"ТОРО {number:02d} · {TITLES[number]}"
        return self.name


def read_metadata(path):
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        key, sep, value = line.partition("=")
        if not sep or not key or key in result:
            raise ValueError(f"{path}: некорректная или повторная запись {line!r}")
        result[key] = value
    return result


def number(meta, key):
    value = float(meta[key])
    if not np.isfinite(value):
        raise ValueError(f"run.txt: {key} должно быть конечным числом")
    return value


def positive_int(meta, key):
    value = int(meta[key])
    if value <= 0:
        raise ValueError(f"run.txt: {key} должно быть положительным целым")
    return value


def read_case(directory):
    """Проверяем слои и сортируем по (time, cell), не полагаясь на порядок строк."""
    meta = read_metadata(directory / "run.txt")
    if meta.get("status") != "complete":
        raise ValueError(f"{directory}: требуется завершённый расчёт (status=complete)")
    if meta.get("model") != "euler_1d":
        raise ValueError(f"{directory}: поддерживается только model=euler_1d")
    cells = positive_int(meta, "cells")
    snapshots = positive_int(meta, "snapshots")
    gamma = number(meta, "gamma")
    xmin, xmax = number(meta, "x_min"), number(meta, "x_max")
    start, end = number(meta, "t_start"), number(meta, "t_end")
    if gamma <= 1 or xmax <= xmin or end < start or number(meta, "final_time") != end:
        raise ValueError(f"{directory}: неверные gamma, границы или конечное время")

    path = directory / "solution.csv"
    with path.open(encoding="utf-8", newline="") as stream:
        header = next(csv.reader(stream), [])
    if len(header) != len(set(header)) or not {"time", "cell", "x"}.issubset(header):
        raise ValueError(f"{path}: нужны уникальные столбцы time,cell,x и поля решения")
    if set(header) - {"time", "cell", "x"} - CSV_FIELDS:
        raise ValueError(f"{path}: неизвестные поля решения")
    if not CSV_FIELDS.intersection(header):
        raise ValueError(f"{path}: нет полей решения")
    rows = np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True, dtype=float,
                                     encoding="utf-8"))
    if len(rows) != cells * snapshots or any(not np.all(np.isfinite(rows[k])) for k in header):
        raise ValueError(f"{path}: неполные слои или нечисловые значения")
    rows = rows[np.lexsort((rows["cell"], rows["time"]))]
    cell_ids = rows["cell"].reshape(snapshots, cells)
    times = rows["time"].reshape(snapshots, cells)
    coordinates = rows["x"].reshape(snapshots, cells)
    if not np.all(cell_ids == np.arange(cells)) or not np.all(times == times[:, :1]):
        raise ValueError(f"{path}: пропущенные/повторные ячейки или смешанные времена слоя")
    times = times[:, 0]
    if times[0] != start or times[-1] != end or np.any(np.diff(times) <= 0):
        raise ValueError(f"{path}: неверные времена слоёв; нужны начальный и конечный")
    dx = (xmax - xmin) / cells
    expected_x = xmin + (np.arange(cells) + 0.5) * dx
    tolerance = max(dx * 1e-8, np.finfo(float).eps * max(abs(xmin), abs(xmax)) * 8)
    if not np.allclose(coordinates, expected_x[None, :], rtol=0, atol=tolerance):
        raise ValueError(f"{path}: координаты не соответствуют равномерной сетке из run.txt")
    values = {k: rows[k].reshape(snapshots, cells) for k in header if k in CSV_FIELDS}
    for key in ("density", "pressure"):
        if key in values and np.any(values[key] <= 0):
            raise ValueError(f"{path}: {key} должно быть положительным")

    def derived(key, value):
        if not np.all(np.isfinite(value)):
            raise ValueError(f"{path}: невозможно вычислить конечное поле {key}")
        if key in values and not np.allclose(values[key], value, rtol=1e-8, atol=1e-12):
            raise ValueError(f"{path}: {key} не согласуется с консервативными величинами")
        values.setdefault(key, value)

    # Восстанавливаем недостающие примитивные поля только при наличии нужных U.
    with np.errstate(over="raise", divide="raise", invalid="raise"):
        if {"density", "momentum"}.issubset(values):
            derived("velocity", values["momentum"] / values["density"])
        if {"density", "velocity", "energy"}.issubset(values):
            derived("pressure", (gamma - 1) * (values["energy"] -
                    0.5 * values["density"] * values["velocity"] ** 2))
        if "pressure" in values and np.any(values["pressure"] <= 0):
            raise ValueError(f"{path}: вычислено неположительное давление")
        if {"density", "pressure"}.issubset(values):
            ratio = values["pressure"] / values["density"]
            derived("internal_energy", ratio / (gamma - 1))
            derived("sound_speed", np.sqrt(gamma * ratio))
            derived("entropy", values["pressure"] / values["density"] ** gamma)
    return Case(directory, meta, times, coordinates[0], values)


def heading(fig, case, subtitle):
    fig.suptitle(case.title, fontsize=16, fontweight="bold", color="#223142")
    fig.text(0.5, 0.925, f"{subtitle} | γ = {number(case.meta, 'gamma'):g} | "
             f"t = {case.times[-1]:g} | {len(case.x)} ячеек", ha="center", fontsize=10)


def profiles(case):
    preferred = ["density", "velocity", "pressure", "internal_energy", "sound_speed", "entropy"]
    fields = [key for key in preferred if key in case.values]
    if not fields:
        fields = [key for key in ("momentum", "energy") if key in case.values]
    nrows = (len(fields) + 1) // 2
    fig, axes = plt.subplots(nrows, 2, figsize=(12, 3 * nrows + 1.2), squeeze=False)
    heading(fig, case, "Численное решение · Годунов, первый порядок")
    for ax, key in zip(axes.flat, fields):
        label, symbol = LABELS[key]
        if len(case.times) > 1:
            ax.plot(case.x, case.values[key][0], color=GRAY, ls="--", lw=1,
                    label=f"Начальное, t={case.times[0]:g}")
        ax.plot(case.x, case.values[key][-1], color=BLUE, lw=1.7,
                label=f"Конечное, t={case.times[-1]:g}")
        ax.set(title=label, xlabel="x", ylabel=symbol,
               xlim=(number(case.meta, "x_min"), number(case.meta, "x_max")))
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4), useOffset=False)
        ax.margins(y=0.1)
    for ax in list(axes.flat)[len(fields):]:
        ax.set_visible(False)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
    fig.subplots_adjust(left=0.08, right=0.97, bottom=0.09, top=0.86, hspace=0.55, wspace=0.3)
    return fig


def time_edges(times):
    """Середины между сохранёнными моментами; крайние границы не выходят за расчёт."""
    return np.concatenate(([times[0]], (times[:-1] + times[1:]) / 2, [times[-1]]))


def xt_diagram(case):
    fields = [k for k in ("density", "velocity", "pressure") if k in case.values]
    if not fields:
        fields = [next(iter(case.values))]
    fig, axes = plt.subplots(1, len(fields), figsize=(5 * len(fields), 6.4), squeeze=False)
    heading(fig, case, "Поля x–t · сохранённые слои")
    xedges = np.linspace(number(case.meta, "x_min"), number(case.meta, "x_max"), len(case.x) + 1)
    for ax, key in zip(axes.flat, fields):
        label, symbol = LABELS[key]
        if len(case.times) == 1:
            ax.text(0.5, 0.5, "Сохранён один момент времени.\nДиаграмма x–t не строится.",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            continue
        values = case.values[key]
        # При адаптивном dt моменты неравномерны: используем реальные границы времени.
        logarithmic = key != "velocity" and values.min() > 0 and values.max() / values.min() > 20
        norm = LogNorm(vmin=values.min(), vmax=values.max()) if logarithmic else None
        ax.grid(False)
        mesh = ax.pcolormesh(xedges, time_edges(case.times), values, shading="flat",
                             cmap="RdBu_r" if key == "velocity" else "YlGnBu",
                             norm=norm, rasterized=True)
        fig.colorbar(mesh, ax=ax, pad=0.03, label=symbol + (" (лог. шкала)" if logarithmic else ""))
        ax.set(title=label, xlabel="x", ylabel="t", ylim=(case.times[0], case.times[-1]))
    fig.subplots_adjust(left=0.075, right=0.86 if len(fields) == 1 else 0.92,
                        bottom=0.12, top=0.84, wspace=0.55)
    return fig


def pu_diagram(case):
    """Точки ячеек в p–u; это НЕ кривые локального решателя Римана."""
    fig, ax = plt.subplots(figsize=(10, 7))
    heading(fig, case, "Диаграмма p–u · состояния ячеек на конечном слое")
    scatter = ax.scatter(case.values["velocity"][-1], case.values["pressure"][-1],
                         c=case.x, cmap="viridis", s=22, alpha=0.8, edgecolors="none",
                         vmin=number(case.meta, "x_min"), vmax=number(case.meta, "x_max"))
    fig.colorbar(scatter, ax=ax, label="Координата ячейки x")
    pressure = case.values["pressure"][-1]
    if pressure.max() / pressure.min() > 100:
        ax.set_yscale("log")
    ax.set(xlabel="Скорость u", ylabel="Давление p")
    fig.text(0.5, 0.035, "Цвет обозначает координату x. Точные волновые кривые здесь не вычисляются.",
             ha="center", fontsize=10)
    fig.subplots_adjust(left=0.1, right=0.96, bottom=0.15, top=0.85)
    return fig


def discover_cases(path):
    if (path / "solution.csv").is_file():
        return [path]
    directories = sorted((p for p in path.glob("toro_*") if p.is_dir()),
                         key=lambda p: (int(p.name[5:]) if p.name[5:].isdigit() else float("inf"), p.name))
    if not directories:
        raise ValueError(f"{path}: нужны каталоги toro_*/solution.csv + run.txt или один каталог расчёта")
    return directories


def render(cases, output):
    output.mkdir(parents=True, exist_ok=True)
    report = ["# Численные расчёты Торо", "",
              "Источник: solution.csv и run.txt программы solver-cli. "
              "Годунов первого порядка; показаны сохранённые слои.", "",
              "Плотность, импульс и энергия — средние по ячейкам. Примитивные и производные "
              "величины восстановлены из этих данных. Точный эталон и ошибки относительно него не вычисляются.", "",
              "На x–t цвет привязан к ближайшему сохранённому моменту; промежуточные шаги не восстанавливаются.", ""]
    summary_rows = []
    figures = 0
    with PdfPages(output / "toro_figures.pdf") as pdf:
        for case in cases:
            directory = output / case.name
            directory.mkdir(parents=True, exist_ok=True)
            drawings = [("profiles", profiles), ("xt", xt_diagram)]
            if {"pressure", "velocity"}.issubset(case.values):
                drawings.append(("pu", pu_diagram))
            for name, draw in drawings:
                fig = draw(case)
                try:
                    fig.savefig(directory / f"{name}.png", dpi=150)
                    pdf.savefig(fig)
                finally:
                    plt.close(fig)
                figures += 1
            report += [f"## {case.title}", "",
                       f"Ячеек: {len(case.x)}. Сохранённых слоёв: {len(case.times)}. "
                       f"Время: {case.times[0]:g} → {case.times[-1]:g}.", "",
                       "| Поле | Минимум в конце | Максимум в конце |", "|---|---:|---:|"]
            for key, values in case.values.items():
                low, high = float(values[-1].min()), float(values[-1].max())
                summary_rows.append([case.name, case.times[-1], len(case.x), len(case.times), key, low, high])
                report.append(f"| {LABELS[key][0]} | {low:.8g} | {high:.8g} |")
            for name, _ in drawings:
                report += ["", f"![{name}]({case.name}/{name}.png)"]
            if len(drawings) == 2:
                report += ["", "Диаграмма p–u пропущена: недостаточно полей для давления и скорости."]
            report += [""]
            print(f"{case.name}: {len(drawings)} рисунка, {len(case.times)} слоёв")
    with (output / "toro_summary.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["case", "final_time", "cells", "snapshots", "field", "min_final", "max_final"])
        writer.writerows(summary_rows)
    (output / "REPORT.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Готово: {figures} PNG, toro_figures.pdf, toro_summary.csv, REPORT.md\n{output.resolve()}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, help="Каталог с toro_* или один каталог с solution.csv/run.txt")
    parser.add_argument("--output-dir", type=Path, help="Каталог рисунков; по умолчанию RESULTS/plots")
    args = parser.parse_args()
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11,
                         "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": "#ced6df", "axes.labelcolor": "#334155", "text.color": "#223142",
                         "axes.grid": True, "grid.color": "#e5eaf0", "grid.linewidth": 0.6,
                         "axes.axisbelow": True, "figure.facecolor": "white", "savefig.facecolor": "white"})
    try:
        # Чтение всех входов завершается до создания выходных файлов.
        cases = [read_case(path) for path in discover_cases(args.results)]
        render(cases, args.output_dir or args.results / "plots")
    except (OSError, ValueError, KeyError, FloatingPointError) as error:
        parser.exit(1, f"Ошибка: {error}\n")


if __name__ == "__main__":
    main()
