#!/usr/bin/env python3
"""Графики Торо: прямое точное решение solver-cli и эталонные тесты."""
import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import numpy as np

BLUE, RED, GRAY = "#2563a6", "#b44b3e", "#8593a3"
TITLES = {1: "Задача Сода", 2: "Две волны разрежения", 3: "Сильный разрыв давления слева",
          4: "Сильный разрыв давления справа", 5: "Столкновение ударных волн"}
LABELS = {"density": "Плотность ρ", "velocity": "Скорость u", "pressure": "Давление p",
          "momentum": "Плотность импульса ρu", "energy": "Плотность полной энергии ρE",
          "internal_energy": "Удельная внутренняя энергия e", "sound_speed": "Скорость звука a",
          "entropy": "Энтропийная функция p / ρᵞ"}
QUANTITIES = {"p_star": "p*", "u_star": "u*", "rho_star_left": "ρ*L", "rho_star_right": "ρ*R"}
GUESS = {"linearized": "акустическая", "two-rarefactions": "две ВР", "two-shocks": "две УВ"}
ERROR_COLUMNS = ["quantity", "actual", "reference", "absolute_error", "relative_error", "tolerance", "passed"]


@dataclass
class Case:
    directory: Path
    meta: dict
    times: np.ndarray
    x: np.ndarray
    values: dict
    xt_x: np.ndarray
    xt_rho: np.ndarray = None
    curves: np.ndarray = None
    errors: list = None

    @property
    def exact(self):
        return self.errors is not None

    @property
    def direct_exact(self):
        return self.meta.get("method") == "exact_riemann"

    @property
    def name(self):
        return self.directory.name

    @property
    def title(self):
        suffix = self.meta.get("test", self.name.rsplit("_", 1)[-1])
        number = int(suffix) if suffix.isdigit() else 0
        return f"Торо {number}: {TITLES[number]}" if number in TITLES else self.name


def number(meta, key):
    value = float(meta[key])
    if not np.isfinite(value):
        raise ValueError(f"{key}: требуется конечное число")
    return value


def table(path):
    data = np.atleast_1d(np.genfromtxt(path, names=True, delimiter=",", encoding="utf-8"))
    if any(not np.all(np.isfinite(data[k])) for k in data.dtype.names):
        raise ValueError(f"{path}: нечисловые значения")
    return data


def grid(data, points):
    index = "point" if "point" in data.dtype.names else "cell"
    data = data[np.lexsort((data[index], data["time"]))]
    cells = data[index].reshape(-1, points)
    times = data["time"].reshape(-1, points)
    x = data["x"].reshape(-1, points)
    if not np.all(cells == np.arange(points)) or not np.all(times == times[:, :1]):
        raise ValueError("Неполный слой, повторные ячейки или смешанные времена")
    if not np.all(x == x[:1]) or np.any(np.diff(times[:, 0]) <= 0):
        raise ValueError("Изменилась сетка или нарушен порядок времён")
    values = {k: data[k].reshape(-1, points) for k in data.dtype.names if k in LABELS}
    return times[:, 0], x[0], values


def physical(values):
    for key in ("density", "pressure"):
        if key in values and np.any(values[key] <= 0):
            raise ValueError(f"{key}: требуются положительные значения")


def read_case(directory):
    if (directory / "profile.csv").exists():
        with (directory / "summary.csv").open(encoding="utf-8") as f:
            meta = dict(list(csv.reader(f))[1:])
        with (directory / "errors.csv").open(encoding="utf-8") as f:
            errors = list(csv.DictReader(f))
        if meta["passed"] != "1" or len(errors) != 4 or any(r["passed"] != "1" for r in errors):
            raise ValueError(f"{directory}: точный тест Торо не пройден")
        profile = table(directory / "profile.csv")
        data = table(directory / "density_xt.csv")
        times, xt_x, fields = grid(data, len(np.unique(data["x"])))
        values = {k: profile[k][None, :] for k in profile.dtype.names if k in LABELS}
        physical(values); physical(fields)
        return Case(directory, meta, times, profile["x"], values, xt_x, fields["density"],
                    table(directory / "wave_curves.csv"), errors)

    meta = dict(line.split("=", 1) for line in (directory / "run.txt").read_text().splitlines() if line.strip())
    if meta.get("status") != "complete":
        raise ValueError(f"{directory}: требуется status=complete")
    if meta.get("model") != "euler_1d":
        raise ValueError("Производные поля рассчитаны для модели идеального газа euler_1d")
    gamma = number(meta, "gamma")
    xmin, xmax = number(meta, "x_min"), number(meta, "x_max")
    start, end = number(meta, "t_start"), number(meta, "t_end")
    cells = int(meta["points"] if "points" in meta else meta["cells"])
    snapshots = int(meta["snapshots"])
    if gamma <= 1 or xmax <= xmin or end < start or cells <= 0 or snapshots <= 0 or number(meta, "final_time") != end:
        raise ValueError("Неверные параметры газа, сетки или конечного времени")
    times, x, values = grid(table(directory / "solution.csv"), cells)
    dx = (xmax - xmin) / cells
    atol = max(dx * 1e-8, np.finfo(float).eps * max(abs(xmin), abs(xmax)) * 8)
    if len(times) != snapshots or times[0] != start or times[-1] != end:
        raise ValueError("Отсутствует начальный или конечный слой")
    if not np.allclose(x, xmin + (np.arange(cells) + .5) * dx, rtol=0, atol=atol):
        raise ValueError("Координаты не соответствуют сетке")
    physical(values)

    def derived(key, value):
        if not np.all(np.isfinite(value)) or (key in values and not np.allclose(values[key], value, rtol=1e-8, atol=1e-12)):
            raise ValueError(f"{key}: физические поля не согласованы с консервативными")
        values.setdefault(key, value)

    with np.errstate(over="raise", divide="raise", invalid="raise"):
        if {"density", "momentum"} <= values.keys():
            derived("velocity", values["momentum"] / values["density"])
        if {"density", "velocity", "energy"} <= values.keys():
            derived("pressure", (gamma - 1) * (values["energy"] - .5 * values["density"] * values["velocity"]**2))
        physical(values)
        if {"density", "pressure"} <= values.keys():
            derived("internal_energy", values["pressure"] / values["density"] / (gamma - 1))
            derived("sound_speed", np.sqrt(gamma * values["pressure"] / values["density"]))
            derived("entropy", values["pressure"] / values["density"]**gamma)
    return Case(directory, meta, times, x, values, x, values.get("density"))


def heading(fig, case, subtitle):
    fig.suptitle(case.title, fontsize=15, fontweight="bold")
    fig.text(.5, .925, f"{subtitle} | γ={number(case.meta, 'gamma'):g} | t={case.times[-1]:g}", ha="center")


def waves(case):
    for side in ("left", "right"):
        yield side, case.meta[side + "_type"] == "shock", number(case.meta, side + "_head"), number(case.meta, side + "_tail")


def profiles(case):
    fields = [k for k in LABELS if k in case.values and k not in ("momentum", "energy")]
    fields = fields or list(case.values)
    fig, axes = plt.subplots((len(fields)+1)//2, 2, figsize=(12, 3*((len(fields)+1)//2)+1.2), squeeze=False)
    heading(fig, case, "Точное решение" if case.exact or case.direct_exact else "Годунов первого порядка")
    for ax, key in zip(axes.flat, fields):
        if not case.exact and len(case.times) > 1:
            ax.plot(case.x, case.values[key][0], color=GRAY, ls="--", label="Начальное")
        ax.plot(case.x, case.values[key][-1], color=BLUE, label="Конечное")
        if case.exact or case.direct_exact:
            for _, shock, head, tail in waves(case):
                for speed in ([head] if shock else [head, tail]):
                    ax.axvline(number(case.meta, "x0") + speed*np.ptp(case.times), color=GRAY, lw=.8, ls="-" if shock else ":")
            ax.axvline(number(case.meta, "x0") + number(case.meta, "contact")*np.ptp(case.times), color=RED, ls="--", lw=.9)
        ax.set(title=LABELS[key], xlabel="x", xlim=(number(case.meta, "x_min"), number(case.meta, "x_max")))
        ax.margins(y=.1)
    for ax in list(axes.flat)[len(fields):]:
        ax.set_visible(False)
    if case.exact or case.direct_exact:
        handles = [Line2D([], [], color=GRAY, label="УВ"), Line2D([], [], color=GRAY, ls=":", label="Границы ВР"),
                   Line2D([], [], color=RED, ls="--", label="Контакт")]
    else:
        handles, _ = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles=handles, loc="lower center", ncol=3)
    fig.subplots_adjust(left=.08, right=.97, bottom=.1, top=.85, hspace=.55, wspace=.3)
    return fig


def time_edges(times):
    return np.r_[times[0], (times[:-1]+times[1:])/2, times[-1]]


def xt_diagram(case):
    fig, ax = plt.subplots(figsize=(9, 6.5))
    heading(fig, case, "Диаграмма x–t")
    if len(case.times) == 1:
        ax.text(.5, .5, "Сохранён один момент: эволюция не строится", ha="center", transform=ax.transAxes)
    else:
        key = "density" if case.xt_rho is not None else next(iter(case.values))
        values = case.xt_rho if case.xt_rho is not None else case.values[key]
        xedges = time_edges(case.xt_x) if case.exact else np.linspace(number(case.meta, "x_min"), number(case.meta, "x_max"), len(case.xt_x)+1)
        norm = LogNorm(values.min(), values.max()) if values.min() > 0 and values.max()/values.min() > 20 else None
        ax.grid(False)
        image = ax.pcolormesh(xedges, time_edges(case.times), values, cmap="YlGnBu", norm=norm, rasterized=True)
        fig.colorbar(image, ax=ax, label=LABELS[key])
        if case.exact:
            x0, elapsed = number(case.meta, "x0"), case.times-case.times[0]
            for _, shock, head, tail in waves(case):
                for speed in ([head] if shock else np.linspace(head, tail, 9)):
                    ax.plot(x0+speed*elapsed, case.times, color=".25", lw=1.2 if shock or speed in (head, tail) else .5)
            ax.plot(x0+number(case.meta, "contact")*elapsed, case.times, color=RED, ls="--", lw=1.3)
            for side, star, level in (("left", False, .82), ("left", True, .68), ("right", True, .48), ("right", False, .82)):
                dt = level*np.ptp(case.times)
                boundary = number(case.meta, side + ("_tail" if star else "_head"))
                if side == "left":
                    a, b = (x0+boundary*dt, x0+number(case.meta, "contact")*dt) if star else (number(case.meta, "x_min"), x0+boundary*dt)
                else:
                    a, b = (x0+number(case.meta, "contact")*dt, x0+boundary*dt) if star else (x0+boundary*dt, number(case.meta, "x_max"))
                label = "$W_{" + ("*" if star else "") + side[0].upper() + "}$"
                ax.text((a+b)/2, case.times[0]+dt, label, ha="center", fontsize=9, bbox=dict(fc="white", ec="none", alpha=.8))
            ax.legend(handles=[Line2D([], [], color=".25", label="УВ / веер ВР"),
                               Line2D([], [], color=RED, ls="--", label="Контакт")], loc="lower right", fontsize=8)
    ax.set(xlabel="x", ylabel="t", xlim=(number(case.meta, "x_min"), number(case.meta, "x_max")), ylim=(case.times[0], case.times[-1]) if len(case.times)>1 else None)
    fig.subplots_adjust(left=.1, right=.88, bottom=.13, top=.85)
    return fig


def pu_diagram(case):
    fig, ax = plt.subplots(figsize=(9, 6.5))
    heading(fig, case, "Волновые кривые p–u" if case.exact else
            "Точный профиль p–u" if case.direct_exact else "Состояния ячеек p–u")
    if case.exact:
        p = case.curves["p"]
        for side, color in (("left", BLUE), ("right", RED)):
            pk = number(case.meta, side + "_p")
            for mask, style, kind in ((p <= pk, "--", "ВР"), (p >= pk, "-", "УВ")):
                label = ("Левая" if side == "left" else "Правая") + ": " + kind
                ax.plot(case.curves["u_"+side][mask], p[mask], color=color, ls=style, label=label)
            u = number(case.meta, side+"_u")
            ax.scatter(u, pk, color=color, s=30)
            ax.annotate("$W_L$" if side == "left" else "$W_R$", (u, pk), xytext=(5, 7), textcoords="offset points")
        ax.scatter(number(case.meta, "u_star"), number(case.meta, "p_star"), marker="*", s=120, color="black", label="$(u_*,p_*)$")
        ax.legend(fontsize=9, ncol=2)
        logarithmic = int(case.meta["test"]) >= 3
    elif case.direct_exact:
        image = ax.scatter(case.values["velocity"][-1], case.values["pressure"][-1],
                           c=case.x, s=18, cmap="viridis")
        fig.colorbar(image, ax=ax, label="x")
        logarithmic = np.ptp(np.log10(case.values["pressure"][-1])) > 2
    else:
        image = ax.scatter(case.values["velocity"][-1], case.values["pressure"][-1], c=case.x, s=18, cmap="viridis")
        fig.colorbar(image, ax=ax, label="x")
        logarithmic = np.ptp(np.log10(case.values["pressure"][-1])) > 2
    if logarithmic:
        ax.set_yscale("log")
    else:
        ax.set_ylim(bottom=0)
    ax.set(xlabel="u", ylabel="p")
    fig.subplots_adjust(left=.1, right=.96, bottom=.13, top=.85)
    return fig


def report_section(case):
    lines = [f"## {case.title}", "", f"γ={number(case.meta, 'gamma'):g}; t={case.times[-1]:g}; точек профиля: {len(case.x)}.", ""]
    if case.exact:
        lines += ["| Величина | Расчёт | Торо | Абс. отклонение | Отн. отклонение | Допуск |", "|---|---:|---:|---:|---:|---:|"]
        for row in case.errors:
            values = [f"{float(row[k]):.9g}" if row[k] else "—" for k in ERROR_COLUMNS[1:-1]]
            lines.append("| " + " | ".join([QUANTITIES[row["quantity"]]]+values) + " |")
        lines += ["", f"GUESSP: {GUESS[case.meta['guess']]}; p0={number(case.meta, 'initial_pressure'):.9g}; итераций: {case.meta['iterations']}.", "",
                  "| Волна | Тип | Голова / скачок | Хвост |", "|---|---|---:|---:|"]
        for side, shock, head, tail in waves(case):
            label = "Левая" if side == "left" else "Правая"
            lines.append(f"| {label} | {'УВ' if shock else 'ВР'} | {head:.9g} | {'—' if shock else f'{tail:.9g}'} |")
        lines += [f"| Контакт | контакт | {number(case.meta, 'contact'):.9g} | — |", ""]
        for side, shock, _, _ in waves(case):
            lines.append(f"{'Слева' if side == 'left' else 'Справа'}: p* {'>' if shock else '<'} pK, поэтому {'УВ' if shock else 'ВР'}.")
        same = np.isclose(number(case.meta, "rho_star_left"), number(case.meta, "rho_star_right"), rtol=1e-7)
        lines += ["", "На контакте p и u непрерывны. " + ("Звёздные состояния совпадают, видимого скачка нет." if same else "Скачки видны у ρ, e, a и s."),
                  "В веере s=p/ρᵞ постоянно; при прохождении ударной волны энтропийная функция возрастает.", ""]
    elif case.direct_exact:
        lines += ["Точный профиль одной задачи Римана; значения взяты в точках x без шагов по времени.", "",
                  f"p*={number(case.meta, 'p_star'):.9g}; u*={number(case.meta, 'u_star'):.9g}; "
                  f"ρ*L={number(case.meta, 'rho_star_left'):.9g}; ρ*R={number(case.meta, 'rho_star_right'):.9g}.", "",
                  "| Волна | Тип | Голова / скачок | Хвост |", "|---|---|---:|---:|"]
        for side, shock, head, tail in waves(case):
            lines.append(f"| {'Левая' if side == 'left' else 'Правая'} | "
                         f"{'УВ' if shock else 'ВР'} | {head:.9g} | "
                         f"{'—' if shock else f'{tail:.9g}'} |")
        lines += [f"| Контакт | контакт | {number(case.meta, 'contact'):.9g} | — |", ""]
    else:
        lines += ["Сеточный расчёт Годунова: средние по ячейкам; точные волновые кривые не вычисляются.", "",
                  "| Поле | Минимум | Максимум |", "|---|---:|---:|"]
        lines += [f"| {LABELS[k]} | {v[-1].min():.8g} | {v[-1].max():.8g} |" for k, v in case.values.items()]
    return lines


def discover_cases(path):
    if (path / "profile.csv").exists() or (path / "solution.csv").exists():
        return [path]
    paths = sorted(p for p in path.glob("toro_*") if p.is_dir())
    if not paths:
        raise ValueError(f"{path}: нет результатов toro_*")
    return paths


def results_page(case):
    fig, ax = plt.subplots(figsize=(11.7, 8.3))
    heading(fig, case, "Сравнение с таблицей Торо")
    ax.set_axis_off()
    rows = [[QUANTITIES[r["quantity"]]] + [f"{float(r[k]):.8g}" if r[k] else "—"
            for k in ERROR_COLUMNS[1:-1]] for r in case.errors]
    error_table = ax.table(cellText=rows, colLabels=["Величина", "Расчёт", "Торо", "Абс. откл.", "Отн. откл.", "Допуск"],
                           bbox=(0, .57, 1, .25), cellLoc="center")
    error_table.auto_set_font_size(False)
    error_table.set_fontsize(10)
    rows = [["Левая" if side == "left" else "Правая", "УВ" if shock else "ВР",
             f"{head:.8g}", "—" if shock else f"{tail:.8g}"] for side, shock, head, tail in waves(case)]
    rows.append(["Контакт", "контакт", f"{number(case.meta, 'contact'):.8g}", "—"])
    ax.table(cellText=rows, colLabels=["Волна", "Тип", "Голова / скачок", "Хвост"],
             bbox=(0, .31, 1, .20), cellLoc="center")
    states = ["(" + ", ".join(f"{number(case.meta, side+'_'+k):g}" for k in ("rho", "u", "p")) + ")"
              for side in ("left", "right")]
    ax.text(0, .88, f"WL={states[0]}; WR={states[1]}; x0={case.meta['x0']}; точек: {len(case.x)}.", fontsize=10)
    ax.text(0, .24, f"GUESSP: {GUESS[case.meta['guess']]}; p0={number(case.meta, 'initial_pressure'):.9g}; итераций: {case.meta['iterations']}.")
    for i, (side, shock, _, _) in enumerate(waves(case)):
        ax.text(0, .17-i*.05, f"{'Слева' if side == 'left' else 'Справа'} p* {'>' if shock else '<'} pK: {'ударная волна' if shock else 'волна разрежения'}.")
    ax.text(0, .02, "Допуск: atol + 2e-6·|эталон|. Прочерк: относительное отклонение при нулевом эталоне не определено.", fontsize=9)
    return fig


def render(cases, output):
    output.mkdir(parents=True, exist_ok=True)
    report, summary = ["# Точные профили Торо" if cases[0].direct_exact else "# Тесты Торо", ""], []
    if cases[0].exact:
        report += ["Эталоны: Торо, таблица 4.3. Допуск: atol + 2e-6·|эталон|.",
                   "atol учитывает округление таблицы и входов теста 5. При нулевом эталоне относительное отклонение не определено (—).", ""]
    with PdfPages(output / ("toro_report.pdf" if cases[0].exact else "toro_figures.pdf")) as pdf:
        for case in cases:
            directory = output / case.name
            directory.mkdir(parents=True, exist_ok=True)
            drawings = [("profiles", profiles)]
            if not case.direct_exact:
                drawings.append(("xt", xt_diagram))
            if {"pressure", "velocity"} <= case.values.keys():
                drawings.append(("pu", pu_diagram))
            report += report_section(case)
            if case.exact:
                fig = results_page(case)
                pdf.savefig(fig)
                plt.close(fig)
            for name, draw in drawings:
                fig = draw(case)
                fig.savefig(directory / f"{name}.png", dpi=150)
                fig.savefig(directory / f"{name}.pdf")
                pdf.savefig(fig)
                plt.close(fig)
                report += ["", f"![{name}]({case.name}/{name}.png)"]
            if "pu" not in [name for name, _ in drawings]:
                report += ["Диаграмма p–u пропущена: недостаточно полей."]
            if case.exact:
                summary += [[case.name]+[r[k] for k in ERROR_COLUMNS] for r in case.errors]
            else:
                summary += [[case.name, k, v[-1].min(), v[-1].max()] for k, v in case.values.items()]
            report += [""]
            print(f"{case.name}: {len(drawings)} рисунка")
    with (output / "toro_summary.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["case"]+ERROR_COLUMNS if cases[0].exact else ["case", "field", "min_final", "max_final"])
        writer.writerows(summary)
    (output / "REPORT.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Готово: {output.resolve()}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.grid": True,
                         "grid.color": "#e5eaf0", "axes.spines.top": False, "axes.spines.right": False})
    try:
        cases = [read_case(p) for p in discover_cases(args.results)]
        render(cases, args.output_dir or args.results/"plots")
    except (OSError, ValueError, KeyError, FloatingPointError) as error:
        parser.exit(1, f"Ошибка: {error}\n")


if __name__ == "__main__":
    main()
