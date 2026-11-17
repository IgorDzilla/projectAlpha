# Точная задача Римана для идеального газа

`solver-cli` читает два постоянных состояния газа слева и справа от разрыва,
один раз решает задачу Римана и отбирает точный профиль в заданный момент времени.
Он не выполняет шаги по времени, не строит сеточную схему и не использует CFL.
Слово «точный» относится к автомодельному решению; давление звёздной области
находится численным поиском корня с заданным допуском.

## Быстрый запуск

Из каталога `project Alpha`:

```sh
cmake -S numerical-solver -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target solver-cli
./build/solver-cli configs/toro_1.toml --output-dir results/sod-1
python3 plot_toro.py results/sod-1
```

Каталог результата должен быть новым или пустым. `solution.csv` содержит
начальные значения и точный профиль при `domain.t[1]` в `grid.points` точках.
`run.txt` содержит найденные `p_star`, `u_star`, типы и скорости волн. Графики
записываются в `results/sod-1/plots/`. Для Python нужны NumPy и Matplotlib.

[Подробная инструкция](docs/running.md).

## Рабочий путь

```text
configs/toro_1.toml
  -> main.cpp -> EulerConfig.cpp -> EulerRun.cpp
  -> ExactEulerRiemannSolver::solve(left, right)       один раз
  -> ExactEulerRiemannSolver::sample(solution, x/t)   для точек профиля
  -> solution.csv + run.txt -> plot_toro.py
```

`EulerEquations` переводит между примитивными и консервативными переменными.
Число `grid.points` определяет только плотность выборки профиля; оно не меняет
само решение Римана. `domain.t` задаёт начальный момент и момент просмотра:
время нужно для положения движущихся волн, но не для интегрирования по времени.

В исходном коде остаются `FiniteVolume`, `ExplicitEuler` и `ExplicitSolver` как
отдельные реализации и учебные проверки. `solver-cli` их не вызывает.
`utils/ConfigProcessor` читает другой, волновой формат TOML v1 и также не
участвует в этом запуске.

Пять точных задач Торо и свойства волн проверяет `riemann-test`:

```sh
cmake -S numerical-solver -B build -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

[Проверки и графики по всем пяти задачам](tests/README.md).
