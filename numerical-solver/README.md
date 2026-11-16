# Решатель Римана

Программа считает, что произойдёт с газом, если слева и справа от одной точки задать разные плотность, скорость и давление.

Ты задаёшь числа в TOML-файле. Программа считает движение газа и записывает результат в CSV. Отдельный Python-скрипт строит графики.

Сейчас работает одномерная задача с идеальным газом. Метод — Годунов первого порядка. Расчёт идёт последовательно.

## Быстрый запуск

Все команды ниже выполняются из **project Alpha**:

```sh
cd "/Users/kozelskijigor/Documents/diffSchemes/project Alpha"
cmake -S numerical-solver -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target solver-cli
./build/solver-cli configs/toro_1.toml --output-dir results/sod-1
python3 plot_toro.py results/sod-1
```

Что получится:

- `results/sod-1/solution.csv` — значения по ячейкам в разные моменты времени.
- `results/sod-1/run.txt` — параметры расчёта и отметка об успешном завершении.
- `results/sod-1/plots/` — картинки, общий PDF и короткий отчёт.

Если папка `sod-1` уже заполнена, укажи другое имя, например `sod-2`. Решатель не затирает старые результаты.

[Подробно: запуск, настройки и частые ошибки](docs/running.md).

## Какие конфиги читаются

| Что | Чем читается | Можно запустить расчёт? |
| --- | --- | --- |
| `configs/toro_*.toml`, версия 2, газ | `EulerConfig` внутри CLI | Да |
| Волновой TOML версии 1 | Твой `utils/ConfigProcessor` | Прочитать можно; сам волновой расчёт пока не написан |

ExprTk находится в `utils` и вычисляет формулы для его ридера. Газовые конфиги CLI содержат обычные числа и ExprTk не используют. [Коротко про utils](../utils/README.md).

## Где какой код

| Файл | Что делает |
| --- | --- |
| [main.cpp](apps/solver-cli/main.cpp) | Получает путь к конфигу из командной строки |
| [EulerConfig.cpp](apps/solver-cli/config/EulerConfig.cpp) | Читает газовый TOML и проверяет параметры |
| [EulerRun.cpp](apps/solver-cli/EulerRun.cpp) | Создаёт нужные объекты и запускает расчёт |
| [RiemannProblem1D.cpp](src/problems/RiemannProblem1D.cpp) | Раскладывает начальные данные по ячейкам |
| [EulerEquations.cpp](src/models/EulerEquations.cpp) | Формулы газа: давление, энергия, скорость звука |
| [ExactEulerRiemannSolver.cpp](src/methods/riemann/ExactEulerRiemannSolver.cpp) | Решает один разрыв между двумя состояниями газа |
| [GodunovFlux.hpp](include/numerics/methods/fluxes/GodunovFlux.hpp) | Вычисляет перенос массы, импульса и энергии через одну границу ячеек |
| [FiniteVolume.cpp](src/methods/FiniteVolume.cpp) | Собирает эти потоки по всей сетке и ограничивает шаг времени |
| [ExplicitEuler.hpp](include/numerics/integrators/ExplicitEuler.hpp) | Считает новые значения за один шаг |
| [ExplicitSolver.hpp](include/numerics/solvers/ExplicitSolver.hpp) | Повторяет шаги до заданного времени |
| [CsvResultWriter.cpp](src/io/CsvResultWriter.cpp) | Записывает результат в CSV |
| [plot_toro.py](../plot_toro.py) | Строит графики из результата |

`Mesh` — сетка. `Field` — значения одной величины на сетке. `State` — все поля расчёта и его текущее время. `Solver` — общий интерфейс решателя; газовый расчёт использует `ExplicitSolver`.

В расчёте газа хранятся средние значения по ячейкам. `Field` также умеет хранить значения в центрах и узлах: `CellCenter`, `Node`. У 10 ячеек 11 узлов. Выбор варианта хранения сам ничего не усредняет и не пересчитывает.

Файлы других моделей и методов оставлены пустыми. Их наличие не означает, что эти методы уже работают.

[Как проверить программу](tests/README.md).
