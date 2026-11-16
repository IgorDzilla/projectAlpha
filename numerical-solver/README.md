# Газодинамический комплекс

Последовательный расчёт одномерных уравнений Эйлера для идеального газа: метод Годунова первого порядка, точный решатель Римана на гранях и явный Эйлер по времени. Начальные и граничные условия задаются TOML-файлами; результат записывается в CSV, графики строятся на Python.

[Проект комплекса по пунктам 2.1–2.10](docs/architecture.md) описывает текущие модули, их интерфейсы, хранение данных, сборку и тестирование, а также развитие до других потоков, 2D и MPI/OpenMP. Реализованные части и запланированные расширения указаны отдельно.

## Быстрый запуск

Все команды ниже выполняются из корня репозитория `projectAlpha`:

```sh
cmake -S numerical-solver -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build
./build/solver-cli configs/toro_1.toml --output-dir results/sod-1
python3 plot_toro.py results/sod-1
```

Результаты:

- `results/sod-1/solution.csv` — значения по ячейкам в разные моменты времени.
- `results/sod-1/run.txt` — параметры расчёта и отметка об успешном завершении.
- `results/sod-1/plots/` — картинки, общий PDF и короткий отчёт.

Папка результата должна быть новой или пустой. Для повторного запуска укажи другое имя.

[Подробно: запуск, настройки и частые ошибки](docs/running.md).

## Какие конфиги читаются

| Что | Чем читается | Можно запустить расчёт? |
| --- | --- | --- |
| `configs/toro_*.toml`, версия 2, газ | `EulerConfig` внутри CLI | Да |
| Волновой TOML версии 1 | `utils/ConfigProcessor` | Чтение реализовано; волновой расчёт пока не написан |

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

## Тесты и графики

После сборки запустить пять точных задач Торо одной командой:

```sh
ctest --test-dir build -R '^riemann$' --output-on-failure
```

Пересчитать их и построить таблицы и графики одной командой:

```sh
cmake --build build --target toro-plots
```

Для графиков нужны Python с NumPy и Matplotlib. Цель доступна при найденном CMake интерпретаторе. Условия берутся из общих `configs/toro_*.toml`; данные и отчёт находятся в `build/toro-results`, PDF — в `plots/toro_report.pdf`.

Все 11 проверок C++: `ctest --test-dir build --output-on-failure`. [Проверки Python и допуски](tests/README.md), [форматы результатов и подробные команды](docs/running.md).
