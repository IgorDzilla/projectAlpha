# Numerical Solver

Каркас библиотеки численных методов и консольного приложения для решения волновой задачи, переноса, теплопроводности и газовой динамики, включая одномерный распад разрыва (задачу Римана).

**Работает полный последовательный расчёт газа в 1D:** TOML v2 → начальные средние → Годунов и конечные объёмы → явный Эйлер → CSV. Есть CLI; пять конфигов Торо проходят полный запуск. [Инструкция запуска](docs/running.md). Остальные физические профили пока не реализованы, старый ридер волнового TOML v1 сохранён. CTest запускает 12 проверок.

**Параллелизация отложена до явной отмашки пользователя.** Каталог `parallel/` резервирует место для будущих интерфейсов. При первой реализации допускается только последовательное выполнение; подключения OpenMP, MPI, потоков и GPU пока не планируются.

В соседнем [utils](../utils/README.md) уже есть пользовательский конфиг-ридер для `wave_1d`. Его место в системе и ограничения описаны в [документе интеграции](docs/configuration.md). Код волнового ридера сохранён на месте. Для газового TOML v2 добавлен отдельный EulerConfig в слое приложения, использующий тот же toml.hpp.

## Быстрый запуск

Из каталога `numerical-solver`:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build build --target solver-cli
./build/solver-cli ../configs/toro_1.toml --output-dir /tmp/sod-new-run
```

Каталог вывода должен быть новым или пустым. Результат: `solution.csv` и `run.txt`.

## Документация

- [Запуск и выходные файлы](docs/running.md).
- [Как соединены все цепочки расчёта](docs/solver-chain.md).
- [Архитектура и ответственность модулей](docs/architecture.md).
- [Одномерный распад разрыва: постановка, метод и конфигурация](docs/riemann-1d.md).
- [Как вызвать адаптированный решатель Римана](docs/riemann-api.md).
- [Точный решатель Римана и поток Годунова — контракт](docs/exact-riemann-solver.md).
- [Интеграция конфиг-ридера](docs/configuration.md).
- [Этапы реализации и открытые решения](docs/roadmap.md).
- [План проверки численных результатов](tests/README.md).
- [План измерения производительности](benchmarks/README.md).
- [Размещение будущих реализаций](src/README.md).

## Структура

```text
numerical-solver/
├── CMakeLists.txt
├── apps/solver-cli/
│   ├── main.cpp
│   ├── EulerRun.hpp, EulerRun.cpp
│   └── config/EulerConfig.hpp, EulerConfig.cpp
├── include/numerics/
│   ├── core/                  # Mesh, Field, State
│   ├── models/                # EulerEquations, EulerState
│   ├── problems/              # RiemannProblem1D
│   ├── methods/               # SpatialOperator, FiniteVolume
│   │   ├── fluxes/            # GodunovFlux
│   │   └── riemann/           # ExactEulerRiemannSolver
│   ├── integrators/           # TimeIntegrator, ExplicitEuler
│   ├── solvers/               # Solver, ExplicitSolver
│   ├── io/                    # ResultWriter, CsvResultWriter
│   └── parallel/              # Пустые резервы
├── src/                       # Реализации нетемплейтных компонентов
├── tests/
├── docs/
└── benchmarks/
```

Соседний каталог `../utils/` содержит `toml-processor.hpp`, `toml-processor.cpp`, `math-expression.cpp`, зависимости `toml.hpp` и `exprtk.hpp`, лицензию, спецификацию, пример и тесты. Его код не дублируется в этом дереве.

Модели кроме EulerEquations, Upwind/CentralDifference, RusanovFlux, общий BoundaryCondition, волновой шаг, RK4, LinearSolver и параллельные компоненты остаются заготовками. Газовый FiniteVolume сейчас поддерживает обе границы transmissive напрямую.

## Сборка и проверка

Из этого каталога:

```sh
cmake -S . -B /tmp/numerical-solver-build -DBUILD_TESTING=ON
cmake --build /tmp/numerical-solver-build
ctest --test-dir /tmp/numerical-solver-build --output-on-failure
```

Для сборки нужен соседний каталог `../utils` с существующим ридером и toml.hpp. Цели `numerics` и `euler-app` разделяют численную библиотеку и сборку газовой задачи; `solver-cli` связывается с `euler-app`. CTest проверяет также исходный `config-reader`. Подключение численной библиотеки: `target_link_libraries(my-app PRIVATE numerics)`.

## Mesh и Field

`CartesianMesh<Dim>` задаёт равномерную сетку ячеек. Первая ось меняется быстрее остальных при переводе индексов в один номер. Доступны центры, объём ячейки, размеры сетки, шаги, границы области, координаты граней и узлов. `nodeCount()` возвращает число узлов, `nodePosition(id)` — координаты узла. Для многомерных индексов узлов используются `flattenNode`/`unflattenNode`, для ячеек — прежние `flatten`/`unflatten`.

`Field<T, Dim>` поддерживает три варианта:

| FieldLocation | Что хранится | Число значений |
| --- | --- | --- |
| `CellCenter` | Значение функции в центре ячейки | `mesh.cellCount()` |
| `CellAverage` | Среднее функции по объёму ячейки; вариант по умолчанию | `mesh.cellCount()` |
| `Node` | Значение функции в узле — вершине ячейки, включая границы | `mesh.nodeCount()` |

В 1D у N ячеек N+1 узел, в нескольких измерениях число узлов равно произведению `(shape[d]+1)`. Среднее по ячейке и значение в её центре в общем случае различаются. Выбор размещения не выполняет интерполяцию или усреднение. Фиктивные ячейки и узлы пока не добавлены.

```cpp
CartesianMesh<1> mesh({0.0}, {0.1}, {10});
Field<double, 1> density(mesh, 1.0);
density[3] = 0.125;
double x = mesh.cellCenter(3)[0];
double leftFace = mesh.faceCoordinate(0, 3);
density.fill(1.0);

Field<double, 1> displacement(mesh, 0.0, FieldLocation::Node);
displacement[mesh.nodeCount() - 1] = 2.0; // Правый граничный узел.
auto rightEnd = mesh.nodePosition(mesh.nodeCount() - 1);
```

`State` также принимает `FieldLocation::Node`: все его компоненты и рабочие буферы получают узловой размер. В одном State размещение у всех компонент одинаково.

`at()` и `operator[]` проверяют индекс. Копия поля получает собственные значения, но ссылается на ту же сетку. Присваивание копирует и значения, и привязку к сетке. Сетка должна жить дольше всех связанных полей; менять её присваиванием или перемещать при живых полях нельзя. Временная сетка в конструкторе поля запрещена.

`isCompatibleWith()` требует тот же объект сетки, одинаковое размещение и размер; отдельные объекты сеток автоматически не приравниваются, даже если их параметры совпадают. Поле не проверяет физический смысл значений: положительность плотности и давления остаётся задачей модели.
