# Как проверить программу

Все команды — из **project Alpha**.

## Проверки C++

```sh
cmake -S numerical-solver -B build-tests -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON
cmake --build build-tests
ctest --test-dir build-tests --output-on-failure
```

Если всё хорошо, в конце будет `100% tests passed`.

Сейчас есть 11 проверок:

| Название | Что проверяет |
| --- | --- |
| `core` | Сетку и хранение значений |
| `solver` | Шаги времени, сохранение слоёв и остановку при ошибке |
| `riemann` | Локальный решатель: волны, контакты, допустимость данных |
| `finite-volume` | Начальные данные, баланс потоков и уменьшение ошибки при более мелкой сетке |
| `euler-config` | Чтение газового TOML, ошибки параметров и запись CSV |
| `cli-toro-1` … `cli-toro-5` | Полный запуск каждого из пяти TOML и проверку результата |
| `config-reader` | Твой utils: волновой TOML v1 и формулы через ExprTk |

Данные Торо берутся из `configs/toro_*.toml`. Отдельных копий пяти задач в C++ нет. Маленькие искусственные примеры проверяют отдельные свойства, например неизменность постоянного состояния или отказ при отрицательной плотности.

Результаты тестов лежат внутри `build-tests`. Твои расчёты в `configs/results` они не трогают.

Запустить только проверку utils:

```sh
ctest --test-dir build-tests -R '^config-reader$' --output-on-failure
```

## Проверки графиков

Нужны NumPy и Matplotlib.

```sh
python3 -m unittest discover -s numerical-solver/tests -p test_plot_toro.py -v
```

Здесь 9 проверок. Они проверяют чтение CSV, восстановление полей, построение графиков и отказ при испорченных данных. В конце должно быть `OK`.

Если используешь окружение `.venv`, вместо `python3` укажи `.venv/bin/python`.
