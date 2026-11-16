#pragma once

#include "../core/State.hpp"

// Конкретный оператор связывает модель, пространственную схему и границы.
template <std::size_t Dim> class SpatialOperator {
public:
  virtual ~SpatialOperator() = default;

  // Заполнить ВСЕ производные. Границы учитываются при state.time().
  // Физические значения входного состояния менять нельзя.
  virtual void evaluate(const State<Dim> &state,
                        typename State<Dim>::Fields &rates) const = 0;

  // Ограничение именно для явного Эйлера. +infinity означает отсутствие ограничения.
  virtual double maxExplicitEulerStep(const State<Dim> &state) const = 0;

  // Проверки физики: например, положительность плотности и давления.
  // Структуру и конечность чисел Solver проверяет отдельно всегда.
  virtual void validate(const State<Dim> &) const {}
};
