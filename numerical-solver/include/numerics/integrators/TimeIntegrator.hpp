#pragma once

#include "../methods/SpatialOperator.hpp"

// Интерфейс явных интеграторов системы первого порядка по времени.
template <std::size_t Dim> class TimeIntegrator {
public:
  virtual ~TimeIntegrator() = default;

  virtual double maxTimeStep(const State<Dim> &state,
                             const SpatialOperator<Dim> &spatial) const = 0;

  // Вернуть новые поля, не меняя вход. Время принимает Solver после проверки.
  virtual typename State<Dim>::Fields advance(
      const State<Dim> &state, const SpatialOperator<Dim> &spatial,
      double dt) const = 0;
};
