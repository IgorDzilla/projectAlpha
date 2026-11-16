#pragma once

#include "TimeIntegrator.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

template <std::size_t Dim> class ExplicitEuler final : public TimeIntegrator<Dim> {
public:
  double maxTimeStep(const State<Dim> &state,
                     const SpatialOperator<Dim> &spatial) const override {
    const double limit = spatial.maxExplicitEulerStep(state);
    if (std::isnan(limit) || limit <= 0.0)
      throw std::runtime_error("Invalid explicit Euler step limit");
    return limit;
  }

  typename State<Dim>::Fields advance(
      const State<Dim> &state, const SpatialOperator<Dim> &spatial,
      double dt) const override {
    state.validate();
    spatial.validate(state);
    if (!std::isfinite(dt) || dt <= 0.0 || dt > maxTimeStep(state, spatial))
      throw std::invalid_argument("Invalid explicit Euler time step");

    // NaN позволяет обнаружить забытые оператором компоненты или ячейки.
    auto rates = state.makeFields();
    for (auto &field : rates)
      field.fill(std::numeric_limits<double>::quiet_NaN());
    spatial.evaluate(state, rates);
    state.validateFields(rates);

    auto result = state.fields();
    for (std::size_t c = 0; c < result.size(); ++c)
      for (std::size_t i = 0; i < result[c].size(); ++i)
        result[c][i] = state.component(c)[i] + dt * rates[c][i];
    state.validateFields(result);
    return result;
  }
};
