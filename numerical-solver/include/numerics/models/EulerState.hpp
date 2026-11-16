#pragma once

#include "EulerEquations.hpp"
#include "../core/State.hpp"

// Единый порядок компонент газодинамического State: rho, rho*u, E.
inline void validateEulerLayout(const State<1>& state) {
  state.validate();
  if (state.location() != FieldLocation::CellAverage ||
      state.names() != std::vector<std::string>{"density", "momentum", "energy"})
    throw std::invalid_argument("Euler requires cell averages: density, momentum, energy");
}

inline EulerConservativeState eulerCell(const State<1>& state, std::size_t i) {
  return {state.component(0)[i], state.component(1)[i], state.component(2)[i]};
}
