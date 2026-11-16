#include "numerics/models/EulerEquations.hpp"

#include <cmath>
#include <stdexcept>

EulerEquations::EulerEquations(double gamma) : gamma_(gamma) {
  if (!std::isfinite(gamma_) || gamma_ <= 1.0) {
    throw std::invalid_argument("gamma must be finite and greater than one");
  }
}

void EulerEquations::validate(const EulerPrimitiveState& state) const {
  if (!std::isfinite(state.density) || state.density <= 0.0 ||
      !std::isfinite(state.velocity) || !std::isfinite(state.pressure) ||
      state.pressure <= 0.0) {
    throw std::invalid_argument("Euler state requires finite velocity and positive finite density/pressure");
  }
}

EulerConservativeState EulerEquations::toConservative(
    const EulerPrimitiveState& state) const {
  validate(state);
  const double momentum = state.density * state.velocity;
  const double energy = state.pressure / (gamma_ - 1.0) +
                        0.5 * momentum * state.velocity;
  if (!std::isfinite(momentum) || !std::isfinite(energy) || energy <= 0.0) {
    throw std::runtime_error("Euler conservative conversion overflow/underflow");
  }
  return {state.density, momentum, energy};
}

EulerPrimitiveState EulerEquations::toPrimitive(
    const EulerConservativeState& state) const {
  if (!std::isfinite(state.density) || state.density <= 0.0 ||
      !std::isfinite(state.momentum) || !std::isfinite(state.energy)) {
    throw std::invalid_argument("invalid Euler conservative state");
  }
  const double velocity = state.momentum / state.density;
  const double pressure = (gamma_ - 1.0) *
                         (state.energy - 0.5 * state.momentum * velocity);
  const EulerPrimitiveState result{state.density, velocity, pressure};
  validate(result); // Отрицательную внутреннюю энергию не исправляем подрезкой.
  return result;
}

double EulerEquations::soundSpeed(const EulerPrimitiveState& state) const {
  validate(state);
  const double speed = std::sqrt(gamma_) * std::sqrt(state.pressure) /
                       std::sqrt(state.density);
  if (!std::isfinite(speed) || speed <= 0.0) {
    throw std::runtime_error("Euler sound speed overflow/underflow");
  }
  return speed;
}

EulerFlux EulerEquations::physicalFlux(const EulerPrimitiveState& state) const {
  const auto conserved = toConservative(state);
  const EulerFlux flux{conserved.momentum,
                       conserved.momentum * state.velocity + state.pressure,
                       (conserved.energy + state.pressure) * state.velocity};
  if (!std::isfinite(flux.mass) || !std::isfinite(flux.momentum) ||
      !std::isfinite(flux.energy)) {
    throw std::runtime_error("Euler physical flux overflow");
  }
  return flux;
}
