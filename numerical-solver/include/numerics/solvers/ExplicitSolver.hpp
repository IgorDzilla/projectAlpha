#pragma once

#include "Solver.hpp"
#include "../integrators/TimeIntegrator.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <utility>

// Последовательный расчёт с явным интегратором и переменным допустимым шагом.
template <std::size_t Dim> class ExplicitSolver final : public Solver<Dim> {
  State<Dim> state_;
  std::unique_ptr<SpatialOperator<Dim>> spatial_;
  std::unique_ptr<TimeIntegrator<Dim>> integrator_;
  double stepLimit_;

  void validate(const State<Dim> &state) const {
    state.validate();
    spatial_->validate(state);
  }

  void advance(double dt, double nextTime) {
    if (!std::isfinite(nextTime) || nextTime <= state_.time())
      throw std::runtime_error("Time step does not advance finite time");
    auto candidate = state_;
    candidate.replaceFields(integrator_->advance(state_, *spatial_, dt));
    candidate.setTime(nextTime);
    validate(candidate);
    // До этой строки исходные поля и время остаются неизменными.
    state_ = std::move(candidate);
  }

public:
  using Observer = typename Solver<Dim>::Observer;

  ExplicitSolver(State<Dim> initial,
                 std::unique_ptr<SpatialOperator<Dim>> spatial,
                 std::unique_ptr<TimeIntegrator<Dim>> integrator,
                 double stepLimit)
      : state_(std::move(initial)), spatial_(std::move(spatial)),
        integrator_(std::move(integrator)), stepLimit_(stepLimit) {
    if (!spatial_ || !integrator_)
      throw std::invalid_argument("Solver requires a spatial operator and integrator");
    if (!std::isfinite(stepLimit_) || stepLimit_ <= 0.0)
      throw std::invalid_argument("Solver step limit must be positive and finite");
    validate(state_);
  }

  const State<Dim> &state() const noexcept override { return state_; }

  double maxTimeStep() const override {
    validate(state_);
    const double limit = integrator_->maxTimeStep(state_, *spatial_);
    if (std::isnan(limit) || limit <= 0.0)
      throw std::runtime_error("Invalid integrator step limit");
    return std::min(stepLimit_, limit);
  }

  void step(double dt) override {
    if (!std::isfinite(dt) || dt <= 0.0 || dt > maxTimeStep())
      throw std::invalid_argument("Requested time step exceeds the allowed limit");
    advance(dt, state_.time() + dt);
  }

  void run(double endTime, const Observer &observer = {},
           const RunOptions &options = {}) override {
    if (!std::isfinite(endTime) || endTime < state_.time())
      throw std::invalid_argument("End time must be finite and not precede current time");
    if (options.saveEvery == 0 || options.maxSteps == 0)
      throw std::invalid_argument("Output interval and step budget must be positive");
    if (observer)
      observer(state_);

    std::size_t steps = 0;
    while (state_.time() < endTime) {
      if (steps == options.maxSteps)
        throw std::runtime_error("Maximum number of time steps reached");
      const double remaining = endTime - state_.time();
      if (!std::isfinite(remaining))
        throw std::invalid_argument("Time interval cannot be represented");
      const double dt = std::min(maxTimeStep(), remaining);
      // Для последнего шага сохраняем в точности запрошенное конечное время.
      const double nextTime = dt == remaining ? endTime : state_.time() + dt;
      advance(dt, nextTime);
      ++steps;
      if (observer && (steps % options.saveEvery == 0 || state_.time() == endTime))
        observer(state_);
    }
  }
};
