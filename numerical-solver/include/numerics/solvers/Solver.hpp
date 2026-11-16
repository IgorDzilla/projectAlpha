#pragma once

#include "../core/State.hpp"

#include <cstddef>
#include <functional>

struct RunOptions {
  std::size_t saveEvery = 1;
  std::size_t maxSteps = 1000000;
};

// Общий интерфейс расчёта. Формулы схем в него не входят.
template <std::size_t Dim> class Solver {
public:
  using Observer = std::function<void(const State<Dim> &)>;

  virtual ~Solver() = default;
  virtual const State<Dim> &state() const noexcept = 0;
  virtual double maxTimeStep() const = 0;
  virtual void step(double dt) = 0;
  virtual void run(double endTime, const Observer &observer = {},
                   const RunOptions &options = {}) = 0;
};
