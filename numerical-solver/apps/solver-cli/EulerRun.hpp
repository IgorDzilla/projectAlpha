#pragma once
#include "config/EulerConfig.hpp"

struct EulerRunResult {
  double finalTime;
  std::size_t snapshots;
};

// Собирает и запускает только поддержанный профиль: Euler 1D + Godunov + Euler.
EulerRunResult runEuler(const EulerConfig& config);
