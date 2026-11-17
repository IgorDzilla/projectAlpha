#pragma once
#include "config/EulerConfig.hpp"

struct EulerRunResult {
  double finalTime;
  std::size_t snapshots;
  double pStar;
  double uStar;
};

// Решает один исходный разрыв и отбирает точный профиль в заданное время.
EulerRunResult runEuler(const EulerConfig& config);
