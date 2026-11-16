#pragma once

#include "numerics/methods/riemann/ExactEulerRiemannSolver.hpp"
#include "toml.hpp"

#include <filesystem>
#include <optional>
#include <vector>
#include <string>

struct EulerOutputConfig {
  std::filesystem::path directory;
  std::size_t saveEvery = 1;
  std::vector<std::string> fields{"density", "momentum", "energy", "velocity", "pressure"};
};

// Проверенный профиль TOML v2. В численную библиотеку TOML не передаётся.
struct EulerConfig {
  double gamma = 1.4;
  double xMin = 0.0, xMax = 1.0, tStart = 0.0, tEnd = 0.0;
  std::size_t cells = 0;
  double discontinuity = 0.5;
  EulerPrimitiveState left, right;
  double cfl = 0.4;
  ExactRiemannOptions exact;
  std::size_t maxSteps = 1000000;
  std::optional<EulerOutputConfig> output;
};

EulerConfig parseEulerConfig(const toml::table& table,
                             const std::filesystem::path& baseDirectory = ".");
EulerConfig readEulerConfig(const std::filesystem::path& path);
