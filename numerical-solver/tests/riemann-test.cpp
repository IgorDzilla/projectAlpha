#include "config/EulerConfig.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace fs = std::filesystem;

struct Reference {
  const char* name;
  std::array<double, 4> values;   // p*, u*, rho*L, rho*R
  std::array<double, 4> rounding; // округление эталона и, для теста 5, входных данных
  RiemannWaveType left;
  RiemannWaveType right;
};

const std::array<Reference, 5> references{{
    {"Sod", {0.30313, 0.92745, 0.42632, 0.26557},
     {5e-6, 5e-6, 5e-6, 5e-6}, RiemannWaveType::Rarefaction, RiemannWaveType::Shock},
    {"Two rarefactions", {0.00189, 0.0, 0.02185, 0.02185},
     {5e-6, 1e-8, 5e-6, 5e-6}, RiemannWaveType::Rarefaction, RiemannWaveType::Rarefaction},
    {"Strong left pressure", {460.894, 19.5975, 0.57506, 5.99924},
     {5e-4, 5e-5, 5e-6, 5e-6}, RiemannWaveType::Rarefaction, RiemannWaveType::Shock},
    {"Strong right pressure", {46.0950, -6.19633, 5.99242, 0.57511},
     {5e-5, 5e-6, 5e-6, 5e-6}, RiemannWaveType::Shock, RiemannWaveType::Rarefaction},
    {"Colliding shocks", {1691.64, 8.68975, 14.2823, 31.0426},
     // uL=19.5975 округлено: точный расчёт этих данных даёт u*=8.6897744116.
     // Вариация входов на половину их разряда даёт u* от 8.6897415 до 8.6898074.
     {5e-3, 5e-5, 5e-5, 5e-5}, RiemannWaveType::Shock, RiemannWaveType::Shock},
}};

std::ofstream csv(const fs::path& path) {
  std::ofstream file;
  file.exceptions(std::ios::failbit | std::ios::badbit);
  file.open(path);
  file << std::setprecision(17);
  return file;
}

double soundSpeed(const EulerPrimitiveState& w, double gamma) {
  return EulerEquations(gamma).soundSpeed(w);
}

void check(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, const char* message) {
  check(std::isfinite(actual) && std::isfinite(expected) &&
        std::abs(actual - expected) <= 2e-8 * std::max({1.0, std::abs(actual), std::abs(expected)}),
        message);
}

void same(const EulerPrimitiveState& actual, const EulerPrimitiveState& expected) {
  near(actual.density, expected.density, "sample density");
  near(actual.velocity, expected.velocity, "sample velocity");
  near(actual.pressure, expected.pressure, "sample pressure");
}

void checkPhysics(const ExactEulerRiemannSolver& solver, const RiemannSolution& r) {
  const EulerEquations gas(solver.gamma());
  for (bool left : {true, false}) {
    const auto& outer = left ? r.left : r.right;
    const auto& star = left ? r.starLeft : r.starRight;
    const auto& wave = left ? r.leftWave : r.rightWave;
    const double sign = left ? -1.0 : 1.0;
    same(solver.sample(r, wave.headSpeed + sign), outer);
    same(solver.sample(r, (wave.tailSpeed + r.contactSpeed) / 2), star);
    if (wave.type == RiemannWaveType::Shock) {
      // Ранкин — Гюгонио: скачок потока равен скорости фронта на скачок U.
      const auto a = gas.toConservative(outer), b = gas.toConservative(star);
      const auto fa = gas.physicalFlux(outer), fb = gas.physicalFlux(star);
      near(fb.mass - fa.mass, wave.headSpeed * (b.density - a.density), "shock mass balance");
      near(fb.momentum - fa.momentum, wave.headSpeed * (b.momentum - a.momentum), "shock momentum balance");
      near(fb.energy - fa.energy, wave.headSpeed * (b.energy - a.energy), "shock energy balance");
      check(star.pressure / std::pow(star.density, gas.gamma()) >=
            outer.pressure / std::pow(outer.density, gas.gamma()), "shock entropy decrease");
    } else {
      const double ray = (wave.headSpeed + wave.tailSpeed) / 2;
      const auto fan = solver.sample(r, ray);
      const double a = gas.soundSpeed(fan), aOuter = gas.soundSpeed(outer);
      // Веер изэнтропичен; сохраняется инвариант, луч совпадает с характеристикой.
      near(fan.pressure / std::pow(fan.density, gas.gamma()),
           outer.pressure / std::pow(outer.density, gas.gamma()), "fan entropy");
      near(fan.velocity - sign * 2 * a / (gas.gamma() - 1),
           outer.velocity - sign * 2 * aOuter / (gas.gamma() - 1), "fan Riemann invariant");
      near(fan.velocity + sign * a, ray, "fan characteristic");
    }
  }
  const auto l = solver.sample(r, (r.leftWave.tailSpeed + r.contactSpeed) / 2);
  const auto right = solver.sample(r, (r.rightWave.tailSpeed + r.contactSpeed) / 2);
  near(l.pressure, right.pressure, "contact pressure continuity");
  near(l.velocity, right.velocity, "contact velocity continuity");
  same(solver.sample(r, r.contactSpeed), r.starLeft);
}

void checkContactAndVacuum(double gamma) {
  const EulerEquations gas(gamma);
  const ExactEulerRiemannSolver solver(gamma);
  for (double u : {-2.0, 0.0, 2.0}) {
    const EulerPrimitiveState left{1, u, 2}, right{3, u, 2};
    const auto contact = solver.solve(left, right);
    same(solver.sample(contact, u - .1), left);
    same(solver.sample(contact, u + .1), right);
    const auto uniform = solver.solve(left, left);
    for (double ray : {u - 10, u, u + 10}) same(solver.sample(uniform, ray), left);
  }
  const auto contact = solver.solve({1, 0, 2}, {3, 0, 2});
  const auto flux = gas.physicalFlux(solver.sample(contact, 0));
  near(flux.mass, 0, "stationary contact mass flux");
  near(flux.momentum, 2, "stationary contact momentum flux");
  near(flux.energy, 0, "stationary contact energy flux");
  const double critical = 2 * gas.soundSpeed({1, 0, 1}) / (gamma - 1);
  const auto expansion = solver.solve({1, -.99 * critical, 1}, {1, .99 * critical, 1});
  const double expected = std::pow(1 - .99, 2 * gamma / (gamma - 1));
  near(expansion.pStar / expected, 1, "near-vacuum analytic pressure");
  check(expansion.starLeft.density > 0 && expansion.starRight.density > 0,
        "near-vacuum density must stay positive");
  for (double u : {critical, 1.01 * critical}) {
    bool vacuum = false;
    try { solver.solve({1, -u, 1}, {1, u, 1}); }
    catch (const std::domain_error&) { vacuum = true; }
    check(vacuum, "vacuum must be detected at and above the expansion limit");
  }
}

double coordinate(double first, double last, std::size_t i, std::size_t count) {
  return first + (last - first) * static_cast<double>(i) / static_cast<double>(count - 1);
}

// Только формула для диаграммы p-u из §15.6; Ньютон здесь не дублируется.
double waveCurve(double p, const EulerPrimitiveState& w, double gamma) {
  if (p > w.pressure)
    return (p - w.pressure) * std::sqrt(2.0 / ((gamma + 1) * w.density) /
                                      (p + (gamma - 1) / (gamma + 1) * w.pressure));
  return 2 * soundSpeed(w, gamma) / (gamma - 1) *
         (std::pow(p / w.pressure, (gamma - 1) / (2 * gamma)) - 1);
}

void exportPlots(const fs::path& directory, const EulerConfig& c,
                 const ExactEulerRiemannSolver& solver, const RiemannSolution& r) {
  auto profile = csv(directory / "profile.csv");
  profile << "x,density,velocity,pressure,internal_energy,sound_speed,entropy\n";
  for (std::size_t i = 0; i < 2001; ++i) {
    const double x = coordinate(c.xMin, c.xMax, i, 2001);
    const auto w = solver.sample(r, (x - c.discontinuity) / (c.tEnd - c.tStart));
    profile << x << ',' << w.density << ',' << w.velocity << ',' << w.pressure << ','
            << w.pressure / ((c.gamma - 1) * w.density) << ','
            << soundSpeed(w, c.gamma) << ',' << w.pressure / std::pow(w.density, c.gamma) << '\n';
  }
  auto xt = csv(directory / "density_xt.csv");
  xt << "time,cell,x,density\n";
  for (std::size_t j = 0; j < 121; ++j) {
    const double t = coordinate(0, (c.tEnd - c.tStart), j, 121);
    for (std::size_t i = 0; i < 301; ++i) {
      const double x = coordinate(c.xMin, c.xMax, i, 301);
      // При t=0 берём начальные данные, в точке разрыва — правое состояние.
      const auto w = j == 0 ? (x < c.discontinuity ? r.left : r.right)
                            : solver.sample(r, (x - c.discontinuity) / t);
      xt << c.tStart + t << ',' << i << ',' << x << ',' << w.density << '\n';
    }
  }
  const double pMin = 0.02 * std::min({r.left.pressure, r.right.pressure, r.pStar});
  const double pMax = 2.0 * std::max({r.left.pressure, r.right.pressure, r.pStar});
  std::vector<double> pressures{r.left.pressure, r.right.pressure, r.pStar};
  for (std::size_t i = 0; i < 801; ++i)
    pressures.push_back(std::exp(coordinate(std::log(pMin), std::log(pMax), i, 801)));
  std::sort(pressures.begin(), pressures.end());
  pressures.erase(std::unique(pressures.begin(), pressures.end()), pressures.end());
  auto curves = csv(directory / "wave_curves.csv");
  curves << "p,u_left,u_right\n";
  for (double p : pressures)
    curves << p << ',' << r.left.velocity - waveCurve(p, r.left, c.gamma) << ','
           << r.right.velocity + waveCurve(p, r.right, c.gamma) << '\n';
}

bool compare(const fs::path& directory, const Reference& reference,
             const RiemannSolution& r) {
  const std::array<double, 4> actual{r.pStar, r.uStar, r.starLeft.density, r.starRight.density};
  const std::array<const char*, 4> names{"p_star", "u_star", "rho_star_left", "rho_star_right"};
  auto errors = csv(directory / "errors.csv");
  errors << "quantity,actual,reference,absolute_error,relative_error,tolerance,passed\n";
  bool passed = true;
  for (std::size_t i = 0; i < actual.size(); ++i) {
    const double error = std::abs(actual[i] - reference.values[i]);
    const double tolerance = reference.rounding[i] + 2e-6 * std::abs(reference.values[i]);
    const bool ok = std::isfinite(actual[i]) && error <= tolerance;
    passed &= ok;
    errors << names[i] << ',' << actual[i] << ',' << reference.values[i] << ',' << error << ',';
    if (reference.values[i] != 0) errors << error / std::abs(reference.values[i]);
    errors << ',' << tolerance << ',' << ok << '\n';
    if (!ok)
      std::cerr << "  " << names[i] << ": expected " << reference.values[i]
                << ", got " << actual[i] << '\n';
  }
  passed &= r.leftWave.type == reference.left && r.rightWave.type == reference.right;
  return passed;
}

void exportSummary(const fs::path& directory, std::size_t id, const Reference& reference,
                   const EulerConfig& c, const RiemannSolution& r, bool passed) {
  auto out = csv(directory / "summary.csv");
  out << "key,value\n";
  out << "test," << id << "\nname," << reference.name << "\npassed," << passed
      << "\ngamma," << c.gamma << "\nx_min," << c.xMin << "\nx_max," << c.xMax
      << "\nx0," << c.discontinuity << "\nt_start," << c.tStart << "\ntout," << c.tEnd
      << "\np_star," << r.pStar << "\nu_star," << r.uStar
      << "\nrho_star_left," << r.starLeft.density << "\nrho_star_right," << r.starRight.density
      << "\nguess," << toString(r.guessType)
      << "\ninitial_pressure," << r.initialPressure << "\niterations," << r.iterations
      << "\npressure_change," << r.relativePressureChange
      << "\nleft_type," << toString(r.leftWave.type)
      << "\nleft_head," << r.leftWave.headSpeed << "\nleft_tail," << r.leftWave.tailSpeed
      << "\ncontact," << r.contactSpeed
      << "\nright_type," << toString(r.rightWave.type)
      << "\nright_head," << r.rightWave.headSpeed << "\nright_tail," << r.rightWave.tailSpeed;
  for (bool left : {true, false}) {
    const auto& w = left ? r.left : r.right;
    const std::string side = left ? "left" : "right";
    out << '\n' << side << "_rho," << w.density << '\n' << side << "_u," << w.velocity
        << '\n' << side << "_p," << w.pressure;
  }
  out << '\n';
}
} // namespace

int main(int argc, char** argv) {
  if (argc > 3) {
    std::cerr << "Usage: riemann-test [config-directory] [output-directory]\n";
    return 1;
  }
  const fs::path configs = argc > 1 ? argv[1] : TORO_CONFIG_DIR;
  const fs::path output = argc > 2 ? argv[2] : TORO_OUTPUT_DIR;
  bool allPassed = true;
  std::cout << std::setprecision(8);
  for (std::size_t i = 0; i < references.size(); ++i) {
    const auto id = std::to_string(i + 1);
    const auto directory = output / ("toro_" + id);
    std::cout << "\nToro " << id << ": " << references[i].name << '\n';
    try {
      fs::create_directories(directory);
      fs::remove(directory / "summary.csv"); // Не оставляем старый PASS после неудачного запуска.
      const auto path = configs / ("toro_" + id + ".toml");
      const auto config = readEulerConfig(path);
      if (config.tEnd <= config.tStart)
        throw std::runtime_error("Toro profiles require a positive duration");
      const ExactEulerRiemannSolver solver(config.gamma, config.exact);
      const auto solution = solver.solve(config.left, config.right);
      if (solution.relativePressureChange >= config.exact.pressureRtol)
        throw std::runtime_error("pressure change exceeds the requested tolerance");
      const bool passed = compare(directory, references[i], solution);
      checkPhysics(solver, solution);
      EulerEquations(config.gamma).physicalFlux(solver.sample(solution, 0.0));
      if (i == 0) {
        checkContactAndVacuum(config.gamma);
        checkContactAndVacuum(5.0 / 3.0);
        auto coarse = config.exact, fine = config.exact;
        coarse.pressureRtol = 1e-3; fine.pressureRtol = 1e-10;
        const auto early = ExactEulerRiemannSolver(config.gamma, coarse).solve(config.left, config.right);
        const auto precise = ExactEulerRiemannSolver(config.gamma, fine).solve(config.left, config.right);
        if (early.relativePressureChange >= coarse.pressureRtol ||
            precise.relativePressureChange >= fine.pressureRtol || early.iterations >= precise.iterations)
          throw std::runtime_error("pressure tolerance does not control convergence");
      }
      exportPlots(directory, config, solver, solution);
      fs::copy_file(path, directory / "config.toml", fs::copy_options::overwrite_existing);
      exportSummary(directory, i + 1, references[i], config, solution, passed);
      std::cout << (passed ? "[PASS]" : "[FAIL]") << " Toro " << id
                << ", GUESSP: " << toString(solution.guessType)
                << ", Iterations: " << solution.iterations << '\n';
      allPassed &= passed;
    } catch (const std::exception& error) {
      std::cerr << "[FAIL] Toro " << id << ": " << error.what() << '\n';
      allPassed = false;
    }
  }
  std::cout << '\n' << (allPassed ? "All five Toro tests passed.\n" : "At least one Toro test failed.\n")
            << "Output: " << fs::absolute(output).string() << '\n';
  return allPassed ? 0 : 1;
}
