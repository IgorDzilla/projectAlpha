#include "EulerRun.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <locale>
#include <stdexcept>

namespace {
void writeProfile(std::ofstream& csv, const EulerConfig& c,
                  const EulerEquations& gas, const ExactEulerRiemannSolver& solver,
                  const RiemannSolution& solution, double time) {
  const double dx = (c.xMax - c.xMin) / static_cast<double>(c.points);
  for (std::size_t i = 0; i < c.points; ++i) {
    const double x = std::fma(static_cast<double>(i) + 0.5, dx, c.xMin);
    const EulerPrimitiveState w = time == c.tStart
        ? (x < c.discontinuity ? c.left : c.right)
        : solver.sample(solution, (x - c.discontinuity) / (time - c.tStart));
    const auto u = gas.toConservative(w);
    csv << time << ',' << i << ',' << x;
    for (const auto& name : c.output->fields) {
      if (name == "density") csv << ',' << u.density;
      else if (name == "momentum") csv << ',' << u.momentum;
      else if (name == "energy") csv << ',' << u.energy;
      else if (name == "velocity") csv << ',' << w.velocity;
      else if (name == "pressure") csv << ',' << w.pressure;
    }
    csv << '\n';
  }
}
}

EulerRunResult runEuler(const EulerConfig& c) {
  if (!std::isfinite(c.tStart) || !std::isfinite(c.tEnd) || c.tEnd < c.tStart ||
      !std::isfinite(c.tEnd - c.tStart) || c.points == 0 ||
      !std::isfinite(c.xMin) || !std::isfinite(c.xMax) || c.xMax <= c.xMin ||
      !(c.discontinuity > c.xMin && c.discontinuity < c.xMax))
    throw std::invalid_argument("Invalid exact Riemann domain, time, or sample count");

  const EulerEquations gas(c.gamma);
  const ExactEulerRiemannSolver solver(c.gamma, c.exact);
  // Один раз решаем исходный разрыв. Все профили — выборки одного решения.
  const auto solution = solver.solve(c.left, c.right);
  if (!c.output) return {c.tEnd, 0, solution.pStar, solution.uStar};

  const auto& directory = c.output->directory;
  if (directory.empty()) throw std::invalid_argument("Output directory is empty");
  if (std::filesystem::exists(directory) && !std::filesystem::is_empty(directory))
    throw std::runtime_error("Output directory is not empty: " + directory.string());
  std::filesystem::create_directories(directory);

  std::ofstream csv(directory / "solution.csv");
  csv.exceptions(std::ios::failbit | std::ios::badbit);
  csv.imbue(std::locale::classic());
  csv << std::setprecision(17) << "time,point,x";
  for (const auto& name : c.output->fields) csv << ',' << name;
  csv << '\n';
  writeProfile(csv, c, gas, solver, solution, c.tStart);
  const std::size_t snapshots = c.tEnd > c.tStart ? 2 : 1;
  if (snapshots == 2) writeProfile(csv, c, gas, solver, solution, c.tEnd);
  csv.flush();

  std::ofstream meta(directory / "run.txt");
  meta.exceptions(std::ios::failbit | std::ios::badbit);
  meta.imbue(std::locale::classic());
  meta << std::setprecision(17)
       << "model=euler_1d\nmethod=exact_riemann\nrepresentation=point_samples\n"
       << "gamma=" << c.gamma << "\npoints=" << c.points
       << "\nx_min=" << c.xMin << "\nx_max=" << c.xMax
       << "\nx0=" << c.discontinuity
       << "\nt_start=" << c.tStart << "\nt_end=" << c.tEnd
       << "\nleft_rho=" << c.left.density << "\nleft_u=" << c.left.velocity
       << "\nleft_p=" << c.left.pressure
       << "\nright_rho=" << c.right.density << "\nright_u=" << c.right.velocity
       << "\nright_p=" << c.right.pressure
       << "\np_star=" << solution.pStar << "\nu_star=" << solution.uStar
       << "\nrho_star_left=" << solution.starLeft.density
       << "\nrho_star_right=" << solution.starRight.density
       << "\nleft_type=" << toString(solution.leftWave.type)
       << "\nleft_head=" << solution.leftWave.headSpeed
       << "\nleft_tail=" << solution.leftWave.tailSpeed
       << "\nright_type=" << toString(solution.rightWave.type)
       << "\nright_head=" << solution.rightWave.headSpeed
       << "\nright_tail=" << solution.rightWave.tailSpeed
       << "\ncontact=" << solution.contactSpeed
       << "\niterations=" << solution.iterations
       << "\nsnapshots=" << snapshots << "\nfinal_time=" << c.tEnd
       << "\nstatus=complete\n";
  meta.flush();
  return {c.tEnd, snapshots, solution.pStar, solution.uStar};
}
