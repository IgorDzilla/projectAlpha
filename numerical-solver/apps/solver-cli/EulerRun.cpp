#include "EulerRun.hpp"
#include "numerics/problems/RiemannProblem1D.hpp"
#include "numerics/methods/FiniteVolume.hpp"
#include "numerics/integrators/ExplicitEuler.hpp"
#include "numerics/solvers/ExplicitSolver.hpp"
#include "numerics/io/CsvResultWriter.hpp"

#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>

EulerRunResult runEuler(const EulerConfig& c) {
  if (!std::isfinite(c.tStart) || !std::isfinite(c.tEnd) || c.tEnd < c.tStart ||
      !std::isfinite(c.tEnd-c.tStart) || c.cells == 0 || c.maxSteps == 0 ||
      (c.output && c.output->saveEvery == 0))
    throw std::invalid_argument("Invalid Euler run interval, cell count or step/output limit");
  CartesianMesh<1> mesh({c.xMin}, {(c.xMax-c.xMin)/static_cast<double>(c.cells)}, {c.cells});
  RiemannProblem1D problem(c.left,c.right,c.discontinuity,c.gamma);
  auto initial=problem.initialState(mesh,c.tStart);
  // Отказ по вакууму/несходимости исходного разрыва — до создания файлов.
  const auto riemann=ExactEulerRiemannSolver(c.gamma,c.exact).solve(c.left,c.right);
  ExplicitSolver<1> solver(std::move(initial),
      std::make_unique<FiniteVolume>(mesh,c.gamma,c.cfl,c.exact),
      std::make_unique<ExplicitEuler<1>>(),std::numeric_limits<double>::max());

  std::unique_ptr<CsvResultWriter> writer;
  std::ofstream metadata;
  if (c.output) {
    writer=std::make_unique<CsvResultWriter>(c.output->directory,c.gamma,c.output->fields);
    metadata.exceptions(std::ios::failbit | std::ios::badbit);
    metadata.open(c.output->directory / "run.txt");
    metadata.imbue(std::locale::classic());
    metadata << std::setprecision(17)
             << "model=euler_1d\nmethod=finite_volume\nflux=godunov_exact\nintegrator=explicit_euler\n"
             << "boundaries=transmissive,transmissive\nexecution=sequential\n"
             << "gamma=" << c.gamma << "\ncfl=" << c.cfl << "\ncells=" << c.cells
             << "\nx_min=" << c.xMin << "\nx_max=" << c.xMax << "\ndx=" << mesh.step()[0]
             << "\nt_start=" << c.tStart << "\nt_end=" << c.tEnd
             << "\nx0=" << c.discontinuity
             << "\nleft_rho_u_p=" << c.left.density << ',' << c.left.velocity << ',' << c.left.pressure
             << "\nright_rho_u_p=" << c.right.density << ',' << c.right.velocity << ',' << c.right.pressure
             << "\npressure_rtol=" << c.exact.pressureRtol << "\nresidual_rtol=" << c.exact.residualRtol
             << "\nmax_iterations=" << c.exact.maxIterations
             << "\nmax_bracket_expansions=" << c.exact.maxBracketExpansions
             << "\nriemann_iterations=" << riemann.iterations
             << "\nriemann_pressure_change=" << riemann.relativePressureChange
             << "\nmax_steps=" << c.maxSteps << "\nsave_every=" << c.output->saveEvery
             << "\nrepresentation=conservative_cell_averages; velocity/pressure derived from U\n";
    metadata.flush();
  }
  RunOptions options;
  options.maxSteps=c.maxSteps;
  options.saveEvery=c.output ? c.output->saveEvery : 1;
  Solver<1>::Observer observer;
  if (writer) observer=[&](const State<1>& state) { writer->write(state); };
  try {
    solver.run(c.tEnd,observer,options);
  } catch (const std::exception& e) {
    if (metadata.is_open()) {
      // Сбой записи диагностики не должен скрывать исходную причину отказа.
      try { metadata << "status=failed\nlast_accepted_time=" << solver.state().time() << '\n'; metadata.flush(); }
      catch (const std::exception&) {}
    }
    throw std::runtime_error(std::string("Euler run failed: ")+e.what());
  }
  if (metadata.is_open()) {
    metadata << "status=complete\nfinal_time=" << solver.state().time()
             << "\nsnapshots=" << writer->snapshots() << '\n';
    metadata.flush();
  }
  return {solver.state().time(),writer ? writer->snapshots() : 0};
}
