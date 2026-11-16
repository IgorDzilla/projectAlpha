#include "config/EulerConfig.hpp"
#include "numerics/methods/FiniteVolume.hpp"
#include "numerics/problems/RiemannProblem1D.hpp"
#include "numerics/integrators/ExplicitEuler.hpp"
#include "numerics/solvers/ExplicitSolver.hpp"

#include <array>
#include <iostream>

namespace {
void check(bool value, const char* message) { if (!value) throw std::runtime_error(message); }
void near(double a,double b,const char* message,double tol=2e-10) {
  check(std::isfinite(a) && std::isfinite(b) && std::abs(a-b)<=tol*std::max({1.0,std::abs(a),std::abs(b)}),message);
}
template<class F> void fails(F action) {
  bool failed=false;
  try { action(); } catch (const std::exception&) { failed=true; }
  check(failed,"expected failure");
}
std::array<double,3> integral(const State<1>& s) {
  std::array<double,3> result{};
  for (std::size_t c=0;c<3;++c)
    for (std::size_t i=0;i<s.mesh().cellCount();++i)
      result[c]+=s.component(c)[i]*s.mesh().cellVolume(i);
  return result;
}
void projection() {
  const EulerEquations gas;
  const EulerPrimitiveState left{2,3,4},right{1,-2,1};
  const auto l=gas.toConservative(left),r=gas.toConservative(right);
  for (double x0 : {0.5,0.37}) {
    CartesianMesh<1> mesh({0},{0.25},{4});
    auto state=RiemannProblem1D(left,right,x0).initialState(mesh,2);
    const auto total=integral(state);
    near(total[0],x0*l.density+(1-x0)*r.density,"projection mass");
    near(total[1],x0*l.momentum+(1-x0)*r.momentum,"projection momentum");
    near(total[2],x0*l.energy+(1-x0)*r.energy,"projection energy");
    check(state.time()==2,"initial time");
    FiniteVolume op(mesh);
    op.validate(state);
  }
}
void constantAndContact() {
  for (std::size_t n : {std::size_t(1),std::size_t(10)}) {
    CartesianMesh<1> mesh({0},{1.0/n},{n});
    for (double velocity : {-4.0,0.0,4.0}) {
      const EulerPrimitiveState w{1.3,velocity,2.1};
      auto state=RiemannProblem1D(w,w,0.5).initialState(mesh);
      const auto initial=state;
      ExplicitSolver<1> solver(state,std::make_unique<FiniteVolume>(mesh),
                              std::make_unique<ExplicitEuler<1>>(),1.0);
      solver.run(0.1);
      for (std::size_t c=0;c<3;++c)
        for (std::size_t i=0;i<n;++i)
          near(solver.state().component(c)[i],initial.component(c)[i],"constant state changed");
    }
  }
  CartesianMesh<1> mesh({0},{0.05},{20});
  auto state=RiemannProblem1D({1,0,1},{0.125,0,1},0.5).initialState(mesh);
  const auto before=state;
  ExplicitSolver<1> solver(state,std::make_unique<FiniteVolume>(mesh),
                          std::make_unique<ExplicitEuler<1>>(),1);
  solver.run(0.15);
  for (std::size_t c=0;c<3;++c)
    for (std::size_t i=0;i<20;++i)
      near(solver.state().component(c)[i],before.component(c)[i],"stationary contact changed");
}
void balanceAndCache(const EulerConfig& config) {
  const double dx=(config.xMax-config.xMin)/40;
  CartesianMesh<1> mesh({config.xMin},{dx},{40});
  auto state=RiemannProblem1D(config.left,config.right,config.discontinuity,config.gamma)
      .initialState(mesh,config.tStart);
  const EulerEquations gas(config.gamma);
  FiniteVolume op(mesh,config.gamma,config.cfl,config.exact);
  const double limit=op.maxExplicitEulerStep(state);
  // Ограничение включает скорость на каждой грани, включая обе внешние.
  double maximum=0;
  GodunovFlux flux(config.gamma,config.exact);
  for (std::size_t f=0;f<=40;++f) {
    auto l=gas.toPrimitive(eulerCell(state,f==0?0:f-1));
    auto r=gas.toPrimitive(eulerCell(state,f==40?39:f));
    maximum=std::max(maximum,flux(l,r).signalSpeed);
  }
  near(limit,config.cfl*dx/maximum,"CFL bound");
  const auto before=integral(state);
  const auto leftFlux=gas.physicalFlux(gas.toPrimitive(eulerCell(state,0)));
  const auto rightFlux=gas.physicalFlux(gas.toPrimitive(eulerCell(state,39)));
  ExplicitEuler<1> euler;
  auto next=state;
  next.replaceFields(euler.advance(state,op,limit));
  const auto after=integral(next);
  const std::array<double,3> boundary{{leftFlux.mass-rightFlux.mass,leftFlux.momentum-rightFlux.momentum,leftFlux.energy-rightFlux.energy}};
  for (std::size_t c=0;c<3;++c) near(after[c]-before[c],limit*boundary[c],"boundary conservation balance");
  near(integral(state)[0],before[0],"input state changed");
  // Меняем значения того же State без изменения адреса и времени.
  auto fast=RiemannProblem1D({1,100,1},{1,100,1},config.discontinuity,config.gamma)
      .initialState(mesh,config.tStart);
  state.replaceFields(fast.fields());
  check(op.maxExplicitEulerStep(state)<limit,"stale cached CFL");
  auto rates=state.makeFields();
  op.evaluate(state,rates);
  for (const auto& field : rates)
    for (std::size_t i=0;i<40;++i) near(field[i],0,"stale cached flux");
  auto wrong=state.makeFields(); wrong.pop_back();
  fails([&]{op.evaluate(state,wrong);});
  fails([&]{euler.advance(state,op,2*op.maxExplicitEulerStep(state));});
  state.component(2)[3]=-1;
  fails([&]{op.maxExplicitEulerStep(state);});
  fails([&]{FiniteVolume bad(mesh,1.4,0.51);});
  CartesianMesh<1> other({config.xMin},{dx},{40});
  fails([&]{op.validate(RiemannProblem1D(config.left,config.right,config.discontinuity,config.gamma).initialState(other));});
  State<1> nodes(mesh,{"density","momentum","energy"},0,FieldLocation::Node);
  fails([&]{op.validate(nodes);});
}
// Точный профиль усредняем отдельно в каждой ячейке, разбив интеграл по волнам.
// Четырёхточечная квадратура Гаусса не интегрирует скачок целиком одним отрезком.
double exactAverageDensity(const ExactEulerRiemannSolver& exact, const RiemannSolution& r,
                           double a,double b,double time,double x0) {
  std::vector<double> points{a,b};
  for (double speed : {r.leftWave.headSpeed,r.leftWave.tailSpeed,r.contactSpeed,
                       r.rightWave.headSpeed,r.rightWave.tailSpeed}) {
    const double x=x0+time*speed;
    if (x>a && x<b) points.push_back(x);
  }
  std::sort(points.begin(),points.end());
  constexpr double nodes[]{-0.8611363115940526,-0.3399810435848563,0.3399810435848563,0.8611363115940526};
  constexpr double weights[]{0.3478548451374538,0.6521451548625461,0.6521451548625461,0.3478548451374538};
  double sum=0;
  for (std::size_t k=1;k<points.size();++k) {
    const double mid=(points[k-1]+points[k])/2, half=(points[k]-points[k-1])/2;
    for (int j=0;j<4;++j) sum+=half*weights[j]*exact.sample(r,(mid+half*nodes[j]-x0)/time).density;
  }
  return sum/(b-a);
}
double sodError(const EulerConfig& config, std::size_t n) {
  CartesianMesh<1> mesh({config.xMin},{(config.xMax-config.xMin)/n},{n});
  auto state=RiemannProblem1D(config.left,config.right,config.discontinuity,config.gamma)
      .initialState(mesh,config.tStart);
  ExplicitSolver<1> solver(state,std::make_unique<FiniteVolume>(mesh,config.gamma,config.cfl,config.exact),
                          std::make_unique<ExplicitEuler<1>>(),std::numeric_limits<double>::max());
  solver.run(config.tEnd);
  check(solver.state().time()==config.tEnd,"final time not exact");
  const ExactEulerRiemannSolver exact(config.gamma,config.exact);
  const auto solution=exact.solve(config.left,config.right);
  const EulerEquations gas(config.gamma);
  double error=0;
  for (std::size_t i=0;i<n;++i) {
    gas.toPrimitive(eulerCell(solver.state(),i));
    const double reference=exactAverageDensity(exact,solution,mesh.faceCoordinate(0,i),mesh.faceCoordinate(0,i+1),config.tEnd-config.tStart,config.discontinuity);
    error+=std::abs(solver.state().component(0)[i]-reference)*mesh.step()[0];
  }
  std::cout << "Sod cells=" << n << " density L1=" << error << '\n';
  return error;
}
void failureIsAtomic(const EulerConfig& config) {
  CartesianMesh<1> mesh({config.xMin},{(config.xMax-config.xMin)/config.cells},{config.cells});
  auto opts=config.exact; opts.maxIterations=1;
  auto state=RiemannProblem1D(config.left,config.right,config.discontinuity,config.gamma)
      .initialState(mesh,config.tStart);
  ExplicitSolver<1> solver(state,std::make_unique<FiniteVolume>(mesh,config.gamma,config.cfl,opts),
                          std::make_unique<ExplicitEuler<1>>(),1);
  fails([&]{solver.run(config.tEnd);});
  check(solver.state().time()==config.tStart,"failed solve advanced time");
  for (std::size_t c=0;c<3;++c)
    for (std::size_t i=0;i<config.cells;++i)
      near(solver.state().component(c)[i],state.component(c)[i],"failed solve changed state");
}
}
int main(int argc, char** argv) {
  try {
    check(argc==2,"expected Sod TOML path");
    const auto config=readEulerConfig(argv[1]);
    check(config.tEnd>config.tStart && config.cells>=2,"refinement requires positive duration and at least two cells");
    projection(); constantAndContact(); balanceAndCache(config); failureIsAtomic(config);
    const double coarse=sodError(config,config.cells/2), fine=sodError(config,config.cells);
    check(fine<0.85*coarse && fine<0.02,"Sod refinement failed");
    std::cout << "Finite volume checks passed\n";
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
