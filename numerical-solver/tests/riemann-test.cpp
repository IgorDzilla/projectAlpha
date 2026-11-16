#include "numerics/methods/riemann/ExactEulerRiemannSolver.hpp"
#include "numerics/core/State.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
void check(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
void near(double a, double b, const char* message, double tolerance = 2.0e-9) {
  check(std::isfinite(a) && std::isfinite(b) &&
        std::abs(a-b) <= tolerance * std::max({1.0, std::abs(a), std::abs(b)}), message);
}
void same(const EulerPrimitiveState& a, const EulerPrimitiveState& b) {
  near(a.density,b.density,"density mismatch");
  near(a.velocity,b.velocity,"velocity mismatch");
  near(a.pressure,b.pressure,"pressure mismatch");
}
template<class Error, class F> void throws(F action) {
  bool caught = false;
  try { action(); } catch (const Error&) { caught = true; }
  check(caught,"expected exception was not thrown");
}
void checkShock(const EulerEquations& gas, const EulerPrimitiveState& a,
                const EulerPrimitiveState& b, double speed) {
  const auto ua=gas.toConservative(a), ub=gas.toConservative(b);
  const auto fa=gas.physicalFlux(a), fb=gas.physicalFlux(b);
  near(fb.mass-fa.mass,speed*(ub.density-ua.density),"mass jump",2e-8);
  near(fb.momentum-fa.momentum,speed*(ub.momentum-ua.momentum),"momentum jump",2e-8);
  near(fb.energy-fa.energy,speed*(ub.energy-ua.energy),"energy jump",2e-8);
}
void checkGasAndState() {
  for (double gamma : {1.4,5.0/3.0}) {
    const EulerEquations gas(gamma);
    const EulerPrimitiveState w{2.0,3.0,4.0};
    const auto u=gas.toConservative(w);
    near(u.density,2.0,"rho"); near(u.momentum,6.0,"momentum");
    near(u.energy,4.0/(gamma-1.0)+9.0,"energy per volume");
    same(gas.toPrimitive(u),w);
    const auto f=gas.physicalFlux(w);
    near(f.mass,6.0,"mass flux"); near(f.momentum,22.0,"momentum flux");
    near(f.energy,3.0*(u.energy+4.0),"energy flux");
    near(gas.soundSpeed(w),std::sqrt(gamma*2.0),"sound speed");
  }
  // Данные из общих ячеечных полей можно передать в решатель без другой сетки.
  CartesianMesh<1> mesh({0.0},{0.5},{2});
  State<1> state(mesh,{"density","momentum","energy"});
  const EulerEquations gas;
  const auto left=gas.toConservative({1.0,0.0,1.0});
  const auto right=gas.toConservative({0.125,0.0,0.1});
  state.component(0)[0]=left.density; state.component(1)[0]=left.momentum;
  state.component(2)[0]=left.energy;
  state.component(0)[1]=right.density; state.component(1)[1]=right.momentum;
  state.component(2)[1]=right.energy;
  const auto at=[&](std::size_t i) {
    return gas.toPrimitive({state.component(0)[i],state.component(1)[i],state.component(2)[i]});
  };
  const ExactEulerRiemannSolver solver(gas.gamma());
  const auto r=solver.solve(at(0),at(1));
  near(r.pStar,0.3031301780506468,"Sod reference pressure");
  near(r.uStar,0.9274526200489499,"Sod reference velocity");
  near(state.component(0)[0],1.0,"local solver must not change fields");
  check(state.time()==0.0,"local solver must not change time");
}
void checkWaves(double gamma) {
  const EulerEquations gas(gamma);
  const ExactEulerRiemannSolver solver(gamma);
  const EulerPrimitiveState left{1,0,1}, right{0.125,0,0.1};
  const auto r=solver.solve(left,right);
  check(r.leftWave.type==RiemannWaveType::Rarefaction,"left fan");
  check(r.rightWave.type==RiemannWaveType::Shock,"right shock");
  check(r.exactPressureRoot || (r.relativeResidual<=1e-10 && r.relativeBracketWidth<=1e-10),
        "both convergence criteria");
  checkShock(gas,r.starRight,right,r.rightWave.headSpeed);
  same(solver.sample(r,-100.0),left); same(solver.sample(r,100.0),right);
  same(solver.sample(r,r.contactSpeed),r.starRight);
  same(solver.sample(r,r.rightWave.headSpeed),right);
  same(solver.sample(r,r.leftWave.headSpeed),left);
  same(solver.sample(r,r.leftWave.tailSpeed),r.starLeft);
  const double ray=0.5*(r.leftWave.headSpeed+r.leftWave.tailSpeed);
  const auto fan=solver.sample(r,ray);
  near(fan.pressure/std::pow(fan.density,gamma),left.pressure/std::pow(left.density,gamma),
       "left fan entropy");
  near(fan.velocity+2*gas.soundSpeed(fan)/(gamma-1),
       left.velocity+2*gas.soundSpeed(left)/(gamma-1),"left fan invariant");
  near(fan.velocity-gas.soundSpeed(fan),ray,"left fan characteristic");
  // Зеркальная задача покрывает левый удар и правый веер.
  const auto mirrored=solver.solve(right,left);
  near(mirrored.pStar,r.pStar,"mirror pressure");
  near(mirrored.uStar,-r.uStar,"mirror velocity");
  checkShock(gas,right,mirrored.starLeft,mirrored.leftWave.headSpeed);
  same(solver.sample(mirrored,mirrored.leftWave.headSpeed),mirrored.starLeft);
  const auto reflected=solver.sample(mirrored,-ray);
  same(reflected,{fan.density,-fan.velocity,fan.pressure});
  same(solver.sample(mirrored,mirrored.rightWave.tailSpeed),mirrored.starRight);
  same(solver.sample(mirrored,mirrored.rightWave.headSpeed),left);
  // Сдвиг всех скоростей помещает неподвижную грань внутрь левого веера.
  const double shift=-ray;
  const auto transonic=solver.solve({1,shift,1},{0.125,shift,0.1});
  const auto face=solver.sample(transonic,0);
  same(face,{fan.density,fan.velocity+shift,fan.pressure});
  near(face.velocity,gas.soundSpeed(face),"sonic face inside fan");
  const auto negative=solver.solve({0.125,-shift,0.1},{1,-shift,1});
  same(solver.sample(negative,0),{face.density,-face.velocity,face.pressure});
  for (double speed : {r.leftWave.headSpeed,r.leftWave.tailSpeed,
                       r.rightWave.headSpeed,r.contactSpeed,
                       gas.soundSpeed(left),gas.soundSpeed(right)}) {
    check(r.signalSpeed>=std::abs(speed),"missing signal speed");
  }
  for (double velocity : {-100.0,100.0}) {
    const auto advected=solver.solve({1,velocity,1},{0.125,velocity,0.1});
    same(solver.sample(advected,0),velocity>0 ? EulerPrimitiveState{1,velocity,1}
                                                         : EulerPrimitiveState{0.125,velocity,0.1});
  }
}
void checkContactsAndExtremes() {
  const ExactEulerRiemannSolver solver;
  const EulerEquations gas;
  for (double velocity : {-2.0,0.0,2.0}) {
    const EulerPrimitiveState l{1,velocity,1}, r{3,velocity,1};
    const auto contact=solver.solve(l,r);
    check(contact.iterations==0,"contact should need no iterations");
    same(solver.sample(contact,velocity-0.1),l);
    same(solver.sample(contact,velocity),r);
    same(solver.sample(contact,velocity+0.1),r);
    const auto uniform=solver.solve(l,l);
    for (double ray : {-100.0,velocity,100.0}) same(solver.sample(uniform,ray),l);
    check(uniform.signalSpeed>=std::abs(velocity)+gas.soundSpeed(l),"zero wave CFL");
  }
  const auto stationary=solver.solve({1,0,2},{3,0,2});
  const auto flux=gas.physicalFlux(solver.sample(stationary,0));
  near(flux.mass,0,"stationary contact mass"); near(flux.momentum,2,"stationary contact momentum");
  near(flux.energy,0,"stationary contact energy");
  const auto collision=solver.solve({1,2,1},{1,-2,1});
  check(collision.leftWave.type==RiemannWaveType::Shock &&
        collision.rightWave.type==RiemannWaveType::Shock,"two shocks");
  near(collision.uStar,0,"symmetric collision");
  checkShock(gas,collision.left,collision.starLeft,collision.leftWave.headSpeed);
  checkShock(gas,collision.starRight,collision.right,collision.rightWave.headSpeed);
  const double a=gas.soundSpeed({1,0,1});
  const double velocity=0.99*2*a/(gas.gamma()-1);
  const auto expansion=solver.solve({1,-velocity,1},{1,velocity,1});
  const double expected=std::pow(1-(gas.gamma()-1)*velocity/(2*a),2*gas.gamma()/(gas.gamma()-1));
  check(expansion.pStar>0,"near vacuum pressure positive");
  check(std::abs(expansion.pStar/expected-1)<2e-8,"analytic two-rarefaction pressure");
  // Изменение единиц плотности и давления не должно менять скорости волн.
  const auto reference=solver.solve({1,0,1},{0.125,0,0.1});
  for (double scale : {1e-100,1e100}) {
    const auto scaled=solver.solve({scale,0,scale},{0.125*scale,0,0.1*scale});
    near(scaled.pStar/scale,reference.pStar,"scaled pressure");
    near(scaled.starLeft.density/scale,reference.starLeft.density,"scaled density");
    near(scaled.uStar,reference.uStar,"scaled velocity");
    near(scaled.signalSpeed,reference.signalSpeed,"scaled signal speed");
  }
  const auto strong=solver.solve({1,0,1000},{1,0,0.01});
  check(strong.signalSpeed>=std::abs(strong.rightWave.headSpeed),"strong shock speed missing");
  // Стационарный нормальный удар: upstream Mach=2, gamma=1.4.
  const double u=2*std::sqrt(1.4), rhoRatio=8.0/3.0;
  const EulerPrimitiveState upstream{1,u,1}, downstream{rhoRatio,u/rhoRatio,4.5};
  const auto shock=solver.solve(upstream,downstream);
  const auto f1=gas.physicalFlux(upstream), f2=gas.physicalFlux(downstream);
  near(f1.mass,f2.mass,"stationary shock mass");
  near(f1.momentum,f2.momentum,"stationary shock momentum");
  near(f1.energy,f2.energy,"stationary shock energy");
  const auto sampled=gas.physicalFlux(solver.sample(shock,0));
  near(sampled.energy,f1.energy,"stationary shock sample flux");
}
void checkFailures() {
  const double nan=std::numeric_limits<double>::quiet_NaN();
  const double inf=std::numeric_limits<double>::infinity();
  const ExactEulerRiemannSolver solver;
  const EulerEquations gas;
  const EulerPrimitiveState w{1,0,1};
  for (double gamma : {0.0,1.0,nan,inf}) throws<std::invalid_argument>([&]{ExactEulerRiemannSolver bad(gamma);});
  for (auto bad : {EulerPrimitiveState{0,0,1},{-1,0,1},{1,0,0},{1,nan,1},{1,0,inf}}) {
    throws<std::invalid_argument>([&]{solver.solve(bad,w);});
    throws<std::invalid_argument>([&]{solver.solve(w,bad);});
  }
  throws<std::invalid_argument>([&]{gas.toPrimitive({1,4,1});});
  throws<std::runtime_error>([&]{gas.physicalFlux({1,1e200,1});});
  for (double velocity : {2*gas.soundSpeed(w)/(gas.gamma()-1),10.0}) {
    throws<std::domain_error>([&]{solver.solve({1,-velocity,1},{1,velocity,1});});
  }
  auto r=solver.solve(w,{0.125,0,0.1});
  throws<std::invalid_argument>([&]{solver.sample(r,nan);});
  throws<std::invalid_argument>([&]{solver.sample(r,inf);});
  throws<std::invalid_argument>([&]{ExactEulerRiemannSolver(5.0/3.0).sample(r,0);});
  throws<std::invalid_argument>([&]{solver.sample(RiemannSolution{},0);});
  r.gamma=nan;
  throws<std::invalid_argument>([&]{solver.sample(r,0);});
  ExactRiemannOptions options;
  options.maxIterations=1;
  throws<std::runtime_error>([&]{ExactEulerRiemannSolver(1.4,options).solve(w,{0.125,0,0.1});});
  options={}; options.maxBracketExpansions=1;
  throws<std::runtime_error>([&]{ExactEulerRiemannSolver(1.4,options).solve({1,100,1},{1,-100,1});});
  options={}; options.pressureRtol=nan;
  throws<std::invalid_argument>([&]{ExactEulerRiemannSolver bad(1.4,options);});
  options={}; options.residualRtol=0;
  throws<std::invalid_argument>([&]{ExactEulerRiemannSolver bad(1.4,options);});
  options={}; options.maxIterations=0;
  throws<std::invalid_argument>([&]{ExactEulerRiemannSolver bad(1.4,options);});
  options={}; options.maxBracketExpansions=0;
  throws<std::invalid_argument>([&]{ExactEulerRiemannSolver bad(1.4,options);});
}
}
int main() {
  try {
    checkGasAndState(); checkWaves(1.4); checkWaves(5.0/3.0);
    checkContactsAndExtremes(); checkFailures();
    std::cout<<"Riemann solver and Euler model checks passed\n";
  } catch (const std::exception& error) {
    std::cerr<<error.what()<<'\n'; return 1;
  }
}
