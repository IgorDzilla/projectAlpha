#include "config/EulerConfig.hpp"
#include "numerics/methods/riemann/ExactEulerRiemannSolver.hpp"
#include "numerics/core/Mesh.hpp"

#include <algorithm>
#include <limits>
#include <fstream>
#include <iostream>
#include <sstream>

namespace {
void check(bool v,const char* message) { if (!v) throw std::runtime_error(message); }
void near(double a,double b,const char* message) {
  check(std::abs(a-b)<1e-10*std::max({1.0,std::abs(a),std::abs(b)}),message);
}
}
int main(int argc,char** argv) {
  try {
    check(argc==3,"expected config and result directory");
    const auto c=readEulerConfig(argv[1]);
    const std::filesystem::path directory=argv[2];
    std::ifstream csv(directory/"solution.csv");
    check(bool(csv),"solution.csv missing");
    std::string line; std::getline(csv,line);
    check(line=="time,point,x,density,momentum,energy,velocity,pressure","unexpected header");
    CartesianMesh<1> mesh({c.xMin},{(c.xMax-c.xMin)/c.points},{c.points});
    const ExactEulerRiemannSolver exact(c.gamma,c.exact);
    const auto solution=exact.solve(c.left,c.right);
    const EulerEquations gas(c.gamma);
    std::size_t rows=0,snapshots=0;
    double time=c.tStart,previous=c.tStart;
    while (std::getline(csv,line)) {
      std::replace(line.begin(),line.end(),',',' ');
      std::istringstream in(line);
      double t,x,rho,momentum,energy,u,p;
      std::size_t cell;
      check(bool(in>>t>>cell>>x>>rho>>momentum>>energy>>u>>p),"invalid CSV row");
      std::string extra; check(!(in>>extra),"extra CSV columns");
      for (double v : {t,x,rho,momentum,energy,u,p}) check(std::isfinite(v),"nonfinite CSV value");
      check(cell==rows%c.points,"point order or layer size");
      if (cell==0) {
        previous=time; time=t;
        if (snapshots==0) check(time==c.tStart,"missing initial layer");
        else check(time>previous,"duplicate or unordered layers");
        ++snapshots;
      }
      check(t==time,"mixed times within layer");
      near(x,mesh.cellCenter(cell)[0],"wrong coordinate");
      const auto w=gas.toPrimitive({rho,momentum,energy});
      near(w.velocity,u,"velocity inconsistent with U");
      near(w.pressure,p,"pressure inconsistent with U");
      const auto expectedW=t==c.tStart ? (x<c.discontinuity ? c.left : c.right)
        : exact.sample(solution,(x-c.discontinuity)/(t-c.tStart));
      const auto expected=gas.toConservative(expectedW);
      near(rho,expected.density,"exact density");
      near(momentum,expected.momentum,"exact momentum");
      near(energy,expected.energy,"exact energy");
      near(u,expectedW.velocity,"exact velocity");
      near(p,expectedW.pressure,"exact pressure");
      ++rows;
    }
    check(snapshots==2 && rows==2*c.points && time==c.tEnd,"expected initial and final exact profiles");
    std::ifstream meta(directory/"run.txt");
    const std::string data{std::istreambuf_iterator<char>(meta),std::istreambuf_iterator<char>()};
    check(data.find("status=complete")!=std::string::npos,"run did not complete");
    check(data.find("method=exact_riemann")!=std::string::npos,"not an exact Riemann result");
    std::cout << "Verified " << snapshots << " snapshots, " << rows << " rows, final t=" << time << '\n';
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
