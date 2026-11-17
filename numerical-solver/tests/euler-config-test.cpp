#include "EulerRun.hpp"
#include "numerics/io/CsvResultWriter.hpp"
#include "numerics/problems/RiemannProblem1D.hpp"

#include <algorithm>
#include <limits>
#include <fstream>
#include <iostream>
#include <sstream>

namespace {
void check(bool v,const char* message) { if (!v) throw std::runtime_error(message); }
template<class F> void fails(F action) {
  bool failed=false;
  try { action(); } catch (const std::exception&) { failed=true; }
  check(failed,"invalid input accepted");
}
std::string read(const std::filesystem::path& p) {
  std::ifstream in(p); check(bool(in),"file missing");
  return {std::istreambuf_iterator<char>(in),std::istreambuf_iterator<char>()};
}
}
int main(int argc,char** argv) {
  try {
    check(argc==3,"expected configs and temporary directory");
    const std::filesystem::path configs=argv[1], scratch=argv[2];
    std::filesystem::remove_all(scratch);
    std::filesystem::create_directories(scratch);
    auto fixture=toml::parse_file((configs/"toro_1.toml").string());
    auto c=parseEulerConfig(fixture,configs);
    check(c.points==200 && c.gamma==1.4 && c.tEnd==0.25,"config values");
    check(c.exact.maxIterations==100 && c.exact.pressureRtol==1e-6,"exact defaults");
    check(c.output->directory==std::filesystem::absolute(configs/"results/toro_1").lexically_normal(),"relative output base");
    check(c.output->fields.size()==5,"output values");
    const auto bad=[&](auto mutate) {
      auto t=fixture; mutate(t); fails([&]{parseEulerConfig(t,configs);});
    };
    bad([](auto& t){t.insert_or_assign("config_version",1);});
    bad([](auto& t){t.insert_or_assign("config_version",2);});
    bad([](auto& t){t.insert_or_assign("config_version",2.0);});
    bad([](auto& t){t.insert_or_assign("typo",1);});
    bad([](auto& t){t.erase("initials");});
    bad([](auto& t){t["problem"].as_table()->insert_or_assign("type","wave_1d");});
    bad([](auto& t){t["parameters"].as_table()->insert_or_assign("gamma",1);});
    bad([](auto& t){t["parameters"].as_table()->insert_or_assign("gamma",std::numeric_limits<double>::quiet_NaN());});
    bad([](auto& t){t["grid"].as_table()->insert_or_assign("points",0);});
    bad([](auto& t){t["grid"].as_table()->insert_or_assign("points",2.5);});
    bad([](auto& t){t["grid"].as_table()->insert_or_assign("dx",0.1);});
    bad([](auto& t){t["domain"].as_table()->insert_or_assign("x",toml::array{1,0});});
    bad([](auto& t){t["domain"].as_table()->insert_or_assign("t",toml::array{1,0});});
    bad([](auto& t){t["initials"].as_table()->insert_or_assign("discontinuity",1);});
    bad([](auto& t){t["initials"]["left"].as_table()->insert_or_assign("density",-1);});
    bad([](auto& t){t["initials"]["right"].as_table()->insert_or_assign("pressure",0);});
    bad([](auto& t){t["solver"].as_table()->insert_or_assign("method","finite_volume");});
    bad([](auto& t){t["solver"].as_table()->insert_or_assign("cfl",0.4);});
    bad([](auto& t){t.insert_or_assign("boundaries",toml::table{});});
    bad([](auto& t){t["solver"]["exact"].as_table()->insert_or_assign("max_iterations",0);});
    bad([](auto& t){t["solver"]["exact"].as_table()->insert_or_assign("residual_rtol",0);});
    bad([](auto& t){t["solver"].as_table()->insert_or_assign("exact",4);});
    bad([](auto& t){t["output"].as_table()->insert_or_assign("format","vtk");});
    bad([](auto& t){t["output"].as_table()->insert_or_assign("save_every",10);});
    bad([](auto& t){t["output"].as_table()->insert_or_assign("fields",toml::array{"density","density"});});
    bad([](auto& t){t["output"].as_table()->insert_or_assign("fields",toml::array{"temperature"});});
    bad([](auto& t){t["output"].as_table()->insert_or_assign("directory","");});
    bad([](auto& t){t.insert_or_assign("output",false);});
    auto minimal=fixture;
    minimal["solver"].as_table()->erase("exact");
    minimal.erase("output");
    auto noOutput=parseEulerConfig(minimal,configs);
    noOutput.tStart=0.3; noOutput.tEnd=0.3;
    const auto noFiles=runEuler(noOutput);
    check(noFiles.snapshots==0 && noFiles.finalTime==0.3,"no output / zero duration");

    c.points=20; c.tStart=0.4; c.tEnd=0.401;
    c.output->directory=scratch/"short";
    c.output->fields={"pressure","density"};
    const auto result=runEuler(c);
    check(result.snapshots==2 && result.finalTime==c.tEnd,"initial/final exact samples");
    const auto csv=read(c.output->directory/"solution.csv");
    check(csv.find("time,point,x,pressure,density\n")==0,"CSV field ordering");
    check(std::count(csv.begin(),csv.end(),'\n')==41,"CSV row count");
    check(read(c.output->directory/"run.txt").find("method=exact_riemann")!=std::string::npos,"exact metadata");
    fails([&]{runEuler(c);});
    check(read(c.output->directory/"solution.csv")==csv,"previous results overwritten");
    c.output->directory=scratch/"zero"; c.tEnd=c.tStart;
    check(runEuler(c).snapshots==1,"zero duration duplicated output");
    c.output->directory=scratch/"invalid_time"; c.tEnd=c.tStart-0.1;
    fails([&]{runEuler(c);});
    check(!std::filesystem::exists(c.output->directory),"invalid time created output");
    c.tEnd=c.tStart+0.1;
    c.output->directory=scratch/"vacuum"; c.left={1,-10,1}; c.right={1,10,1};
    fails([&]{runEuler(c);});
    check(!std::filesystem::exists(c.output->directory),"vacuum created output");
    std::cout << "Euler config and output checks passed\n";
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
