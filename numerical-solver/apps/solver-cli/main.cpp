#include "EulerRun.hpp"

#include <iomanip>
#include <iostream>

namespace {
void help() {
  std::cout << "Usage: solver-cli CONFIG.toml [--output-dir DIRECTORY]\n"
            << "Runs Euler 1D: finite_volume + godunov_exact + explicit_euler.\n"
            << "Config paths are relative to the TOML file; --output-dir is relative to the current directory.\n"
            << "The output directory must be empty or new. Without [output], no files are written.\n";
}
}
int main(int argc, char** argv) {
  try {
    if (argc==2 && (std::string(argv[1])=="--help" || std::string(argv[1])=="-h")) {
      help(); return 0;
    }
    if (argc!=2 && argc!=4) { help(); return 2; }
    if (argc==4 && std::string(argv[2])!="--output-dir") {
      std::cerr << "Unknown option: " << argv[2] << '\n'; return 2;
    }
    auto config=readEulerConfig(argv[1]);
    if (argc==4) {
      if (std::string(argv[3]).empty()) throw std::invalid_argument("Empty output directory");
      if (!config.output) config.output=EulerOutputConfig{};
      config.output->directory=std::filesystem::absolute(argv[3]).lexically_normal();
    }
    const auto result=runEuler(config);
    std::cout << std::setprecision(17) << "Completed: t=" << result.finalTime
              << ", cells=" << config.cells << ", saved snapshots=" << result.snapshots << '\n';
    if (config.output) std::cout << "Results: " << config.output->directory.string() << '\n';
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << '\n'; return 1;
  }
}
