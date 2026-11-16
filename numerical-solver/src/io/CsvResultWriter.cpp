#include "numerics/io/CsvResultWriter.hpp"

#include <iomanip>
#include <locale>
#include <set>

CsvResultWriter::CsvResultWriter(const std::filesystem::path& directory,
                                 double gamma, std::vector<std::string> fields)
    : gas_(gamma), fields_(std::move(fields)) {
  const std::set<std::string> allowed{"density","momentum","energy","velocity","pressure"};
  std::set<std::string> seen;
  if (fields_.empty()) throw std::invalid_argument("CSV requires at least one field");
  for (const auto& name : fields_)
    if (!allowed.count(name) || !seen.insert(name).second)
      throw std::invalid_argument("Unknown or duplicate CSV field: " + name);
  if (directory.empty()) throw std::invalid_argument("CSV directory is empty");
  // Не затираем результаты предыдущего запуска.
  if (std::filesystem::exists(directory) && !std::filesystem::is_empty(directory))
    throw std::runtime_error("Output directory is not empty: " + directory.string());
  std::filesystem::create_directories(directory);
  stream_.exceptions(std::ios::failbit | std::ios::badbit);
  stream_.open(directory / "solution.csv");
  stream_.imbue(std::locale::classic());
  stream_ << std::setprecision(17) << "time,cell,x";
  for (const auto& name : fields_) stream_ << ',' << name;
  stream_ << '\n';
}

void CsvResultWriter::write(const State<1>& state) {
  validateEulerLayout(state);
  // Вначале проверяем весь слой, чтобы ошибка физики не записала его половину.
  std::vector<EulerPrimitiveState> primitive;
  primitive.reserve(state.mesh().cellCount());
  for (std::size_t i=0; i<state.mesh().cellCount(); ++i)
    primitive.push_back(gas_.toPrimitive(eulerCell(state,i)));
  for (std::size_t i=0; i<state.mesh().cellCount(); ++i) {
    const auto u=eulerCell(state,i);
    const auto& w=primitive[i];
    stream_ << state.time() << ',' << i << ',' << state.mesh().cellCenter(i)[0];
    for (const auto& name : fields_) {
      double value=0;
      if (name=="density") value=u.density;
      else if (name=="momentum") value=u.momentum;
      else if (name=="energy") value=u.energy;
      else if (name=="velocity") value=w.velocity;
      else if (name=="pressure") value=w.pressure;
      stream_ << ',' << value;
    }
    stream_ << '\n';
  }
  stream_.flush();
  ++snapshots_;
}
