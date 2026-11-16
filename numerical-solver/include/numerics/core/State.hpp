#pragma once

#include "Field.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Поля одной задачи и время, которому они соответствуют.
template <std::size_t Dim> class State {
public:
  using Fields = std::vector<Field<double, Dim>>;

private:
  const Mesh<Dim> *mesh_;
  std::vector<std::string> names_;
  FieldLocation location_;
  Fields fields_;
  double time_;

public:
  State(const Mesh<Dim> &mesh, std::vector<std::string> names,
        double time = 0.0, FieldLocation location = FieldLocation::CellAverage)
      : mesh_(&mesh), names_(std::move(names)), location_(location), time_(time) {
    if (names_.empty())
      throw std::invalid_argument("State must have at least one component");
    for (std::size_t i = 0; i < names_.size(); ++i) {
      if (names_[i].empty())
        throw std::invalid_argument("State component name must not be empty");
      for (std::size_t j = 0; j < i; ++j)
        if (names_[i] == names_[j])
          throw std::invalid_argument("State component names must be unique");
    }
    setTime(time);
    fields_ = makeFields();
    validate();
  }

  State(const Mesh<Dim> &&, std::vector<std::string>, double = 0.0,
        FieldLocation = FieldLocation::CellAverage) = delete;

  const Mesh<Dim> &mesh() const noexcept { return *mesh_; }
  std::size_t componentCount() const noexcept { return names_.size(); }
  const std::vector<std::string> &names() const noexcept { return names_; }
  FieldLocation location() const noexcept { return location_; }
  double time() const noexcept { return time_; }

  void setTime(double time) {
    if (!std::isfinite(time))
      throw std::invalid_argument("State time must be finite");
    time_ = time;
  }

  Field<double, Dim> &component(std::size_t id) { return fields_.at(id); }
  const Field<double, Dim> &component(std::size_t id) const { return fields_.at(id); }
  const Fields &fields() const noexcept { return fields_; }

  // Такие же поля без времени: подходят для правой части и рабочих буферов.
  Fields makeFields(double value = 0.0) const {
    if (!std::isfinite(value))
      throw std::invalid_argument("Initial field value must be finite");
    Fields result;
    result.reserve(componentCount());
    for (std::size_t i = 0; i < componentCount(); ++i)
      result.emplace_back(*mesh_, value, location_);
    return result;
  }

  void validateFields(const Fields &fields) const {
    if (fields.size() != componentCount())
      throw std::invalid_argument("State component count mismatch");
    // Все компоненты состояния имеют одно размещение: узлы либо ячейки.
    const std::size_t expectedSize = fieldSize(*mesh_, location_);
    for (const auto &field : fields) {
      if (&field.mesh() != mesh_ || field.location() != location_ ||
          field.size() != expectedSize)
        throw std::invalid_argument("State field mesh, size or location mismatch");
      for (std::size_t i = 0; i < field.size(); ++i)
        if (!std::isfinite(field[i]))
          throw std::runtime_error("Non-finite field value at index " + std::to_string(i));
    }
  }

  void validate() const {
    if (!std::isfinite(time_) || mesh_->cellCount() == 0 || names_.empty())
      throw std::invalid_argument("Invalid state time, empty mesh or missing components");
    validateFields(fields_);
  }

  void replaceFields(Fields fields) {
    validateFields(fields);
    fields_ = std::move(fields);
  }
};
