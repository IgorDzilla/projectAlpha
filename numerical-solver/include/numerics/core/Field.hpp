#pragma once

#include "Mesh.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

// Что означает каждое число в поле. Смена типа хранения сама по себе
// не выполняет интерполяцию или усреднение: это задача численного метода.
enum class FieldLocation {
  CellCenter,  // Значение функции в центре ячейки; одно число на ячейку.
  CellAverage, // Среднее функции по объёму ячейки; одно число на ячейку.
  Node         // Значение функции в узле (вершине), включая узлы на границе.
};

// В 1D: N значений для CellCenter/CellAverage, N+1 для Node.
// В нескольких измерениях число узлов равно произведению (shape[d] + 1).
template <std::size_t Dim>
std::size_t fieldSize(const Mesh<Dim> &mesh, FieldLocation location) {
  switch (location) {
  case FieldLocation::CellCenter:
  case FieldLocation::CellAverage:
    return mesh.cellCount();
  case FieldLocation::Node:
    return mesh.nodeCount();
  }
  throw std::invalid_argument("Unknown field location");
}

template <typename T, std::size_t Dim> class Field {
public:
  using Container = std::vector<T>;
  using Reference = typename Container::reference;
  using ConstReference = typename Container::const_reference;

private:
  // Сетка должна жить дольше поля; её нельзя менять или перемещать при живом поле.
  const Mesh<Dim> *mesh_;
  FieldLocation location_;
  Container values_;

public:
  // По умолчанию храним средние по ячейкам, как в методе конечных объёмов.
  explicit Field(const Mesh<Dim> &mesh, const T &initialValue = T{},
                 FieldLocation location = FieldLocation::CellAverage)
      : mesh_(&mesh), location_(location), values_(fieldSize(mesh, location), initialValue) {}

  // Не позволяем привязать поле к временной сетке.
  Field(const Mesh<Dim> &&, const T & = T{},
        FieldLocation = FieldLocation::CellAverage) = delete;

  std::size_t size() const noexcept { return values_.size(); }

  const Mesh<Dim> &mesh() const noexcept { return *mesh_; }

  FieldLocation location() const noexcept { return location_; }

  // Для Node id — номер узла, для остальных вариантов — номер ячейки.
  // Оба способа доступа проверяют индекс.
  Reference at(std::size_t id) { return values_.at(id); }
  ConstReference at(std::size_t id) const { return values_.at(id); }

  Reference operator[](std::size_t id) { return at(id); }
  ConstReference operator[](std::size_t id) const { return at(id); }

  void fill(const T &value) { std::fill(values_.begin(), values_.end(), value); }

  // Одинаковых размеров недостаточно: сетка и смысл значений должны совпадать.
  bool isCompatibleWith(const Field &other) const noexcept {
    return mesh_ == other.mesh_ && location_ == other.location_ &&
           size() == other.size();
  }
};
