#include "numerics/core/Field.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace {
void check(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

template <typename Exception, typename Action>
void rejects(Action action, const char *message) {
  try {
    action();
  } catch (const Exception &) {
    return;
  }
  throw std::runtime_error(message);
}

void checkMesh() {
  CartesianMesh<1> line({-1.0}, {0.5}, {4});
  check(line.cellCount() == 4, "1D cell count");
  check(line.cellCenter(0)[0] == -0.75 && line.cellCenter(3)[0] == 0.75,
        "1D cell centers");
  check(line.cellVolume(0) == 0.5, "1D cell length");
  check(line.origin()[0] == -1.0 && line.upperBound()[0] == 1.0,
        "Domain bounds");
  check(line.faceCoordinate(0, 0) == -1.0 && line.faceCoordinate(0, 4) == 1.0,
        "First and last faces");
  check(line.nodeCount() == 5 && line.nodeShape()[0] == 5, "1D node count");
  check(line.nodePosition(0)[0] == -1.0 && line.nodePosition(4)[0] == 1.0,
        "Boundary nodes must be included");
  check(line.nodePosition(1)[0] == -0.5, "Nodes must not be shifted to cell centers");

  CartesianMesh<2> plane({-1.0, 2.0}, {0.5, 2.0}, {3, 4});
  check(plane.cellCount() == 12 && plane.flatten({2, 1}) == 5,
        "First axis must vary fastest");
  check(plane.cellCenter(5) == CartesianMesh<2>::Point{0.25, 5.0},
        "2D center");
  check(plane.cellVolume(11) == 1.0, "2D cell area");
  for (std::size_t id = 0; id < plane.cellCount(); ++id)
    check(plane.flatten(plane.unflatten(id)) == id, "2D indexing roundtrip");
  check(plane.nodeCount() == 20 && plane.nodeShape() == CartesianMesh<2>::Index{4, 5},
        "2D nodes use the product of axis node counts");
  check(plane.flattenNode({2, 3}) == 14 &&
        plane.nodePosition(14) == CartesianMesh<2>::Point{0.0, 8.0}, "2D node coordinates");
  for (std::size_t id = 0; id < plane.nodeCount(); ++id)
    check(plane.flattenNode(plane.unflattenNode(id)) == id, "2D node indexing roundtrip");

  std::unique_ptr<Mesh<3>> box = std::make_unique<CartesianMesh<3>>(
      CartesianMesh<3>::Point{0.0, 0.0, 0.0},
      CartesianMesh<3>::Point{0.5, 2.0, 3.0},
      CartesianMesh<3>::Index{2, 3, 4});
  check(box->cellCount() == 24 && box->cellVolume(23) == 3.0,
        "Public polymorphic interface");
  check(box->nodeCount() == 60 && box->nodePosition(59) == Mesh<3>::Point{1.0, 6.0, 12.0},
        "3D nodes through the abstract mesh interface");

  rejects<std::out_of_range>([&] { line.cellCenter(4); }, "Invalid center id");
  rejects<std::out_of_range>([&] { line.cellVolume(4); }, "Invalid volume id");
  rejects<std::out_of_range>([&] { plane.flatten({3, 0}); }, "Invalid index");
  rejects<std::out_of_range>([&] { plane.unflatten(12); }, "Invalid flat id");
  rejects<std::out_of_range>([&] { line.faceCoordinate(1, 0); }, "Invalid axis");
  rejects<std::out_of_range>([&] { line.faceCoordinate(0, 5); }, "Invalid face");
  rejects<std::out_of_range>([&] { line.nodePosition(5); }, "Invalid node id");
  rejects<std::out_of_range>([&] { plane.flattenNode({4, 0}); }, "Invalid node index");
  rejects<std::out_of_range>([&] { plane.unflattenNode(20); }, "Invalid flat node id");
}

void checkInvalidGeometry() {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  for (double step : {0.0, -1.0, nan, inf})
    rejects<std::invalid_argument>([&] { CartesianMesh<1> m({0.0}, {step}, {2}); },
                                   "Invalid step accepted");
  for (double origin : {nan, inf, -inf})
    rejects<std::invalid_argument>([&] { CartesianMesh<1> m({origin}, {1.0}, {2}); },
                                   "Invalid origin accepted");
  rejects<std::invalid_argument>([] { CartesianMesh<1> m({0.0}, {1.0}, {0}); },
                                 "Zero size accepted");
  rejects<std::overflow_error>([] {
    CartesianMesh<2> m({0.0, 0.0}, {1.0, 1.0},
                      {std::numeric_limits<std::size_t>::max() / 2 + 1, 2});
  }, "Cell count overflow accepted");
  rejects<std::overflow_error>([] {
    CartesianMesh<1> m({0.0}, {1.0}, {std::numeric_limits<std::size_t>::max()});
  }, "Node axis count overflow accepted");
  rejects<std::overflow_error>([] {
    CartesianMesh<2> m({0.0, 0.0}, {1.0, 1.0},
                      {std::numeric_limits<std::size_t>::max() / 2, 1});
  }, "Node count product overflow accepted despite valid cell count");
  rejects<std::invalid_argument>([] {
    CartesianMesh<2> m({0.0, 0.0}, {1e200, 1e200}, {1, 1});
  }, "Cell volume overflow accepted");
  rejects<std::invalid_argument>([] {
    CartesianMesh<2> m({0.0, 0.0}, {1e-200, 1e-200}, {1, 1});
  }, "Cell volume underflow accepted");
  rejects<std::invalid_argument>([] {
    CartesianMesh<1> m({1e308}, {1e308}, {2});
  }, "Domain overflow accepted");
  rejects<std::invalid_argument>([] {
    CartesianMesh<1> m({1e20}, {1.0}, {2});
  }, "Collapsed domain accepted");
}

void checkField() {
  CartesianMesh<1> mesh({0.0}, {0.25}, {4});
  Field<double, 1> density(mesh, 1.0);
  check(density.size() == 4 && &density.mesh() == &mesh, "Field mesh and size");
  check(density.location() == FieldLocation::CellAverage, "Field default location");
  density[2] = 3.0;
  const auto &readOnly = density;
  check(readOnly[2] == 3.0 && readOnly.at(0) == 1.0, "Const field access");
  rejects<std::out_of_range>([&] { density[4] = 0.0; }, "Invalid field write");
  rejects<std::out_of_range>([&] { (void)readOnly[4]; }, "Invalid field read");

  auto copy = density;
  copy[2] = 7.0;
  check(density[2] == 3.0 && copy[2] == 7.0, "Field copy must own its values");
  check(copy.isCompatibleWith(density), "Copy compatibility");
  copy.fill(5.0);
  for (std::size_t i = 0; i < copy.size(); ++i)
    check(copy[i] == 5.0, "Field fill");

  Field<double, 1> centers(mesh, 0.0, FieldLocation::CellCenter);
  check(centers.size() == 4, "Cell center field size must not change");
  check(!centers.isCompatibleWith(density), "Different field locations");
  Field<double, 1> nodes(mesh, 2.0, FieldLocation::Node);
  check(nodes.size() == 5 && nodes.location() == FieldLocation::Node, "Nodal field size");
  nodes[4] = 3.0;
  check(nodes[0] == 2.0 && nodes[4] == 3.0, "Nodal boundary values");
  rejects<std::out_of_range>([&] { nodes.at(5); }, "Invalid nodal field access");
  check(!nodes.isCompatibleWith(density) && !nodes.isCompatibleWith(centers),
        "Nodal and cell fields must be incompatible");
  auto nodeCopy = nodes;
  nodeCopy.fill(7.0);
  check(nodes[4] == 3.0 && nodeCopy[4] == 7.0 && nodeCopy.isCompatibleWith(nodes),
        "Nodal field copy and fill");
  rejects<std::invalid_argument>([&] {
    Field<double, 1> invalid(mesh, 0.0, static_cast<FieldLocation>(-1));
  }, "Unknown field location accepted");
  CartesianMesh<1> otherMesh({0.0}, {0.25}, {4});
  Field<double, 1> other(otherMesh);
  check(!other.isCompatibleWith(density), "Separate meshes are not interchangeable");
  other = density;
  other[0] = 9.0;
  check(other.isCompatibleWith(density) && density[0] == 1.0,
        "Assignment copies values and mesh association");

  CartesianMesh<2> plane({0.0, 0.0}, {1.0, 1.0}, {2, 3});
  Field<int, 2> labels(plane);
  labels[plane.flatten({1, 2})] = 42;
  check(labels.at(5) == 42 && labels.at(0) == 0, "Multidimensional field indexing");
  Field<int, 2> nodeLabels(plane, 0, FieldLocation::Node);
  nodeLabels[plane.flattenNode({2, 3})] = 12;
  check(nodeLabels.size() == 12 && nodeLabels[11] == 12, "2D nodal field indexing");
}

static_assert(!std::is_constructible_v<Field<double, 1>, CartesianMesh<1> &&>,
              "A field must not accept a temporary mesh");
static_assert(!std::is_constructible_v<Field<double, 1>, const CartesianMesh<1> &&>,
              "A field must not accept a const temporary mesh");
static_assert(std::is_same_v<decltype(std::declval<const Field<double, 1> &>()[0]),
                             const double &>, "Const field access must be read-only");
} // namespace

int main() {
  try {
    checkMesh();
    checkInvalidGeometry();
    checkField();
    std::cout << "Mesh and Field tests passed\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
