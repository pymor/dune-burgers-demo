#ifndef SUBGRID_HH
#define SUBGRID_HH

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/alugrid.hh>

template <int DIM, class GV, class MAPPER>
std::tuple<Dune::ALUGrid<DIM, DIM, Dune::cube, Dune::nonconforming> *, std::vector<int> *, std::vector<int> *>
  makeSubgrid(const GV& gv, const MAPPER& mapper, const std::vector<int>& dofs)
{
  typedef Dune::ALUGrid<DIM, DIM, Dune::cube, Dune::nonconforming> RG;
  typedef Dune::FieldVector<double, DIM> FV;

  // Accumulate data
  const double diameter = pow(gv.template begin<0>()->geometry().volume(), 1./DIM);
  const auto geometryType = gv.template begin<0>()->geometry().type();

  const int stencil = 2 * DIM + 1;
  std::vector<int> sourceDofs(dofs.size() * stencil, 0);
  std::vector<FV> centers(dofs.size());
  int dofCounter = 0;

  auto iend = gv.template end<0>();
  for (auto it = gv.template begin<0>(); it != iend; ++it) {
    const int index = mapper.map(*it);
    auto pos = std::find(dofs.begin(), dofs.end(), index);
    if (pos != dofs.end()) {
      auto sourceDofIndex = pos - dofs.begin();
      sourceDofs[sourceDofIndex * stencil] = index;
      centers[sourceDofIndex] = it->geometry().center();

      auto isend = gv.iend(*it);
      for (auto is = gv.ibegin(*it); is != isend; ++is) {
        const int localIndex = is->indexInInside();
        const int outsideIndex = mapper.map(*(is->outside()));
        sourceDofs[sourceDofIndex * stencil + localIndex + 1] = outsideIndex;
      }
      dofCounter++;
    }
  }

  // make sure we have found all dofs
  assert (dofCounter == dofs.size());

  Dune::GridFactory<RG> factory;

  int vertexCounter = 0;
  for (const auto& center: centers) {
    auto origin = center + FV(-diameter * 1.5);

    // insert all vertices of stencil
    std::vector<int> vertexNumbers((1 << (2 * DIM)), -1);
    for (int vertexNumber = 0; vertexNumber < (1 << (2 * DIM)); vertexNumber++) {
      int vn = vertexNumber;
      bool corner = true;
      for (int d = 0; d < DIM; d++) {
        if ((vn % 4 == 1) || (vn % 4 == 2)) {
          corner = false;
          break;
        }
        vn /= 4;
      }
      if (corner) {
        continue;
      }
      vertexNumbers[vertexNumber] = vertexCounter;
      vertexCounter++;

      FV vertex(origin);
      vn = vertexNumber;
      for (int dim = 0; dim < DIM; dim++) {
          vertex[dim] += diameter * (vn % 4);
          vn /= 4;
      }
      factory.insertVertex(vertex);
    }

    // insert center element
    std::vector<unsigned int> vertexIds((1 << DIM));
    int centerllc = 0;
    for (int d = 0; d < DIM; d++) {
      centerllc += 1 << (2 * d);
    }
    for (int vertexNumber = 0; vertexNumber < (1 << DIM); vertexNumber++) {
      int pos = centerllc;
      int vn = vertexNumber;
      for (int d = 0; d < DIM; d++) {
        pos += (1 << (2 * d)) * (vn % 2);
        vn /= 2;
      }
      vertexIds[vertexNumber] = vertexNumbers[pos];
    }
    factory.insertElement(geometryType, vertexIds);

    // insert neighbour elements
    for (int dim = 0; dim < DIM; dim++) {
      for (int side = 0; side < 2; side++) {
        int llc = centerllc;
        llc += (1 << (2 * dim)) * (side * 2 - 1);

        for (int vertexNumber = 0; vertexNumber < (1 << DIM); vertexNumber++) {
          int pos = llc;
          int vn = vertexNumber;
          for (int d = 0; d < DIM; d++) {
            pos += (1 << (2 * d)) * (vn % 2);
            vn /= 2;
          }
          vertexIds[vertexNumber] = vertexNumbers[pos];
        }
        factory.insertElement(geometryType, vertexIds);
      }
    }
  }

  auto subgrid = factory.createGrid();
  auto sgv = subgrid->leafView();
  Dune::SingleCodimSingleGeomTypeMapper<typename RG::LeafGridView, 0> smapper(sgv);
  auto reorderedSourceDofs = new std::vector<int>(sourceDofs.size());
  auto rangeDofs = new std::vector<int>(dofs.size());

  {
    auto iend = sgv.template end<0>();
    for (auto it = sgv.template begin<0>(); it != iend; ++it) {
     auto index = smapper.map(*it);
     auto insertionIndex = factory.insertionIndex(*it);
     (*reorderedSourceDofs)[index] = sourceDofs[insertionIndex];
     if (insertionIndex % stencil == 0) {
       (*rangeDofs)[insertionIndex / stencil] = index;
     }
    }
  }
  return std::make_tuple(subgrid, reorderedSourceDofs, rangeDofs);
}




#endif // SUBGRID_HH
