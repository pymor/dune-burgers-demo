#ifndef SUBGRID_HH
#define SUBGRID_HH

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/alugrid.hh>
#include <dune/grid/common/scsgmapper.hh>

template <class GV, class M>
class NeighbourExchange
    : public Dune::CommDataHandleIF<NeighbourExchange<GV, M>, int>
{
public:
  typedef int DataType;
  typedef std::map<int, Dune::FieldVector<int, DIM * 2> > MapType;

  NeighbourExchange (const GV& gv, const M & mapper, MapType& map, int offset)
    : gv_(gv), mapper_(mapper), map_(map), offset_(offset)
  {}

  bool contains(int dim, int codim) const
  {
    return (codim == 0);
  }

  bool fixedsize (int dim, int codim) const
  {
    return true;
  }

  template < class EntityType >
  std::size_t size (EntityType & e) const
  {
    return DIM * 2;
  }

  template <class MessageBuffer, class EntityType, int CODIM>
  class gather_
  {
  public:
    static void doIt(MessageBuffer& buff, const EntityType& e,
                     const GV & gv, const M & mapper, double offset)
    { }
  };

  template <class MessageBuffer, class EntityType>
  class gather_<MessageBuffer, EntityType, 0>
  {
  public:
    static void doIt(MessageBuffer& buff, const EntityType& e,
                     const GV & gv, const M & mapper, double offset)
    {
      MapType::mapped_type neighbours;
      auto isend = gv.iend(e);
      for (auto is = gv.ibegin(e); is != isend; ++is) {
        const int localIndex = is->indexInInside();
        const int outsideIndex = mapper.map(*(is->outside()));
        neighbours[localIndex] = outsideIndex + offset;
      }
      for (int i = 0; i < DIM * 2; i++) {
        buff.write(neighbours[i]);
      }
    }
  };

  template <class MessageBuffer, class EntityType >
  void gather(MessageBuffer& buff, const EntityType& e ) const
  {
    gather_<MessageBuffer, EntityType, EntityType::codimension>::doIt(buff, e, gv_, mapper_, offset_);
  }

  template <class MessageBuffer, class EntityType >
  void scatter(MessageBuffer& buff, const EntityType& e , std::size_t n)
  {
    MapType::mapped_type neighbours;
    for (int i = 0; i < DIM * 2; i++) {
      buff.read(neighbours[i]);
    }
    map_[mapper_.map(e)] = neighbours;
  }

private :
  const GV & gv_;
  const M & mapper_;
  MapType & map_;
  int offset_;
};


template <int DIM, class GV, class MAPPER>
std::tuple<double, std::vector<double>, std::vector<int> >
  getSubgridData(const GV& gv, const MAPPER& mapper, const std::vector<int>& dofs, int offset)
{

  typedef Dune::FieldVector<double, DIM> FV;

  // Accumulate data
  const double diameter = pow(gv.template begin<0>()->geometry().volume(), 1./DIM);

  typedef NeighbourExchange<GV, MAPPER> NE;
  typename NE::MapType overlapMap;
  NE ne(gv, mapper, overlapMap, offset);
  gv.template communicate<NE>(ne, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);

  const int stencil = 2 * DIM + 1;
  std::vector<double> centers(dofs.size() * DIM, 0.);
  std::vector<int> sourceDofs(dofs.size() * stencil, 0);

  auto iend = gv.template end<0>();
  for (auto it = gv.template begin<0>(); it != iend; ++it) {
    const int index = mapper.map(*it);
    auto pos = std::find(dofs.begin(), dofs.end(), index + offset);
    if (pos != dofs.end()) {
      auto sourceDofIndex = pos - dofs.begin();
      sourceDofs[sourceDofIndex * stencil] = index + offset;

      if (it->partitionType() == Dune::InteriorEntity) {
        for (int d = 0; d < DIM; d++) {
          auto center = it->geometry().center();
          centers[sourceDofIndex * DIM + d] = center[d];
        }

        auto isend = gv.iend(*it);
        for (auto is = gv.ibegin(*it); is != isend; ++is) {
          const int localIndex = is->indexInInside();
          const int outsideIndex = mapper.map(*(is->outside()));
          sourceDofs[sourceDofIndex * stencil + localIndex + 1] = outsideIndex + offset;
        }
      } else {
        Dune::FieldVector<int, DIM * 2> neighbours = overlapMap[index];
        for (int i = 0; i < DIM * 2; i++) {
          sourceDofs[sourceDofIndex * stencil + i + 1] = neighbours[i];
        }
      }

    }
  }

  return std::make_tuple(diameter, std::move(centers), std::move(sourceDofs));
}



template <int DIM>
std::tuple<std::shared_ptr<Dune::ALUGrid<DIM, DIM, Dune::cube, Dune::nonconforming, Dune::No_Comm> >, std::vector<int>, std::vector<int> >
  makeSubgrid(double diameter, const std::vector<double>& centers, const std::vector<int>& sourceDofs)
{
  Dune::GeometryType geometryType;
  geometryType.makeCube(DIM);
  const int stencil = 2 * DIM + 1;
  typedef Dune::ALUGrid<DIM, DIM, Dune::cube, Dune::nonconforming, Dune::No_Comm> RG;
  typedef Dune::FieldVector<double, DIM> FV;

  const int sourceDofCount = centers.size() / DIM;

  Dune::GridFactory<RG> factory;

  int vertexCounter = 0;
  for (int sourceDofIndex = 0; sourceDofIndex < sourceDofCount; sourceDofIndex++) {
    FV center;
    for (int d = 0; d < DIM; d++) {
      center[d] = centers[sourceDofIndex * DIM + d];
    }
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

  std::shared_ptr<RG> subgrid(factory.createGrid());

  auto sgv = subgrid->leafView();
  Dune::SingleCodimSingleGeomTypeMapper<typename RG::LeafGridView, 0> smapper(sgv);
  std::vector<int> reorderedSourceDofs(sourceDofs.size());
  std::vector<int> rangeDofs(sourceDofCount);

  {
    auto iend = sgv.template end<0>();
    for (auto it = sgv.template begin<0>(); it != iend; ++it) {
     auto index = smapper.map(*it);
     auto insertionIndex = factory.insertionIndex(*it);
     reorderedSourceDofs[index] = sourceDofs[insertionIndex];
     if (insertionIndex % stencil == 0) {
       rangeDofs[insertionIndex / stencil] = index;
     }
    }
  }
  return std::make_tuple(subgrid, std::move(reorderedSourceDofs), std::move(rangeDofs));
}

#endif // SUBGRID_HH
