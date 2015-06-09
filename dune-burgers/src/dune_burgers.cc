#if HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <iomanip>
#include <tuple>

#include <dune/common/array.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/static_assert.hh>

#include <dune/grid/common/mcmgmapper.hh>

#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#pragma GCC diagnostic pop

#include <dune/grid/spgrid.hh>
#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wnew-returns-null"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
#include <dune/grid/alugrid.hh>
#pragma GCC diagnostic pop




#include <boost/get_pointer.hpp>

// the following is needed for clang
#if defined( BOOST_NO_CXX11_SMART_PTR )

namespace boost {
       template<typename T>
       T *get_pointer(std::shared_ptr<T> p)
       {
               return p.get();
       }
}

#endif

const int DIM = 3;
bool rank0;


#include "coefficients.hh"
#include "vector.hh"
#include "operator.hh"
#include "subgrid.hh"


#ifdef AS_LIB
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#endif



class Discretization
{
public:
  typedef Dune::SPGrid<double, DIM> G;
  typedef G::LeafGridView GV;
  static const int dim = GV::dimension;
  typedef GV::Grid::ctype Coord;

  typedef SpaceOperator<GV> OP;
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV, Dune::MCMGElementLayout> Mapper;

  // Restricted operator types
  typedef Dune::ALUGrid<GV::dimension, GV::dimensionworld, Dune::cube, Dune::nonconforming, Dune::No_Comm> RG;
  typedef typename RG::LeafGridView RGV;
  typedef SpaceOperator<RGV, std::shared_ptr<RG> > ROP;

  std::shared_ptr<G> g;
  std::shared_ptr<GV> gv;
  std::shared_ptr<Coefficients<DIM> > coefficients;
  std::shared_ptr<OP> go;
  std::shared_ptr<Mapper> mapper;
  double T;
  int nt;

  Discretization(const std::string paramFile)
  {
    // Instatiate grid
    typedef typename G::Traits::Domain Domain;
    typedef typename Domain::Cube Cube;
    typedef typename Domain::Topology Topology;

    Dune::ParameterTree pt;
    Dune::ParameterTreeParser::readINITree(paramFile, pt);

    Dune::FieldVector<double, DIM> O(0.);
    Dune::FieldVector<double, DIM> L(-1.);
    Dune::array<int,DIM> N;

    int voxelPerUnit = pt.get("grid.voxelPerUnit", -1);
    auto spt = pt.sub("geometry.size");
    for (int i = 0; i < DIM; i++) {
      L[i] = spt.get(std::to_string(i), -1.);
    }

    spt = pt.sub("grid.intervals");
    for (int i = 0; i < DIM; i++) {
      N[i] = spt.get(std::to_string(i), -1);
    }

    Cube cube(O, L);
    std::vector<Cube> cubes(1);
    cubes[0] = cube;
    int periodic = 0;
    for (int i = 0; i < DIM; i++) {
      periodic += (1 << i);
    }
    Domain domain(cubes, (Topology) periodic);

    Dune::array<int,DIM> overlap;
    for (int i = 0; i < DIM; i++) {
      overlap[i] = 1;
    }

    if (rank0) std::cout << "Instantiating grid:  " << std::flush;
    g = std::make_shared<G>(domain, N, overlap);
    if (rank0) std::cout << "done" << std::endl;

    gv = std::make_shared<GV>(g->leafView());
    coefficients = std::make_shared<Coefficients<DIM> >(pt);
    go = std::make_shared<OP>(gv, coefficients);
    mapper = std::make_shared<Mapper>(*gv);

    T = pt.get("timestepping.end", -1e99);
    nt = pt.get("timestepping.nt", -1);
  }

  void initialProjection(Vector& u)
  {
    auto iend = gv->template end<0>();
    for (auto it = gv->template begin<0>(); it != iend; ++it) {
      int ind = mapper->map(*it);
      auto x = it->geometry().center();
      double y = 1.;
      for (int i = 0; i < 2; i++) {
      // for (int i = 0; i < GV::dimension; i++) {
        y *= sin(2 * M_PI * x[i]);
      }
      y = 0.5 * (y + 1.);
      u[ind] = y;      
    }
  }

  void visualize(Vector& u, const std::string& filename)
  {
    Dune::VTKWriter<GV> vtkwriter(*gv);
    vtkwriter.addCellData(u, "solution");
    vtkwriter.write(filename, Dune::VTK::appendedraw);
  }

  std::string filename(int step)
  {
    std::ostringstream s;
    s << "burgers-" << std::setfill('0') << std::setw(5) << step;
    return s.str();
  }

  void solve()
  {
    if (rank0) std::cout << "Computing initial values:  " << std::flush;
    Vector u(mapper->size(), 0.);
    initialProjection(u);
    visualize(u, filename(0));
    if (rank0) std::cout << "done" << std::endl;

    double dt = T / nt;
    Vector utmp(mapper->size(), 0.);

    for(int t = 0; t < nt; t++) {
      if (rank0) std::cout << "\rComputing time steps:  " << std::setw(3) << t << std::setw(0) << "/" << nt << " " << std::flush;
      go->apply(u, utmp, coefficients->exponent, -dt);
      u += utmp;
      utmp = 0.;
      visualize(u, filename(t+1));
    }

    if (rank0) std::cout << "\rComputing time steps:  done       " << std::endl;
  }

  std::tuple<double, std::vector<double>, std::vector<int> >
    getSubgridData(const std::vector<int>& dofs, int offset)
  {
    return ::getSubgridData<DIM, GV, Mapper>(*gv, *mapper, dofs, offset);
  }


  std::tuple< std::shared_ptr<ROP>, std::vector<int>, std::vector<int> >
    makeRestrictedSpaceOperator(double diameter, const std::vector<double>& centers, const std::vector<int>& sourceDofs)
  {
    auto retval = makeSubgrid<DIM>(diameter, centers, sourceDofs);

    auto subgrid = std::get<0>(retval);
    auto rgv = std::make_shared<RGV>(subgrid->leafView());
    auto rgo = std::make_shared<ROP>(rgv, coefficients, subgrid);

    return std::make_tuple(rgo,
                           std::move(std::get<1>(retval)),
                           std::move(std::get<2>(retval)));
  }

  std::size_t dimSolution()
  {
    return mapper->size();
  }

#ifdef AS_LIB


  template<typename T> std::vector<T> toStdVector(const boost::python::object& iterable)
  {
      return std::vector<T>(boost::python::stl_input_iterator<T>(iterable),
                            boost::python::stl_input_iterator<T>());
  }

  template <class T> boost::python::list toPythonList(const std::vector<T>& vector) {
      typename std::vector<T>::iterator iter;
      boost::python::list list;
      for (const T& item: vector) {
          list.append(item);
      }
      return list;
  }

  boost::python::object getSubgridDataWrapper(const boost::python::list& dofs, int offset)
  {
    auto dofVec = toStdVector<int>(dofs);
    auto retval = getSubgridData(dofVec, offset);
    return boost::python::make_tuple(std::get<0>(retval),
                                     toPythonList<double>(std::get<1>(retval)),
                                     toPythonList<int>(std::get<2>(retval)));
  }

  boost::python::object makeRestrictedSpaceOperatorWrapper(const double diameter,
                                                           const boost::python::list& centers,
                                                           const::boost::python::list& sourceDofs)
  {
    auto centersVec = toStdVector<double>(centers);
    auto sourceDofsVec = toStdVector<int>(sourceDofs);
    auto retval = makeRestrictedSpaceOperator(diameter, centersVec, sourceDofsVec);
    return boost::python::make_tuple(std::get<0>(retval),
                                     toPythonList<int>(std::get<1>(retval)),
                                     toPythonList<int>(std::get<2>(retval)));
  }


  static void export_()
  {
    using boost::python::class_;
    using boost::python::init;
    using boost::python::make_getter;
    using boost::python::return_value_policy;
    using boost::python::return_by_value;

    class_<Discretization, boost::noncopyable>("Discretization", init<std::string>())
        .def("solve", &Discretization::solve)
        .def("getSubgridData", &Discretization::getSubgridDataWrapper)
        .def("makeRestrictedSpaceOperator", &Discretization::makeRestrictedSpaceOperatorWrapper)
        .def("visualize", &Discretization::visualize)
        .def("initialProjection", &Discretization::initialProjection)
        .def_readonly("T", &Discretization::T)
        .def_readonly("nt", &Discretization::nt)
        .add_property("go", make_getter(&Discretization::go, return_value_policy<return_by_value>()))
        .add_property("dimSolution", &Discretization::dimSolution)
    ;
  }
#endif

};


#ifdef AS_LIB

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using boost::python::class_;
using boost::python::def;
using boost::python::extract;
using boost::python::incref;
using boost::python::register_exception_translator;
using boost::python::to_python_converter;
using boost::python::vector_indexing_suite;

void exceptionTranslator(Dune::Exception const& x)
{
  PyErr_SetString(PyExc_UserWarning, x.what().c_str());
}

BOOST_PYTHON_MODULE(libdune_burgers)
{
  register_exception_translator<Dune::Exception>(exceptionTranslator);

  class_<std::vector<int>, std::shared_ptr<std::vector<int> > >("std_vector_int")
      .def(vector_indexing_suite<std::vector<int> >());

  Vector::export_();
  Discretization::OP::export_("SpaceOperator");
  Discretization::ROP::export_("RestrictedSpaceOperator");
  Discretization::export_();

}

#else

int main(int argc, char** argv)
{
  auto& helper = Dune::MPIHelper::instance(argc, argv);
  rank0 = (helper.rank() == 0);

  if (argc!=2) {
    if (rank0) std::cout << "usage: ./burgers <parameter file>" << std::endl;
    return 1;
  }

  try{
    Discretization discretization(argv[1]);
    discretization.solve();
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
  }
}

#endif
