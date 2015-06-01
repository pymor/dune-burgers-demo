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
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/static_assert.hh>


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


#include "coefficients.hh"
#include "vector.hh"
#include "operator.hh"
#include "subgrid.hh"


#ifdef AS_LIB
#include <boost/python.hpp>
#endif


const int DIM = 3;


class Discretization
{
public:
  typedef Dune::SPGrid<double, DIM> G;
  typedef G::LeafGridView GV;
  static const int dim = GV::dimension;
  typedef GV::Grid::ctype Coord;

  typedef SpaceOperator<GV> OP;
  typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> Mapper;

  // Restricted operator types
  typedef Dune::ALUGrid<GV::dimension, GV::dimensionworld, Dune::cube, Dune::nonconforming> RG;
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
      int size = spt.get(std::to_string(i), -1);
      L[i] = size;
      N[i] = size * voxelPerUnit;
    }

    Cube cube(O, L);
    std::vector<Cube> cubes(1);
    cubes[0] = cube;
    int periodic = 0;
    for (int i = 0; i < DIM; i++) {
      periodic += (1 << i);
    }
    Domain domain(cubes, (Topology) periodic);

    std::cout << "Instantiating grid:  " << std::flush;
    g = std::make_shared<G>(domain, N);
    std::cout << "done" << std::endl;

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
      for (int i = 0; i < GV::dimension; i++) {
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
    std::cout << "Computing initial values:  " << std::flush;
    Vector u(mapper->size(), 0.);
    initialProjection(u);
    visualize(u, filename(0));
    std::cout << "done" << std::endl;

    double dt = T / nt;
    Vector utmp(mapper->size(), 0.);

    for(int t = 0; t < nt; t++) {
      std::cout << "\rComputing time steps:  " << std::setw(3) << t << std::setw(0) << "/" << nt << " " << std::flush;
      go->apply(u, utmp, coefficients->exponent, -dt);
      u += utmp;
      utmp = 0.;
      visualize(u, filename(t+1));
    }

    std::cout << "\rComputing time steps:  done       " << std::endl;
  }


  std::tuple< std::shared_ptr<ROP>, std::shared_ptr<std::vector<int> >, std::shared_ptr<std::vector<int> > > makeRestrictedSpaceOperator(const std::vector<int>& dofs)
  {
    auto retval = makeSubgrid<DIM, GV, Mapper>(*gv, *mapper, dofs);
    std::shared_ptr<RG> subgrid(std::get<0>(retval));
    std::shared_ptr<std::vector<int> > sourceDofs(std::get<1>(retval));
    std::shared_ptr<std::vector<int> > rangeDofs(std::get<2>(retval));

    auto rgv = std::make_shared<RGV>(subgrid->leafView());
    auto rgo = std::make_shared<ROP>(rgv, coefficients, subgrid);

    return std::make_tuple(rgo, sourceDofs, rangeDofs);
  }

  std::size_t dimSolution()
  {
    return mapper->size();
  }

#ifdef AS_LIB
  boost::python::object makeRestrictedSpaceOperatorWrapper(const boost::python::list& dofs)
  {
    using boost::python::len;
    using boost::python::extract;
    using boost::python::make_tuple;

    std::vector<int> dofVec(len(dofs));
    for (std::size_t i = 0; i < len(dofs); i++) {
      dofVec[i] = boost::python::extract<int>(dofs[i]);
    }
    auto retval = makeRestrictedSpaceOperator(dofVec);
    return boost::python::make_tuple(std::get<0>(retval), std::get<1>(retval), std::get<2>(retval));
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
  if (argc!=2) {
    std::cout << "usage: ./burgers <parameter file>" << std::endl;
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
