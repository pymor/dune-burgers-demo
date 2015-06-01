#ifndef VECTOR_HH
#define VECTOR_HH

#include <dune/common/dynvector.hh>

#ifdef AS_LIB
#include <boost/python.hpp>
#endif


class Vector : public Dune::DynamicVector<double> {
public:

  Vector(int dim, double val) : Dune::DynamicVector<double>(dim, val) {}
  Vector(int dim) : Dune::DynamicVector<double>(dim) {}
  Vector(const Vector& other) : Dune::DynamicVector<double>(other) {}

  using Dune::DynamicVector<double>::operator =;

  Vector* copy()
  {
    return new Vector(*this);
  }

  void assign(Vector& other)
  {
    (*this) = other;
  }

  void setitem(std::size_t i, double value)
  {
    (*this)[i] = value;
  }

  double getitem(std::size_t i)
  {
    return (*this)[i];
  }

  double dot(const Vector& other)
  {
    return Dune::DynamicVector<double>::dot(other);
  }

  void scal(const double a)
  {
    (*this) *= a;
  }

  void axpy(const double a, const Vector& other)
  {
    Dune::DynamicVector<double>::axpy(a, other);
  }

  bool almostEqual(const Vector& other, const double rtol = 1e-14, const double atol = 1e-14) const {
    const auto sz = size();
    for (std::size_t i = 0; i < sz; i++) {
      const auto x = (*this)[i];
      const auto y = other[i];
      if (fabs(x-y) > atol + fabs(y) * rtol) {
        return false;
      }
    }
    return true;
  }

  std::size_t amaxInd() {
    std::size_t currentIndex = 0;
    double currentValue = fabs((*this)[0]);
    auto sz = size();
    for (std::size_t i = 1; i < sz; i++) {
      double newValue = fabs((*this)[i]);
      if (newValue > currentValue) {
        currentValue = newValue;
        currentIndex = i;
      }
    }
    return currentIndex;
  }

#ifdef AS_LIB
  static boost::python::object buffer(Vector& self)
  {
    using boost::python::object;
    using boost::python::handle;
    PyObject* py_buf = PyBuffer_FromReadWriteMemory(&(self[0]), self.dim() * sizeof(double));
    object retval = object(handle<>(py_buf));
    return retval;
  }

  static void export_()
  {
    using boost::python::object;
    using boost::python::class_;
    using boost::python::init;
    using boost::python::return_value_policy;
    using boost::python::manage_new_object;
    using boost::python::return_internal_reference;

    class_<Vector>("Vector", init<int, double>())
        .add_property("dim", &Vector::size)
        .def("l1Norm", &Vector::one_norm)
        .def("l2Norm", &Vector::two_norm)
        .def("supNorm", &Vector::infinity_norm)
        .def("axpy", &Vector::axpy)
        .def("__setitem__", &Vector::setitem)
        .def("__getitem__", &Vector::getitem)
        .def("dot", &Vector::dot)
        .def("copy", &Vector::copy, return_value_policy<manage_new_object>())
        .def("assign", &Vector::assign)
        .def("almostEqual", &Vector::almostEqual)
        .def("scal", &Vector::scal)
        .def("amaxInd", &Vector::amaxInd)
        .def("buffer", &buffer)
    ;
  }
#endif
};


#endif  // VECTOR_HH
