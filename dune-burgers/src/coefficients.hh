#ifndef COEFFICIENTS_HH
#define COEFFICIENTS_HH

#include <dune/common/parametertree.hh>

template <int DIM>
struct Coefficients {

  double lambda;
  Dune::FieldVector<double, DIM> v;
  double exponent;

  Coefficients(const Dune::ParameterTree& pt)
  {
    lambda = pt.get("coefficients.lambda", -1e99);
    exponent = pt.get("coefficients.exponent", -1e99);

    auto spt = pt.sub("coefficients.v");
    for (int i = 0; i < DIM; i++) {
     v[i] = spt.get(std::to_string(i), -1e99);
    }
  }

};

#endif // COEFFICIENTS_HH

