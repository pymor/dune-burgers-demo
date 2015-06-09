#ifndef OPERATOR_HH
#define OPERATOR_HH


#include<dune/grid/common/scsgmapper.hh>

#include "coefficients.hh"
#include "vector.hh"

#ifdef AS_LIB
#include <boost/python.hpp>
#endif


template <class GV, typename D = char>
class SpaceOperator
{
public:
  typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> Mapper;

  SpaceOperator(std::shared_ptr<GV> gv, std::shared_ptr<Coefficients<GV::dimension> >& coeffs)
      : gv_(gv), coefficients_(coeffs), mapper_(*gv)
  {}

  SpaceOperator(std::shared_ptr<GV> gv, std::shared_ptr<Coefficients<GV::dimension> >& coeffs, D d)
      : gv_(gv), coefficients_(coeffs), mapper_(*gv), d_(d)
  {}

  SpaceOperator(const SpaceOperator& that) = delete;
  SpaceOperator& operator=(const SpaceOperator& that) = delete;

  void apply(const Vector& source, Vector& range, double exponent, double weight=1.) const
  {
    auto endit = gv_->template end<0>();
    for (auto it = gv_->template begin<0>(); it != endit; ++it) {

      int indInside = mapper_.map(*it);

      auto isend = gv_->iend(*it);
      for (auto is = gv_->ibegin(*it); is != isend; ++is) {

        if (is->neighbor()) {

          auto outside = is->outside();
          int indOutside = mapper_.map(*outside);

          if (indInside < indOutside) {

            const double u_s = source[indInside];
            const double u_n = source[indOutside];

            const double faceVolume = is->geometry().volume();
            const auto normal = is->centerUnitOuterNormal();

            const double pow_u_s = ((u_s > 0) - (u_s < 0)) * pow(fabs(u_s), exponent);
            const double pow_u_n = ((u_n > 0) - (u_n < 0)) * pow(fabs(u_n), exponent);
            double convectiveFlux = (normal * coefficients_->v) * (pow_u_s + pow_u_n) / 2;

            const double interfaceFlux = ((u_s - u_n) / (2 * coefficients_->lambda) + convectiveFlux) * faceVolume;

            range[indInside] += interfaceFlux / it->geometry().volume() * weight;
            range[indOutside] -= interfaceFlux / outside->geometry().volume() * weight;

          }
        }
      }
    }
  }

  std::size_t dimSource()
  {
    return mapper_.size();
  }

private:
  std::shared_ptr<GV> gv_;
  std::shared_ptr<Coefficients<GV::dimension> > coefficients_;
  Mapper mapper_;
  D d_;

#ifdef AS_LIB
public:
  static void export_(const char* classname)
  {
    using boost::python::class_;
    using boost::python::no_init;

    class_<SpaceOperator, std::shared_ptr<SpaceOperator>, boost::noncopyable>(classname, no_init)
        .def("apply", &SpaceOperator::apply)
        .add_property("dimSource", &SpaceOperator::dimSource)
        .add_property("dimRange", &SpaceOperator::dimSource)
    ;
  }
#endif

};


#endif
