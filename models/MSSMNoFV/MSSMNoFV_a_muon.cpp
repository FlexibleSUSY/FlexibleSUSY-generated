// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================


/**
 * @file MSSMNoFV_a_muon.cpp
 *
 * This file was generated with FlexibleSUSY 2.6.2 and SARAH 4.14.5 .
 */

#include "MSSMNoFV_a_muon.hpp"
#include "MSSMNoFV_mass_eigenstates.hpp"

#include "cxx_qft/MSSMNoFV_qft.hpp"
#include "MSSMNoFV_FFV_form_factors.hpp"

#include "lowe.h"
#include "wrappers.hpp"
#include "numerics2.hpp"
#include "logger.hpp"
#include "dilog.hpp"

#define DERIVEDPARAMETER(p) context.model.p()

using namespace flexiblesusy;
using namespace MSSMNoFV_cxx_diagrams;

using Muon = fields::Fm;
using VZ = fields::VZ;

namespace {

double get_QED_2L(context_base&, const softsusy::QedQcd&);
double BarZeeLoopFPS(double);
double BarZeeLoopFS(double);
double BarZeeLoopS(double);
double BarZeeLoopV(double);

/**
 * @class AMuonBarZeeFermionLoop
 * @brief A template that calculate contributions to the
 *        anomalous magnetic dipole moment of a given particle in
 *        a BarZee diagram with a fermion loop and scalar/Photon exchange
 *        particle. (for diagram type see arXiv:1502.04199 Figure 2, (1))
 * @tparam Args Specifies in order the field of which to
 *              calculate the electric magnetic moment,
 *              the photon emitter and the exchange particle
 *              in a BarZee diagram where the photon emitter
 *              is a Fermion and the exchange particle a scalar particle.
 *
 * This template evaluates the contribution to the magnetic
 * dipole moment of a BarZee diagram with fields given by
 * \a Args.
 */
template<class PhotonEmitter, class ExchangeParticle>
struct AMuonBarZeeFermionLoop {
   static double value(const typename field_indices<Muon>::type& indices,
                       const context_base& context,
                       const softsusy::QedQcd& qedqcd);
};

/**
 * @class AMuonBarZeeFermionZLoop
 * @brief A template that calculate contributions to the
 *        anomalous magnetic dipole moment of a given particle in
 *        a BarZee diagram with a fermion loop and scalar/Z-boson exchange
 *        particle.
 * @tparam Args Specifies in order the field of which to
 *              calculate the electric magnetic moment,
 *              the photon emitter and the exchange particle
 *              in a BarZee diagram where the photon emitter
 *              is a Fermion and the exchange particle a scalar particle.
 *
 * This template evaluates the contribution to the magnetic
 * dipole moment of a BarZee diagram with fields given by
 * \a Args.
 */
template<class PhotonEmitter, class ExchangeParticle>
struct AMuonBarZeeFermionLoopZ {
   static double value(const typename field_indices<Muon>::type& indices,
                       const context_base& context,
                       const softsusy::QedQcd& qedqcd);
};

/**
 * @class AMuonBarZeeScalarLoop
 * @brief A template that calculate contributions to the
 *        anomalous dipole moment of a given particle in
 *        a BarZee diagram with a scalar loop and scalar exchange
 *        particle. (for diagram type see arXiv:1502.04199 Figure 2, (2))
 * @tparam Args Specifies in order the field of which to
 *              calculate the magnetic dipole moment,
 *              the photon emitter and the exchange particle
 *              in a BarZee diagram where the photon emitter
 *              is a scalar particle and the exchange particle a scalar particle.
 *
 * This template evaluates the contribution to the magnetic
 * dipole moment of a BarZee diagram with fields given by
 * \a Args.
 */
template<class PhotonEmitter, class ExchangeParticle>
struct AMuonBarZeeScalarLoop {
   static double value(const typename field_indices<Muon>::type& indices,
                       const context_base& context,
                       const softsusy::QedQcd& qedqcd);
};

/**
 * @class AMuonBarZeeVectorLoop
 * @brief A template that calculate contributions to the
 *        anomalous dipole moment of a given particle in
 *        a BarZee diagram with a scalar loop and scalar exchange
 *        particle. (for diagram type see arXiv:1502.04199 Figure 2, (3))
 * @tparam Args Specifies in order the field of which to
 *              calculate the magnetic dipole moment,
 *              the photon emitter and the exchange particle
 *              in a BarZee diagram where the photon emitter
 *              is a scalar particle and the exchange particle a scalar particle.
 *
 * This template evaluates the contribution to the magnetic
 * dipole moment of a BarZee diagram with fields given by
 * \a Args.
 */
template<class PhotonEmitter, class ExchangeParticle>
struct AMuonBarZeeVectorLoop {
   static double value(const typename field_indices<Muon>::type& indices,
                       const context_base& context,
                       const softsusy::QedQcd& qedqcd);
};

/**
    Loop function needed for the calculation of the BarZeeDiagram with
    fermion loop. (arXiv:1502.04199 Equation 26)

    @param respective squared mass ratio.
*/

double BarZeeLoopFPS(double m) {
   // @todo(alex): check stability for small values
   if (m == 0) {
      return 0;
   }

   if (m == 0.25) {
      return ln2;
   }

   double r1, theta1, r2, theta2;
   if (m <= 0.25) {
      const double y = std::sqrt(1.0-4.0*m);
      r1 = 1.0 - (1.0 - y)/(2*m);
      r2 = 1.0 - (1.0 + y)/(2*m);
      if (r1 > 0) {
         theta1 = 0.0;
      }
      else {
         r1 *= -1;
         theta1 = Pi;
      }
      if (r2 > 0) {
         theta2 = 0.0;
      }
      else {
         r2 *= -1;
         theta2 = Pi;
      }
      return m / y * (dilog(std::polar(r1, theta1)).real() - dilog(std::polar(r2, theta2)).real());
   }
   else {
      const double y = std::sqrt(-1.0+4.0*m);
      const double real = 1 - 1/(2*m);
      r1 = r2 = 1.0;

      theta1 = std::atan2(y/(2*m), real);
      theta2 = std::atan2(-y/(2*m), real);

      return m / y * (dilog(std::polar(r1, theta1)).imag() - dilog(std::polar(r2, theta2)).imag());
   }
}

/**
    Loop function needed for the calculation of the BarZeeDiagram with
    fermion loop. (arXiv:1502.04199 Equation 25)

    @param respective squared mass ratio.
*/

double BarZeeLoopFS(double m)
{
    if (m == 0) {
        return 0;
    }

    return (2*m - 1) * BarZeeLoopFPS(m) - m * (2 + std::log(m));
}

/**
    Loop function needed for the calculation of the BarZeeDiagram with
    fermion loop. (arXiv:1502.04199 Equation 27)

    @param respective squared mass ratio.
*/

double BarZeeLoopS(const double m)
{
   return 1 + std::log(m)/2 - BarZeeLoopFPS(m);
}

/**
    Loop function needed for the calculation of the BarZeeDiagram with
    fermion loop. (arXiv:1502.04199 Equation 28)

    @param respective squared mass ratio.
*/

double BarZeeLoopV(const double m) {
   if (m == 0.25) {
      return 4.75;
   }

   double j;

   double r1, theta1, r2, theta2;

   if (m < 0.25) {
      const double y = std::sqrt(1.0-4.0*m);
      r1 = 2 / (1-y);
      theta1 = std::atan2(0, r1);
      r1 = std::abs(r1);
      r2 = 2 / (1+y);
      theta2 = std::atan2(0, r2);
      r2 = std::abs(r2);

      j = 1/y * (std::log(std::abs(y-1) / (y+1)) * std::log(m) + dilog(std::polar(r1, theta1)).real() - dilog(std::polar(r2, theta2)).real());
   }
   else {
      const double y = std::sqrt(-1.0+4.0*m);
      r1 = r2 = std::sqrt(1/m);
      const double real = 1/(2*m);
      const double imag = y/(2*m);
      theta1 = std::atan2(imag, real);
      theta2 = std::atan2(-imag, real);
      j = 1/y * ( (std::atan2(y, -1) - std::atan2(y, 1)) * std::log(m) + dilog(std::polar(r1, theta1)).imag() - dilog(std::polar(r2, theta2)).imag());
   }

   return BarZeeLoopS(m) + 15.0/2.0 * m * (2.0 + std::log(m)) + m/2 * (19 - 12*m) * j - 9*m * BarZeeLoopFPS(m);
}

double get_MSUSY(const MSSMNoFV_mass_eigenstates_interface& model)
{
   return Min(model.get_MSd().tail<2>().minCoeff(), model.get_MSu().tail<2>().
      minCoeff(), model.get_MSe().tail<2>().minCoeff(), model.get_MSm().tail<2>().
      minCoeff(), model.get_MStau().tail<2>().minCoeff(), model.get_MSs().tail<2>(
      ).minCoeff(), model.get_MSc().tail<2>().minCoeff(), model.get_MSb().tail<2>(
      ).minCoeff(), model.get_MSt().tail<2>().minCoeff(), model.get_MHpm().tail<1>
      ().minCoeff(), model.get_MCha().tail<2>().minCoeff());
}

void run_to_MSUSY(MSSMNoFV_mass_eigenstates& model)
{
   const double precision_goal = model.get_precision();
   const int Nmax = static_cast<int>(-std::log10(precision_goal) * 10);
   int it = 0;
   double precision = 0.;
   double MSUSY_old = 0., MSUSY_new = 0.;

   VERBOSE_MSG("MSSMNoFV_a_muon: iterate to run to MSUSY ...");

   do {
      MSUSY_new = get_MSUSY(model);

      if (MSUSY_new <= 0.)
         return;

      VERBOSE_MSG("MSSMNoFV_a_muon:    running to new MSUSY = " << MSUSY_new);

      model.run_to(MSUSY_new);
      model.calculate_DRbar_masses();

      precision = MaxRelDiff(MSUSY_old, MSUSY_new);
      MSUSY_old = MSUSY_new;
   } while (precision > precision_goal && ++it < Nmax);

   VERBOSE_MSG("MSSMNoFV_a_muon: iteration finished. MSUSY = " << MSUSY_new);

   if (precision > precision_goal) {
      WARNING("MSSMNoFV_a_muon: run_to_MSUSY(): "
              "Iteration did not converge after " << Nmax
              << " iterations (precision goal: " << precision_goal << ")");
   }
}

double muonPhysicalMass(const softsusy::QedQcd& qedqcd)
{
   return qedqcd.displayPoleMmuon();
}

double calculate_a_muon_impl(const MSSMNoFV_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)
{
   VERBOSE_MSG("MSSMNoFV_a_muon: calculating a_mu at Q = " << model.get_scale());

   context_base context{ model };

   using namespace MSSMNoFV_cxx_diagrams::fields;

   const auto form_factors = MSSMNoFV_FFV_form_factors::calculate_Fm_Fm_VP_form_factors(model, true);

   if (!is_zero((form_factors[2] + form_factors[3]).imag())) {
      ERROR("Error in the g-2 calculation! Form factor F2 should be real");
      return std::numeric_limits<double>::quiet_NaN();
   }

   double val =
      // vector form factor
      1./2.*(form_factors[2] + form_factors[3]).real()
      // reinstate the mass that was factored out in definition of form factor
      * context.mass<Muon>({})
      // factor out e/(2*muon pole mass)
      / (unit_charge(context_base{model})/(2.*muonPhysicalMass(qedqcd)))
      // definition of photon momentum flow for g-2 is opposite than for form factors
      * (-1.);

   // add 2-loop QED logarithms
   val *= 1. + get_QED_2L(context, qedqcd);

   std::complex<double> valBarZee {0., 0.};

   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Cha>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Fb>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Fc>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Fd>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Fe>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Fm>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Fs>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Ft>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Ftau>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Fu>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Hpm>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Sb>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Sc>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Sd>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Se>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Sm>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Ss>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<St>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Stau>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Su>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Cha>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Fb>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Fc>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Fd>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Fe>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Fm>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Fs>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Ft>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Ftau>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoop<typename bar<Fu>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Hpm>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Sb>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Sc>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Sd>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Se>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Sm>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Ss>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<St>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Stau>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeScalarLoop<typename conj<Su>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeVectorLoop<typename conj<VWm>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Cha>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Fb>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Fc>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Fd>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Fe>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Fm>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Fs>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Ft>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Ftau>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Fu>::type, Ah>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Cha>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Fb>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Fc>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Fd>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Fe>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Fm>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Fs>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Ft>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {1.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Ftau>::type, hh>::value({}, context, qedqcd);
   valBarZee += std::complex<double> {3.000000000000000, 0} * AMuonBarZeeFermionLoopZ<typename bar<Fu>::type, hh>::value({}, context, qedqcd);


   if (!is_zero(valBarZee.imag())) {
      ERROR("Error in the g-2 calculation! Form factor F2 should be real");
      return std::numeric_limits<double>::quiet_NaN();
   } else {
      // @todo: enable once the calculation of Barr-Zee contributions is validated
      // val += valBarZee.real();
   }

   return val;
}

/// generates array with N scales from mean/factor to mean*factor
template <int N>
std::array<double,N> generate_scales(double mean, double factor)
{
   static_assert(N > 1, "N must be greater than 1!");

   const double start = mean / factor, stop = mean * factor;
   std::array<double,N> scales;

   scales[0] = start;

   for (int i = 1; i < (N-1); i++)
      scales[i] = std::exp(std::log(start) + (std::log(stop) - std::log(start))*i / N);

   scales[N-1] = stop;

   return scales;
}

/// returns minimum and maximum a_mu when scale is varied by a factor 2
std::pair<double,double> vary_scale(const MSSMNoFV_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)
{
   auto scales = generate_scales<7>(model.get_scale(), 2.);

   std::transform(scales.begin(), scales.end(), scales.begin(),
                  [&model,&qedqcd] (double scale) {
                     double amu = 0.;
                     try {
                        auto m = model;
                        m.run_to(scale);
                        m.get_physical().clear();
                        m.calculate_DRbar_masses();
                        m.solve_ewsb();
                        m.calculate_MFm_pole();
                        amu = calculate_a_muon_impl(m, qedqcd);
                     }
                     catch(const Error& e) {
                        ERROR("MSSMNoFV_a_muon: scale variation: " << e.what_detailed());
                     }
                     return amu;
                  });

   const auto minmax = std::minmax_element(scales.cbegin(), scales.cend());

   return std::make_pair(*(minmax.first), *(minmax.second));
}

} // anonymous namespace

double MSSMNoFV_a_muon::calculate_a_muon(const MSSMNoFV_mass_eigenstates& model_, const softsusy::QedQcd& qedqcd)
{
   MSSMNoFV_mass_eigenstates model(model_);

   VERBOSE_MSG("MSSMNoFV_a_muon: starting calculation of a_mu ...");

   try {
      run_to_MSUSY(model);
      model.get_physical().clear();
   } catch (const Error& e) {
      ERROR("MSSMNoFV_a_muon:" << e.what_detailed());
      return std::numeric_limits<double>::quiet_NaN();
   }

   double m_muon_pole = muonPhysicalMass(qedqcd);

   if (m_muon_pole == 0.0) {
      model.solve_ewsb();
      model.calculate_MFm_pole();
   }

   return calculate_a_muon_impl(model, qedqcd);
}

double MSSMNoFV_a_muon::calculate_a_muon_uncertainty(const MSSMNoFV_mass_eigenstates& model_, const softsusy::QedQcd& qedqcd)
{
   MSSMNoFV_mass_eigenstates model(model_);

   VERBOSE_MSG("MSSMNoFV_a_muon: starting calculation of a_mu uncertainty ...");

   try {
      run_to_MSUSY(model);
      model.get_physical().clear();
   } catch (const Error& e) {
      ERROR("MSSMNoFV_a_muon uncertainty: " << e.what_detailed());
      return std::numeric_limits<double>::quiet_NaN();
   }

   const auto delta_amu_scale_minmax = vary_scale(model, qedqcd);
   const auto delta_amu_scale = std::abs(delta_amu_scale_minmax.second - delta_amu_scale_minmax.first);

   return delta_amu_scale;
}

namespace {
double get_QED_2L(context_base& context, const softsusy::QedQcd& qedqcd)
{
   const double MSUSY = Abs(get_MSUSY(context.model));
   const double m_muon = muonPhysicalMass(qedqcd);
   const double alpha_em = Sqr(Muon::electric_charge * unit_charge(context))/(4*Pi);
   const double qed_2L = alpha_em/(4*Pi) * 16 * FiniteLog(m_muon/MSUSY);

   return qed_2L;
}

template<class Fermion, class Scalar>
double AMuonBarZeeFermionLoop<
Fermion, Scalar
>::value(const typename field_indices<Muon>::type& indices, const context_base& context, const softsusy::QedQcd& qedqcd)
{
   using MuonVertex = Vertex<
                      Scalar,
                      typename Muon::lorentz_conjugate,
                      Muon
                      >;

   using FermionVertex = Vertex<
                         Scalar,
                         Fermion,
                         typename Fermion::lorentz_conjugate
                         >;

   double res = 0.0;
   for (const auto& index1: index_range<MuonVertex>()) {
      for (const auto& index2: index_range<FermionVertex>()) {
          const auto muonIndices = MuonVertex::template indices_of_field<2>(index1);

         if (muonIndices != indices)
            continue;

         const auto scalarIndices = MuonVertex::template indices_of_field<0>(index1);
         const auto fermionIndices = FermionVertex::template indices_of_field<1>(index2);

         if (isSMField<Fermion>(fermionIndices) &&
             isSMField<Scalar>(scalarIndices))
               continue;

         if (scalarIndices != FermionVertex::template indices_of_field<0>(index2))
            continue;

         if (fermionIndices != FermionVertex::template indices_of_field<2>(index2))
            continue;

         const auto muonVertex = MuonVertex::evaluate(index1, context);
         const auto fermionVertex = FermionVertex::evaluate(index2, context);

         const auto fermionMass = context.mass<Fermion>(fermionIndices);
         const auto scalarMass = context.mass<Scalar>(scalarIndices);
         const auto muonMass = context.mass<Muon>(muonIndices);

         const double fermionChargeCount =
            Fermion::electric_charge / Muon::electric_charge;

         const std::complex<double> zrm = muonVertex.right();
         const std::complex<double> zlm = muonVertex.left();
         const std::complex<double> zrf = fermionVertex.right();
         const std::complex<double> zlf = fermionVertex.left();

         const std::complex<double> coeffA = (zrf + zlf) * (zrm + zlm);
         const std::complex<double> coeffB = (zlf - zrf) * (zrm - zlm);

         const double massRatioSquared = Sqr(fermionMass / scalarMass);

         const double part1 = coeffA.real() * BarZeeLoopFS(massRatioSquared);
         const double part2 = coeffB.real() * BarZeeLoopFPS(massRatioSquared);

         const double preFactor = Sqr(unit_charge(context))/(64*Power4(Pi)) * fermionChargeCount * fermionChargeCount
            * muonPhysicalMass(qedqcd) / fermionMass;

         res += preFactor * (part1 + part2);
      }
   }

   return res;
}

template<class Fermion, class Scalar>
double AMuonBarZeeFermionLoopZ<
Fermion, Scalar
>::value(const typename field_indices<Muon>::type& indices, const context_base& context, const softsusy::QedQcd& qedqcd)
{
   using MuonVertex = Vertex<
                      Scalar,
                      typename Muon::lorentz_conjugate,
                      Muon
                      >;

   using MuonZVertex = Vertex<
                      VZ,
                      Muon,
                      typename Muon::lorentz_conjugate
                      >;

   using FermionVertex = Vertex<
                         Scalar,
                         Fermion,
                         typename Fermion::lorentz_conjugate
                         >;

   using FermionZVertex = Vertex<
                         VZ,
                         typename Fermion::lorentz_conjugate,
                         Fermion
                         >;

   double res = 0.0;
   for (const auto& index1: index_range<MuonVertex>()) {
      for (const auto& index2: index_range<FermionVertex>()) {
         for (const auto& indexMuonZ: index_range<MuonZVertex>()) {
            for (const auto& indexFermionZ: index_range<FermionZVertex>()) {

               const auto muonIndices = MuonVertex::template indices_of_field<2>(index1);

               if (muonIndices != indices)
                  continue;

               if ((muonIndices != MuonZVertex::template indices_of_field<2>(indexMuonZ)) ||
                   (muonIndices != MuonZVertex::template indices_of_field<1>(indexMuonZ)))
                  continue;

               const auto scalarIndices = MuonVertex::template indices_of_field<0>(index1);
               const auto fermionIndices = FermionVertex::template indices_of_field<1>(index2);

               if (isSMField<Fermion>(fermionIndices) &&
                   isSMField<Scalar>(scalarIndices))
                  continue;

               if (scalarIndices != FermionVertex::template indices_of_field<0>(index2))
                  continue;

               if (fermionIndices != FermionVertex::template indices_of_field<2>(index2))
                  continue;

               if ((fermionIndices != FermionZVertex::template indices_of_field<2>(indexFermionZ)) ||
                   (fermionIndices != FermionZVertex::template indices_of_field<1>(indexFermionZ)))
                  continue;

               const auto muonVertex = MuonVertex::evaluate(index1, context);
               const auto fermionVertex = FermionVertex::evaluate(index2, context);
               const auto muonZVertex = MuonZVertex::evaluate(indexMuonZ, context);
               const auto fermionZVertex = FermionZVertex::evaluate(indexFermionZ, context);

               const auto fermionMass = context.mass<Fermion>(fermionIndices);
               const auto scalarMass = context.mass<Scalar>(scalarIndices);
               const auto muonMass = context.mass<Muon>(muonIndices);
               const auto zMass = context.mass<VZ>({ });
               const auto ThetaW = DERIVEDPARAMETER(ThetaW);
               const double cw = Cos(ThetaW);
               const double sw = Sin(ThetaW);

               if (is_zero(zMass - scalarMass, 3e-13)) {
                  continue;
               }

               const double fermionChargeCount =
                  Fermion::electric_charge / Muon::electric_charge;

               const std::complex<double> zrm = muonVertex.right();
               const std::complex<double> zlm = muonVertex.left();
               const std::complex<double> zrf = fermionVertex.right();
               const std::complex<double> zlf = fermionVertex.left();
               const std::complex<double> zmzr = muonZVertex.right();
               const std::complex<double> zfzr = fermionZVertex.right();

               const double gvm = (0.5) * (zmzr.real()/unit_charge(context) + sw/cw);
               const double gvf = (0.5) * (zfzr.real()/unit_charge(context) - fermionChargeCount*sw/cw);
               const std::complex<double> coeffA = (zrf + zlf) * (zrm + zlm);
               const std::complex<double> coeffB = (zlf - zrf) * (zrm - zlm);

               const double massRatioSquared = Sqr(fermionMass / scalarMass);

               const double part1 = coeffA.real() * Sqr(scalarMass) / (Sqr(scalarMass) - Sqr(zMass)) * (BarZeeLoopFS(massRatioSquared) - BarZeeLoopFS(Sqr(fermionMass / zMass)));
               const double part2 = coeffB.real() * Sqr(scalarMass) / (Sqr(scalarMass) - Sqr(zMass)) * (BarZeeLoopFPS(massRatioSquared) - BarZeeLoopFPS(Sqr(fermionMass / zMass)));

               const double preFactor = - Sqr(unit_charge(context))/(64*Power4(Pi)) * fermionChargeCount * gvm * gvf
                  * muonPhysicalMass(qedqcd) / fermionMass;

               res += preFactor * (part1 + part2);
            }
         }
      }
   }

   return res;
}

template<class ChargedScalar, class NeutralScalar>
double AMuonBarZeeScalarLoop<
ChargedScalar, NeutralScalar
>::value(const typename field_indices<Muon>::type& indices, const context_base& context, const softsusy::QedQcd& qedqcd)
{
   using MuonVertex = Vertex<
                      NeutralScalar,
                      typename Muon::lorentz_conjugate,
                      Muon
                      >;

   using ScalarVertex = Vertex<
                         NeutralScalar,
                         ChargedScalar,
                         typename ChargedScalar::lorentz_conjugate
                         >;

   double res = 0.0;
   for (const auto& index1: index_range<MuonVertex>()) {
      for (const auto& index2: index_range<ScalarVertex>()) {
          const auto muonIndices = MuonVertex::template indices_of_field<2>(index1);

         if (muonIndices != indices)
            continue;

         const auto neutralScalarIndices = MuonVertex::template indices_of_field<0>(index1);
         const auto chargedScalarIndices = ScalarVertex::template indices_of_field<1>(index2);

         if (isSMField<ChargedScalar>(chargedScalarIndices) &&
             isSMField<NeutralScalar>(neutralScalarIndices))
            {
               continue;
            }

         if (neutralScalarIndices != ScalarVertex::template indices_of_field<0>(index2))
            continue;

         if (chargedScalarIndices != ScalarVertex::template indices_of_field<2>(index2))
            continue;

         const auto muonVertex = MuonVertex::evaluate(index1, context);
         const auto scalarVertex = ScalarVertex::evaluate(index2, context);

         const auto chargedScalarMass = context.mass<ChargedScalar>(chargedScalarIndices);
         const auto neutralScalarMass = context.mass<NeutralScalar>(neutralScalarIndices);
         const auto muonMass = context.mass<Muon>(muonIndices);

         const double scalarChargeCount =
            ChargedScalar::electric_charge / Muon::electric_charge;

         const std::complex<double> zrm = muonVertex.right();
         const std::complex<double> zlm = muonVertex.left();
         const std::complex<double> zss = scalarVertex.value();

         const std::complex<double> coeff = (zrm + zlm) * zss;

         const double massRatioSquared = Sqr(chargedScalarMass / neutralScalarMass);

         const double preFactor = Sqr(unit_charge(context))/(64*Power4(Pi)) * scalarChargeCount * scalarChargeCount
            * muonPhysicalMass(qedqcd) / (neutralScalarMass * neutralScalarMass);

         res += preFactor * coeff.real() * BarZeeLoopS(massRatioSquared);
      }
   }

   return res;
}

template<class Vector, class Scalar>
double AMuonBarZeeVectorLoop<
Vector, Scalar
>::value(const typename field_indices<Muon>::type& indices, const context_base& context, const softsusy::QedQcd& qedqcd)
{

   using MuonVertex = Vertex<
                      Scalar,
                      typename Muon::lorentz_conjugate,
                      Muon
                      >;

   using VectorVertex = Vertex<
                         Scalar,
                         Vector,
                         typename Vector::lorentz_conjugate
                         >;

   double res = 0.0;
   for (const auto& index1: index_range<MuonVertex>()) {
      for (const auto& index2: index_range<VectorVertex>()) {
          const auto muonIndices = MuonVertex::template indices_of_field<2>(index1);

         if (muonIndices != indices)
            continue;

         const auto scalarIndices = MuonVertex::template indices_of_field<0>(index1);
         const auto vectorIndices = VectorVertex::template indices_of_field<1>(index2);

         if (isSMField<Vector>(vectorIndices) &&
             isSMField<Scalar>(scalarIndices))
               continue;

         if (scalarIndices != VectorVertex::template indices_of_field<0>(index2))
            continue;

         if (vectorIndices != VectorVertex::template indices_of_field<2>(index2))
            continue;

         const auto muonVertex = MuonVertex::evaluate(index1, context);
         const auto vectorVertex = VectorVertex::evaluate(index2, context);

         const auto vectorMass = context.mass<Vector>(vectorIndices);
         const auto scalarMass = context.mass<Scalar>(scalarIndices);
         const auto muonMass = context.mass<Muon>(muonIndices);

         const double vectorChargeCount =
            Vector::electric_charge / Muon::electric_charge;

         const std::complex<double> zrm = muonVertex.right();
         const std::complex<double> zlm = muonVertex.left();
         const std::complex<double> zsv = vectorVertex.value();

         const std::complex<double> coeff = (zrm + zlm) * zsv;

         const double massRatioSquared = Sqr(vectorMass / scalarMass);

         const double preFactor = - Sqr(unit_charge(context))/(128*Power4(Pi)) * vectorChargeCount * vectorChargeCount
            * muonPhysicalMass(qedqcd) / (vectorMass * vectorMass);

         res += preFactor * coeff.real() * BarZeeLoopV(massRatioSquared);
      }
   }

   return res;
}

} // anonymous namespace
