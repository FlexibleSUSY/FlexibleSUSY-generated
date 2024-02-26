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
 * @file NMSSM_amm.cpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "NMSSM_amm.hpp"
#include "NMSSM_mass_eigenstates.hpp"

#include "cxx_qft/NMSSM_qft.hpp"
#include "cxx_qft/NMSSM_particle_aliases.hpp"
#include "NMSSM_FFV_form_factors.hpp"

#include "amm_loop_functions.hpp"
#include "lowe.h"
#include "wrappers.hpp"
#include "numerics2.hpp"
#include "logger.hpp"
#include "scan.hpp"
#include "spectrum_generator_settings.hpp"

using namespace flexiblesusy;
using namespace NMSSM_cxx_diagrams;

namespace {

constexpr double amm_threshold_1l = 0.5;          ///< threshold for inclusion of 1-loop contributions
constexpr double amm_threshold_2l_qed = 1.25;     ///< threshold for inclusion of 2-loop QED logarithmic contributions
constexpr double amm_threshold_2l_barrzee = 1.75; ///< threshold for inclusion of 2-loop Barr-Zee contributions
constexpr double amm_threshold_3l = 2.25;         ///< threshold for inclusion of 3-loop contributions

template <typename Lepton>
double calculate_qed_1l_log(const context_base&, double);

/**
 * @class AMMBarrZeeFermionLoop
 *
 * @brief A template that calculates the contributions to the
 *        anomalous magnetic dipole moment of a given particle in a
 *        Barr-Zee diagram with a fermion loop and scalar and photon
 *        exchange particles (see arXiv:1502.04199 Figure 2, (1)).
 *
 * @tparam Args Specifies in order the field of which to calculate the
 *              electric magnetic moment, the photon emitter and the
 *              exchange particle in a Barr-Zee diagram where the
 *              photon emitter is a Fermion and one exchange particle
 *              a scalar particle (the other is a photon).
 *
 * This template evaluates the contribution to the magnetic
 * dipole moment of a Barr-Zee diagram with fields given by
 * \a Args.
 */
template<typename Lepton, class PhotonEmitter, class ExchangeParticle>
struct AMMBarrZeeFermionLoop {
   static double value(const typename field_indices<Lepton>::type& indices,
                       const context_base& context,
                       double ml_pole);
};

/**
 * @class AMMBarrZeeFermionLoopZ
 *
 * @brief A template that calculates the contributions to the
 *        anomalous magnetic dipole moment of a given particle in a
 *        Barr-Zee diagram with a fermion loop and scalar and Z-boson
 *        exchange particles.
 *
 * @tparam Args Specifies in order the field of which to calculate the
 *              electric magnetic moment, the photon emitter and the
 *              exchange particle in a Barr-Zee diagram where the
 *              photon emitter is a Fermion and one exchange particle
 *              a scalar particle (the other is a Z-boson).
 *
 * This template evaluates the contribution to the magnetic
 * dipole moment of a Barr-Zee diagram with fields given by
 * \a Args.
 */
template<typename Lepton, class PhotonEmitter, class ExchangeParticle>
struct AMMBarrZeeFermionLoopZ {
   static double value(const typename field_indices<Lepton>::type& indices,
                       const context_base& context,
                       double ml_pole);
};

/**
 * @class AMMBarrZeeScalarLoop
 *
 * @brief A template that calculates the contributions to the
 *        anomalous dipole moment of a given particle in a Barr-Zee
 *        diagram with a scalar loop and scalar and photon exchange
 *        particles (see arXiv:1502.04199 Figure 2, (2)).
 *
 * @tparam Args Specifies in order the field of which to calculate the
 *              magnetic dipole moment, the photon emitter and the
 *              exchange particle in a Barr-Zee diagram where the
 *              photon emitter is a scalar particle and one exchange
 *              particle a scalar particle (the other is a photon).
 *
 * This template evaluates the contribution to the magnetic
 * dipole moment of a Barr-Zee diagram with fields given by
 * \a Args.
 */
template<typename Lepton, class PhotonEmitter, class ExchangeParticle>
struct AMMBarrZeeScalarLoop {
   static double value(const typename field_indices<Lepton>::type& indices,
                       const context_base& context,
                       double ml_pole);
};

/**
 * @class AMMBarrZeeVectorLoop
 *
 * @brief A template that calculates the contributions to the
 *        anomalous dipole moment of a given particle in a Barr-Zee
 *        diagram with a vector boson loop and scalar and photon
 *        exchange particles (see arXiv:1502.04199 Figure 2, (3)).
 *
 * @tparam Args Specifies in order the field of which to calculate the
 *              magnetic dipole moment, the photon emitter and the
 *              exchange particle in a Barr-Zee diagram where the
 *              photon emitter is a scalar particle and one exchange
 *              particle a scalar particle (the other is a photon).
 *
 * This template evaluates the contribution to the magnetic
 * dipole moment of a Barr-Zee diagram with fields given by
 * \a Args.
 */
template<typename Lepton, class PhotonEmitter, class ExchangeParticle>
struct AMMBarrZeeVectorLoop {
   static double value(const typename field_indices<Lepton>::type& indices,
                       const context_base& context,
                       double ml_pole);
};

/// returns mass of lightest electrically charged BSM particle
double get_MLCP(const NMSSM_mass_eigenstates_interface& model)
{
   return Min(model.get_MSd().tail<6>().minCoeff(), model.get_MSu().tail<6>().
      minCoeff(), model.get_MSe().tail<6>().minCoeff(), model.get_MHpm().tail<1>()
      .minCoeff(), model.get_MCha().tail<2>().minCoeff());
}

void run_to_MSUSY(NMSSM_mass_eigenstates& model)
{
   const double precision_goal = model.get_precision();
   const int Nmax = static_cast<int>(-std::log10(precision_goal) * 10);
   int it = 0;
   double precision = 0.;
   double MSUSY_old = 0., MSUSY_new = 0.;

   VERBOSE_MSG("NMSSM_amm: iterate to run to MSUSY ...");

   do {
      MSUSY_new = get_MLCP(model);

      if (MSUSY_new <= 0.)
         return;

      VERBOSE_MSG("NMSSM_amm:    running to new MSUSY = " << MSUSY_new);

      model.run_to(MSUSY_new);
      model.calculate_DRbar_masses();
      model.reorder_DRbar_masses();

      precision = MaxRelDiff(MSUSY_old, MSUSY_new);
      MSUSY_old = MSUSY_new;
   } while (precision > precision_goal && ++it < Nmax);

   VERBOSE_MSG("NMSSM_amm: iteration finished. MSUSY = " << MSUSY_new);

   if (precision > precision_goal) {
      WARNING("NMSSM_amm: run_to_MSUSY(): "
              "Iteration did not converge after " << Nmax
              << " iterations (precision goal: " << precision_goal << ")");
   }
}

template <typename Lepton>
double lepton_pole_mass(const softsusy::QedQcd& qedqcd, int idx)
{
   return qedqcd.displayLeptonPoleMass(idx);
}

/// 1-loop scalar/fermionic contributions
template <typename Lepton>
double calculate_amm_impl_1l(const context_base& context, const NMSSM_mass_eigenstates& model, double ml_pole, int idx)
{
   using namespace NMSSM_cxx_diagrams::fields;

   const auto form_factors = NMSSM_FFV_form_factors::calculate_form_factors<Lepton,Lepton,VP>(idx, idx, model, true);

   // Im[F2]/Re[F2] ~ 0
   if (!is_zero((form_factors[2] + form_factors[3]).real())) {
      const double ratio = (form_factors[2] + form_factors[3]).imag()/(form_factors[2] + form_factors[3]).real();
      if (!is_zero(ratio)) {
         ERROR("Error in the g-2 calculation! Form factor F2 should be real: Im[Fe]/Re[F2] = " + std::to_string(ratio));
         return std::numeric_limits<double>::quiet_NaN();
      }
   }
   // Im[F2] ~ 0
   else {
      if (!is_zero((form_factors[2] + form_factors[3]).imag())) {
         ERROR("Error in the g-2 calculation! Form factor F2 should be real: Im[F2] = " + std::to_string((form_factors[2] + form_factors[3]).imag()));
         return std::numeric_limits<double>::quiet_NaN();
      }
   }

   const double result =
      // vector form factor
      0.5*(form_factors[2] + form_factors[3]).real()
      // reinstate the mass that was factored out in definition of form factor
      * context.mass<Lepton>({idx})
      // factor out e/(2*lepton pole mass)
      / (unit_charge(context)/(2.*ml_pole))
      // definition of photon momentum flow for g-2 is opposite than for form factors
      * (-1.);

   return result;
}

/// 2-loop QED logarithmic contributions
template <typename Lepton>
double calculate_amm_impl_2l_qed(const context_base& context, double ml_pole, double amu_1l)
{
   return amu_1l*calculate_qed_1l_log<Lepton>(context, ml_pole);
}

/// 2-loop Barr-Zee contributions
template <typename Lepton>
double calculate_amm_impl_2l_barrzee(const context_base& context, double ml_pole, int idx)
{
   using namespace NMSSM_cxx_diagrams::fields;

   double val = 0.;

   val += 1.000000000000000 * AMMBarrZeeFermionLoop<Lepton,typename bar<Cha>::type, Ah>::value({idx}, context, ml_pole);
   val += 3.000000000000000 * AMMBarrZeeFermionLoop<Lepton,typename bar<Fd>::type, Ah>::value({idx}, context, ml_pole);
   val += 1.000000000000000 * AMMBarrZeeFermionLoop<Lepton,typename bar<Fe>::type, Ah>::value({idx}, context, ml_pole);
   val += 3.000000000000000 * AMMBarrZeeFermionLoop<Lepton,typename bar<Fu>::type, Ah>::value({idx}, context, ml_pole);
   val += 1.000000000000000 * AMMBarrZeeScalarLoop<Lepton,typename conj<Hpm>::type, Ah>::value({idx}, context, ml_pole);
   val += 3.000000000000000 * AMMBarrZeeScalarLoop<Lepton,typename conj<Sd>::type, Ah>::value({idx}, context, ml_pole);
   val += 1.000000000000000 * AMMBarrZeeScalarLoop<Lepton,typename conj<Se>::type, Ah>::value({idx}, context, ml_pole);
   val += 3.000000000000000 * AMMBarrZeeScalarLoop<Lepton,typename conj<Su>::type, Ah>::value({idx}, context, ml_pole);
   val += 1.000000000000000 * AMMBarrZeeFermionLoop<Lepton,typename bar<Cha>::type, hh>::value({idx}, context, ml_pole);
   val += 3.000000000000000 * AMMBarrZeeFermionLoop<Lepton,typename bar<Fd>::type, hh>::value({idx}, context, ml_pole);
   val += 1.000000000000000 * AMMBarrZeeFermionLoop<Lepton,typename bar<Fe>::type, hh>::value({idx}, context, ml_pole);
   val += 3.000000000000000 * AMMBarrZeeFermionLoop<Lepton,typename bar<Fu>::type, hh>::value({idx}, context, ml_pole);
   val += 1.000000000000000 * AMMBarrZeeScalarLoop<Lepton,typename conj<Hpm>::type, hh>::value({idx}, context, ml_pole);
   val += 3.000000000000000 * AMMBarrZeeScalarLoop<Lepton,typename conj<Sd>::type, hh>::value({idx}, context, ml_pole);
   val += 1.000000000000000 * AMMBarrZeeScalarLoop<Lepton,typename conj<Se>::type, hh>::value({idx}, context, ml_pole);
   val += 3.000000000000000 * AMMBarrZeeScalarLoop<Lepton,typename conj<Su>::type, hh>::value({idx}, context, ml_pole);
   val += 1.000000000000000 * AMMBarrZeeVectorLoop<Lepton,typename conj<VWm>::type, hh>::value({idx}, context, ml_pole);
   val += 1.000000000000000 * AMMBarrZeeFermionLoopZ<Lepton,typename bar<Cha>::type, Ah>::value({idx}, context, ml_pole);
   val += 3.000000000000000 * AMMBarrZeeFermionLoopZ<Lepton,typename bar<Fd>::type, Ah>::value({idx}, context, ml_pole);
   val += 1.000000000000000 * AMMBarrZeeFermionLoopZ<Lepton,typename bar<Fe>::type, Ah>::value({idx}, context, ml_pole);
   val += 3.000000000000000 * AMMBarrZeeFermionLoopZ<Lepton,typename bar<Fu>::type, Ah>::value({idx}, context, ml_pole);
   val += 1.000000000000000 * AMMBarrZeeFermionLoopZ<Lepton,typename bar<Cha>::type, hh>::value({idx}, context, ml_pole);
   val += 3.000000000000000 * AMMBarrZeeFermionLoopZ<Lepton,typename bar<Fd>::type, hh>::value({idx}, context, ml_pole);
   val += 1.000000000000000 * AMMBarrZeeFermionLoopZ<Lepton,typename bar<Fe>::type, hh>::value({idx}, context, ml_pole);
   val += 3.000000000000000 * AMMBarrZeeFermionLoopZ<Lepton,typename bar<Fu>::type, hh>::value({idx}, context, ml_pole);


   return val;
}

template <typename Lepton>
double calculate_amm_impl(const NMSSM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Spectrum_generator_settings& settings, int idx)
{
   VERBOSE_MSG("NMSSM_amm: calculating a_mu at Q = " << model.get_scale());

   using namespace NMSSM_cxx_diagrams::fields;

   const context_base context{model};
   const double ml_pole = lepton_pole_mass<Lepton>(qedqcd, idx);
   const double loops = settings.get(Spectrum_generator_settings::calculate_amm);

   double val = 0.;

   if (loops > amm_threshold_1l) {
      const double amu_1l = calculate_amm_impl_1l<Lepton>(context, model, ml_pole, idx);

      // 1L
      val += amu_1l;

      // 2L QED logarithms
      if (loops > amm_threshold_2l_qed) {
         val += calculate_amm_impl_2l_qed<Lepton>(context, ml_pole, amu_1l);
      }
   }

   // 2L Barr-Zee contributions
   if (loops > amm_threshold_2l_barrzee) {
      val += calculate_amm_impl_2l_barrzee<Lepton>(context, ml_pole, idx);
   }

   return val;
}

} // anonymous namespace

template <typename Lepton>
double NMSSM_amm::calculate_amm(const NMSSM_mass_eigenstates& model_, const softsusy::QedQcd& qedqcd, const Spectrum_generator_settings& settings, int idx)
{
   NMSSM_mass_eigenstates model(model_);

   VERBOSE_MSG("NMSSM_amm: starting calculation of a_mu ...");

   try {
      run_to_MSUSY(model);
      model.get_physical().clear();
   } catch (const Error& e) {
      ERROR("NMSSM_amm:" << e.what_detailed());
      return std::numeric_limits<double>::quiet_NaN();
   }

   return calculate_amm_impl<Lepton>(model, qedqcd, settings, idx);
}

// spliting templated function into header and cpp
// https://isocpp.org/wiki/faq/templates#separate-template-fn-defn-from-decl
template double NMSSM_amm::calculate_amm<NMSSM_cxx_diagrams::fields::Fe>(const NMSSM_mass_eigenstates&, const softsusy::QedQcd&, const Spectrum_generator_settings&, int);

template <typename Lepton>
double NMSSM_amm::calculate_amm_uncertainty(const NMSSM_mass_eigenstates& model_, const softsusy::QedQcd& qedqcd, const Spectrum_generator_settings& settings, int idx)
{
   NMSSM_mass_eigenstates model(model_);

   VERBOSE_MSG("NMSSM_amm: starting calculation of a_mu uncertainty ...");

   try {
      run_to_MSUSY(model);
      model.get_physical().clear();
   } catch (const Error& e) {
      ERROR("NMSSM_amm uncertainty: " << e.what_detailed());
      return std::numeric_limits<double>::quiet_NaN();
   }

   const context_base context{model};
   const double ml_pole = lepton_pole_mass<Lepton>(qedqcd, idx);
   const double loops = settings.get(Spectrum_generator_settings::calculate_amm);

   const double amu_1l = calculate_amm_impl_1l<Lepton>(context, model, ml_pole, idx);
   const double amu_2l_qed = calculate_amm_impl_2l_qed<Lepton>(context, ml_pole, amu_1l);
   const double amu_2l_barrzee = calculate_amm_impl_2l_barrzee<Lepton>(context, ml_pole, idx);
   const double amu_2l = amu_2l_qed + amu_2l_barrzee;

   double damu = 0.;

   if (loops <= amm_threshold_1l) {
      // add 1-loop uncertainty from missing 1-loop contributions
      damu += Abs(amu_1l);
   }

   if (loops <= amm_threshold_2l_qed) {
      // add 2-loop uncertainty from missing 2-loop QED contributions
      damu += Abs(amu_2l_qed);
   }

   if (loops <= amm_threshold_2l_barrzee) {
      // add 2-loop uncertainty from missing 2-loop Barr-Zee contributions
      damu += Abs(amu_2l_barrzee);
   }

   if (loops <= amm_threshold_3l) {
      // add 3-loop uncertainty from missing 3-loop QED contributions
      damu += Abs(amu_2l*calculate_qed_1l_log<Lepton>(context, ml_pole));
   }

   return damu;
}

// spliting templated function into header and cpp
// https://isocpp.org/wiki/faq/templates#separate-template-fn-defn-from-decl
template double NMSSM_amm::calculate_amm_uncertainty<NMSSM_cxx_diagrams::fields::Fe>(const NMSSM_mass_eigenstates&, const softsusy::QedQcd&, const Spectrum_generator_settings&, int);

namespace {

/// returns 1-loop logarithmic term, necessary to calculate the 2-loop QED logarithmic contributions
template <typename Lepton>
double calculate_qed_1l_log(const context_base& context, double m_lepton)
{
   const double MSUSY = Abs(get_MLCP(context.model));
   const double alpha_em = Sqr(Lepton::electricCharge * unit_charge(context))/(4*Pi);
   const double result = alpha_em/(4*Pi) * 16 * FiniteLog(m_lepton/MSUSY);

   return result;
}

template<typename Lepton, class Fermion, class Scalar>
double AMMBarrZeeFermionLoop<
Lepton, Fermion, Scalar
>::value(const typename field_indices<Lepton>::type& indices, const context_base& context, double ml_pole)
{
   using namespace flexiblesusy::amm_loop_functions::two_loop;

   using LeptonVertex = Vertex<
                      Scalar,
                      typename Lepton::lorentz_conjugate,
                      Lepton
                      >;

   using FermionVertex = Vertex<
                         Scalar,
                         Fermion,
                         typename Fermion::lorentz_conjugate
                         >;

   double res = 0.0;

   for (const auto& index1: index_range<LeptonVertex>()) {
      for (const auto& index2: index_range<FermionVertex>()) {
         const auto leptonIndices = LeptonVertex::template indices_of_field<2>(index1);

         if (leptonIndices != indices) {
            continue;
         }

         const auto scalarIndices = LeptonVertex::template indices_of_field<0>(index1);
         const auto fermionIndices = FermionVertex::template indices_of_field<1>(index2);

         if (cxx_diagrams::isSMField<Fermion>(fermionIndices) &&
             cxx_diagrams::isSMField<Scalar>(scalarIndices)) {
            continue;
         }

         if (scalarIndices != FermionVertex::template indices_of_field<0>(index2)) {
            continue;
         }

         if (fermionIndices != FermionVertex::template indices_of_field<2>(index2)) {
            continue;
         }

         const auto leptonVertex = LeptonVertex::evaluate(index1, context);
         const auto fermionVertex = FermionVertex::evaluate(index2, context);

         const auto fermionMass = context.mass<Fermion>(fermionIndices);
         const auto scalarMass = context.mass<Scalar>(scalarIndices);

         static constexpr double fermionChargeCount =
            Fermion::electricCharge / Lepton::electricCharge;

         const std::complex<double> zrm = leptonVertex.right();
         const std::complex<double> zlm = leptonVertex.left();
         const std::complex<double> zrf = fermionVertex.right();
         const std::complex<double> zlf = fermionVertex.left();

         const double coeffA = std::real((zrf + zlf) * (zrm + zlm));
         const double coeffB = std::real((zlf - zrf) * (zrm - zlm));

         const double massRatioSquared = Sqr(fermionMass / scalarMass);

         const double part1 = coeffA == 0 ? 0 : coeffA * BarrZeeLoopFS(massRatioSquared);
         const double part2 = coeffB == 0 ? 0 : coeffB * BarrZeeLoopFP(massRatioSquared);

         const double preFactor = 4 * twoLoop * Sqr(unit_charge(context) * fermionChargeCount)
            * ml_pole / fermionMass;

         res += preFactor * (part1 + part2);
      }
   }

   return res;
}

template<typename Lepton, class Fermion, class Scalar>
double AMMBarrZeeFermionLoopZ<
Lepton, Fermion, Scalar
>::value(const typename field_indices<Lepton>::type& indices, const context_base& context, double ml_pole)
{
   using namespace flexiblesusy::amm_loop_functions::two_loop;

   using LeptonVertex = Vertex<
                      Scalar,
                      typename Lepton::lorentz_conjugate,
                      Lepton
                      >;

   using LeptonZVertex = Vertex<
                      ZBoson,
                      Lepton,
                      typename Lepton::lorentz_conjugate
                      >;

   using FermionVertex = Vertex<
                         Scalar,
                         Fermion,
                         typename Fermion::lorentz_conjugate
                         >;

   using FermionZVertex = Vertex<
                         ZBoson,
                         typename Fermion::lorentz_conjugate,
                         Fermion
                         >;

   double res = 0.0;

   for (const auto& index1: index_range<LeptonVertex>()) {
      for (const auto& index2: index_range<FermionVertex>()) {
         for (const auto& indexLeptonZ: index_range<LeptonZVertex>()) {
            for (const auto& indexFermionZ: index_range<FermionZVertex>()) {
               const auto leptonIndices = LeptonVertex::template indices_of_field<2>(index1);

               if (leptonIndices != indices) {
                  continue;
               }

               if ((leptonIndices != LeptonZVertex::template indices_of_field<2>(indexLeptonZ)) ||
                   (leptonIndices != LeptonZVertex::template indices_of_field<1>(indexLeptonZ))) {
                  continue;
               }

               const auto scalarIndices = LeptonVertex::template indices_of_field<0>(index1);
               const auto fermionIndices = FermionVertex::template indices_of_field<1>(index2);

               if (cxx_diagrams::isSMField<Fermion>(fermionIndices) &&
                   cxx_diagrams::isSMField<Scalar>(scalarIndices)) {
                  continue;
               }

               if (scalarIndices != FermionVertex::template indices_of_field<0>(index2)) {
                  continue;
               }

               if (fermionIndices != FermionVertex::template indices_of_field<2>(index2)) {
                  continue;
               }

               if ((fermionIndices != FermionZVertex::template indices_of_field<2>(indexFermionZ)) ||
                   (fermionIndices != FermionZVertex::template indices_of_field<1>(indexFermionZ))) {
                  continue;
               }

               const auto leptonVertex = LeptonVertex::evaluate(index1, context);
               const auto fermionVertex = FermionVertex::evaluate(index2, context);
               const auto leptonZVertex = LeptonZVertex::evaluate(indexLeptonZ, context);
               const auto fermionZVertex = FermionZVertex::evaluate(indexFermionZ, context);

               const auto fermionMass = context.mass<Fermion>(fermionIndices);
               const auto scalarMass = context.mass<Scalar>(scalarIndices);
               const auto zMass = context.mass<ZBoson>({ });

               static constexpr double fermionChargeCount =
                  Fermion::electricCharge / Lepton::electricCharge;

               const std::complex<double> zrm = leptonVertex.right();
               const std::complex<double> zlm = leptonVertex.left();
               const std::complex<double> zrf = fermionVertex.right();
               const std::complex<double> zlf = fermionVertex.left();
               const std::complex<double> zmzr = leptonZVertex.right();
               const std::complex<double> zmzl = leptonZVertex.left();
               const std::complex<double> zfzr = fermionZVertex.right();
               const std::complex<double> zfzl = fermionZVertex.left();

               const double gvm = 0.5*(zmzr.real() + zmzl.real());
               const double gvf = 0.5*(zfzr.real() + zfzl.real());
               const double coeffA = std::real((zrf + zlf) * (zrm + zlm));
               const double coeffB = std::real((zlf - zrf) * (zrm - zlm));

               const double x = Sqr(fermionMass / scalarMass);
               const double y = Sqr(fermionMass / zMass);

               const double part1 = coeffA == 0 ? 0 : coeffA * BarrZeeLoopFSZ(x, y);
               const double part2 = coeffB == 0 ? 0 : coeffB * BarrZeeLoopFPZ(x, y);

               const double preFactor = 4 * twoLoop * fermionChargeCount * gvm * gvf
                  * ml_pole / fermionMass;

               res += preFactor * (part1 + part2);
            }
         }
      }
   }

   return res;
}

template<typename Lepton, class ChargedScalar, class NeutralScalar>
double AMMBarrZeeScalarLoop<
Lepton, ChargedScalar, NeutralScalar
>::value(const typename field_indices<Lepton>::type& indices, const context_base& context, double ml_pole)
{
   using namespace flexiblesusy::amm_loop_functions::two_loop;

   using LeptonVertex = Vertex<
                      NeutralScalar,
                      typename Lepton::lorentz_conjugate,
                      Lepton
                      >;

   using ScalarVertex = Vertex<
                         NeutralScalar,
                         ChargedScalar,
                         typename ChargedScalar::lorentz_conjugate
                         >;

   double res = 0.0;

   for (const auto& index1: index_range<LeptonVertex>()) {
      for (const auto& index2: index_range<ScalarVertex>()) {
         const auto leptonIndices = LeptonVertex::template indices_of_field<2>(index1);

         if (leptonIndices != indices) {
            continue;
         }

         const auto neutralScalarIndices = LeptonVertex::template indices_of_field<0>(index1);
         const auto chargedScalarIndices = ScalarVertex::template indices_of_field<1>(index2);

         if (cxx_diagrams::isSMField<ChargedScalar>(chargedScalarIndices) &&
             cxx_diagrams::isSMField<NeutralScalar>(neutralScalarIndices)) {
            continue;
         }

         if (neutralScalarIndices != ScalarVertex::template indices_of_field<0>(index2)) {
            continue;
         }

         if (chargedScalarIndices != ScalarVertex::template indices_of_field<2>(index2)) {
            continue;
         }

         const auto leptonVertex = LeptonVertex::evaluate(index1, context);
         const auto scalarVertex = ScalarVertex::evaluate(index2, context);

         const auto chargedScalarMass = context.mass<ChargedScalar>(chargedScalarIndices);
         const auto neutralScalarMass = context.mass<NeutralScalar>(neutralScalarIndices);

         static constexpr double scalarChargeCount =
            ChargedScalar::electricCharge / Lepton::electricCharge;

         const std::complex<double> zrm = leptonVertex.right();
         const std::complex<double> zlm = leptonVertex.left();
         const std::complex<double> zss = scalarVertex.value();

         const double coeff = std::real((zrm + zlm) * zss);

         const double massRatioSquared = Sqr(chargedScalarMass / neutralScalarMass);

         const double preFactor = 4 * twoLoop * Sqr(unit_charge(context) * scalarChargeCount)
            * ml_pole / Sqr(neutralScalarMass);

         res += preFactor * coeff * BarrZeeLoopS(massRatioSquared);
      }
   }

   return res;
}

template<typename Lepton, class Vector, class Scalar>
double AMMBarrZeeVectorLoop<
Lepton, Vector, Scalar
>::value(const typename field_indices<Lepton>::type& indices, const context_base& context, double ml_pole)
{
   using namespace flexiblesusy::amm_loop_functions::two_loop;

   using LeptonVertex = Vertex<
                      Scalar,
                      typename Lepton::lorentz_conjugate,
                      Lepton
                      >;

   using VectorVertex = Vertex<
                         Scalar,
                         Vector,
                         typename Vector::lorentz_conjugate
                         >;

   double res = 0.0;

   for (const auto& index1: index_range<LeptonVertex>()) {
      for (const auto& index2: index_range<VectorVertex>()) {
         const auto leptonIndices = LeptonVertex::template indices_of_field<2>(index1);

         if (leptonIndices != indices) {
            continue;
         }

         const auto scalarIndices = LeptonVertex::template indices_of_field<0>(index1);
         const auto vectorIndices = VectorVertex::template indices_of_field<1>(index2);

         if (cxx_diagrams::isSMField<Vector>(vectorIndices) &&
             cxx_diagrams::isSMField<Scalar>(scalarIndices)) {
            continue;
         }

         if (scalarIndices != VectorVertex::template indices_of_field<0>(index2)) {
            continue;
         }

         if (vectorIndices != VectorVertex::template indices_of_field<2>(index2)) {
            continue;
         }

         const auto leptonVertex = LeptonVertex::evaluate(index1, context);
         const auto vectorVertex = VectorVertex::evaluate(index2, context);

         const auto vectorMass = context.mass<Vector>(vectorIndices);
         const auto scalarMass = context.mass<Scalar>(scalarIndices);

         static constexpr double vectorChargeCount =
            Vector::electricCharge / Lepton::electricCharge;

         const std::complex<double> zrm = leptonVertex.right();
         const std::complex<double> zlm = leptonVertex.left();
         const std::complex<double> zsv = vectorVertex.value();

         const double coeff = std::real((zrm + zlm) * zsv);

         const double massRatioSquared = Sqr(vectorMass / scalarMass);

         const double preFactor = - 2 * twoLoop * Sqr(unit_charge(context) * vectorChargeCount)
            * ml_pole / Sqr(vectorMass);

         res += preFactor * coeff * BarrZeeLoopV(massRatioSquared);
      }
   }

   return res;
}

} // anonymous namespace
