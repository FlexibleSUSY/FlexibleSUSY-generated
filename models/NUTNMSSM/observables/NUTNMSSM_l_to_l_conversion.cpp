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
 * @file NUTNMSSM_l_to_l_conversion.cpp
 *
 * This file was generated at Sun 25 Feb 2024 20:57:40 with FlexibleSUSY
 * 2.8.0 and SARAH 4.15.1
 */

#include <valarray>
#include <complex>

#include "NUTNMSSM_mass_eigenstates.hpp"
#include "cxx_qft/NUTNMSSM_qft.hpp"

#include "NUTNMSSM_l_to_l_conversion.hpp"
#include "observables/l_to_l_conversion/settings.hpp"
#include "NUTNMSSM_FFV_form_factors.hpp"

#define F4F(quark, g) embed_photon<quark, Photon>(model, form_factors, g)
#define NPF(name, g) fix_tensors_sign<name>(model, in, out, g)
#define GS(N,Q) parameters.get(LToLConversion_settings::scalar_##N##Q)*m.N/m.Q
#define GV(N,Q) parameters.get(LToLConversion_settings::vector_##N##Q)
#define GT(N,Q) parameters.get(LToLConversion_settings::tensor_##N##Q)

#include "loop_libraries/loop_library.hpp"
#include "cxx_qft/NUTNMSSM_npointfunctions_wilsoncoeffs.hpp"
#include "concatenate.hpp"
#include <limits>
#include <type_traits>
#include <boost/fusion/include/at_key.hpp>
#include "wrappers.hpp"

namespace flexiblesusy {

namespace NUTNMSSM_cxx_diagrams {
namespace npointfunctions {



} // namespace npointfunctions
} // namespace NUTNMSSM_cxx_diagrams

namespace {

std::valarray<std::complex<double>> zero(
   int generationIndex1,
   int generationIndex2,
   const NUTNMSSM_mass_eigenstates&,
   bool){
   std::valarray<std::complex<double>> res {0.0, 0.0, 0.0, 0.0};
   return res;
}

double get_MSUSY(const NUTNMSSM_mass_eigenstates_interface& model)
{
   return Min(model.get_MSd().tail<6>().minCoeff(), model.get_MSu().tail<6>().
      minCoeff(), model.get_MSe().tail<6>().minCoeff(), model.get_MHpm().tail<1>()
      .minCoeff(), model.get_MCha().tail<2>().minCoeff());
}

void run_to_MSUSY(NUTNMSSM_mass_eigenstates& model)
{
   const double precision_goal = model.get_precision();
   const int Nmax = static_cast<int>(-std::log10(precision_goal) * 10);
   int it = 0;
   double precision = 0.;
   double MSUSY_old = 0., MSUSY_new = 0.;

   VERBOSE_MSG("NUTNMSSM_l_to_l_conversion: iterate to run to MSUSY ...");

   do {
      MSUSY_new = get_MSUSY(model);

      if (MSUSY_new <= 0.)
         return;

      VERBOSE_MSG("NUTNMSSM_l_to_l_conversion:    running to new MSUSY = " << MSUSY_new);

      model.run_to(MSUSY_new);
      model.calculate_DRbar_masses();

      precision = MaxRelDiff(MSUSY_old, MSUSY_new);
      MSUSY_old = MSUSY_new;
   } while (precision > precision_goal && ++it < Nmax);

   VERBOSE_MSG("NUTNMSSM_l_to_l_conversion: iteration finished. MSUSY = " << MSUSY_new);

   if (precision > precision_goal) {
      WARNING("NUTNMSSM_l_to_l_conversion: run_to_MSUSY(): "
              "Iteration did not converge after " << Nmax
              << " iterations (precision goal: " << precision_goal << ")");
   }
}

} // anonymous namespace

using namespace NUTNMSSM_cxx_diagrams;
using namespace NUTNMSSM_FFV_form_factors;

namespace NUTNMSSM_l_to_l_conversion {

struct overlap_integrals {
   double D;
   double Vp;
   double Vn;
   double Sp;
   double Sn;
};

template <class Lepton, class Up, class Down>
struct Mass {
   // PDG 2018 data
   double p = 0.938272081, n = 0.939565413;
   // Other masses
   double l, u, c, d, s, b;
   Mass(const context_base& context, int in) {
      l = context.mass<Lepton>({in});
      u = context.mass<Up>({0});
      c = context.mass<Up>({1});
      d = context.mass<Down>({0});
      s = context.mass<Down>({1});
      b = context.mass<Down>({2});
   }
};

overlap_integrals get_overlap_integrals(const Nucleus N, const softsusy::QedQcd& qedqcd) {
   overlap_integrals res;

   // get muon pole mass from input slha file
   const auto muon_pole_mass_5o2 = pow(qedqcd.displayMass(softsusy::mMuon), 5./2.);

   // Tab. 2 of hep-ph/0203110
   switch (N) {
      case Nucleus::Au:
         res.D  = 0.1670 * muon_pole_mass_5o2;
         res.Vp = 0.0859 * muon_pole_mass_5o2;
         res.Vn = 0.1080 * muon_pole_mass_5o2;
         res.Sp = 0.0523 * muon_pole_mass_5o2;
         res.Sn = 0.0610 * muon_pole_mass_5o2;
         return res;
      case Nucleus::Al:
         res.D  = 0.0357 * muon_pole_mass_5o2;
         res.Vp = 0.0159 * muon_pole_mass_5o2;
         res.Vn = 0.0169 * muon_pole_mass_5o2;
         res.Sp = 0.0153 * muon_pole_mass_5o2;
         res.Sn = 0.0163 * muon_pole_mass_5o2;
         return res;
      default:
         throw std::invalid_argument("Unknown nucleus");
   }
}

double get_capture_rate (const Nucleus N) {
   switch (N) {
      case Nucleus::Au:
         return 8.84868e-18;
      case Nucleus::Al:
         return 4.64079e-19;
      default:
         throw std::invalid_argument("Unknown nucleus");
   }
}

typedef std::valarray<std::complex<double>>  (*ffv_function)
   (int, int, const NUTNMSSM_mass_eigenstates&, bool);

typedef std::array<std::complex<double>, 10> (*npf_function)
   (const NUTNMSSM_mass_eigenstates&, const std::array<int,4>&, const std::array<Eigen::Vector4d,0>&);

/**
 * @tparam    Q     Type of a quark field.
 * @tparam    A     Type of a photon field.
 * @param[in] model Mass eigenstates.
 * @param[in] g     Generation index for quarks.
 * @return Left part of qqv coupling, multiplied by -i.
 */
template <class Q, class A>
std::complex<double> left(const NUTNMSSM_mass_eigenstates& model, int g) {
    context_base context {model};
    using vertex = Vertex<typename Q::lorentz_conjugate, Q, A>;
    std::array<int, 2> indices {g, g};
    const auto value =  vertex::evaluate(indices, context);
    return value.left();
}

/**
 * @tparam    Q     Type of a quark field.
 * @tparam    A     Type of a photon field.
 * @param[in] model Mass eigenstates.
 * @param[in] g     Generation index for quarks.
 * @return Right part of qqv coupling, multiplied by -i.
 */
template <class Q, class A>
std::complex<double> right(const NUTNMSSM_mass_eigenstates& model, int g) {
    context_base context {model};
    using vertex = Vertex<typename Q::lorentz_conjugate, Q, A>;
    std::array<int, 2> indices {g, g};
    const auto value =  vertex::evaluate(indices, context);
    return value.right();
}

/**
 * @tparam    Q     Type of a quark field.
 * @tparam    A     Type of a photon field.
 * @tparam    T     Type of a form factors.
 * @param[in] model Mass eigenstates.
 * @param[in] ff    Lepton-photon form factors.
 * @param[in] g     Generation index for quarks.
 * @return Set of four-fermion coefficients from photon penguin amplitudes
 *         (without overall i; with appropriate embedding).
 */
template <class Q, class A, class T>
Eigen::Array<std::complex<double>,10,1> embed_photon(
   const NUTNMSSM_mass_eigenstates& model, const T& ff, int g) {
   // Get quark-photon couplings (without i, as everywhere):
   const auto qL =  left<Q, A>(model, g);
   const auto qR = right<Q, A>(model, g);
   Eigen::Array<std::complex<double>,10,1> res{};
   // Term from eq. (t.1) will contribute to four-fermion vector coefficients.
   // Minus comes from the form_factors embedding into four-fermion amplitude:
   res[4] = - ff[0] * qL;
   res[5] = - ff[0] * qR;
   res[6] = - ff[1] * qL;
   res[7] = - ff[1] * qR;
   return res;
};

/**
 * @tparam    Name  Function name for npf function.
 * @param[in] model Mass eigenstates.
 * @param[in] in    Generation index for incoming lepton.
 * @param[in] out   Generation index for outgoing lepton.
 * @param[in] g     Generation index for quarks.
 * @return Set of four-fermion coefficients for non-photonic amplitudes
 *         (without overall i; with fixed signs for tensor operators).
 */
template <npf_function Name>
Eigen::Array<std::complex<double>,10,1> fix_tensors_sign(
   const NUTNMSSM_mass_eigenstates& model, int in, int out, int g) {
   const auto npf = Name(model,
         std::array<int,4>{in, g, out, g},
         std::array<Eigen::Vector4d, 0>{});
   Eigen::Array<std::complex<double>,10,1> res(npf.data());
   res[8] = - res[8];
   res[9] = - res[9];
   return res;
};

/**
 * Form factors are defined via the following formula (q = pj - pi > 0; see
 * eq. (3.4) of 1902.06650 (up to electric charge);
 * pi - momenta[going from blob] of outgoing lepton,
 * pj - momenta[going into blob] of incoming lepton):
 *    <pi, q| T exp(i L_NUTNMSSM dx)|pj> =
 *    i * ubari
 *          q^2 gamma_mu (A1_X * P_X)                                  (t.1)
 *          + (zero after embedding term)                              (t.2)
 *          + i mj sigma_munu q^nu (A2_X * P_X)                        (t.3)
 *        uj e^*mu
 * @note Form factors below are ordered as A1_L, A1_R, A2_L, A2_R.
 * @note NPF function returns such result, that G4 = i * npf.
 * @note Match to a low energy effective model with the same covariant derivative
 *       for photon, as in NUTNMSSM, using eq. (t.3).
 *       Minus comes from the form_factors embedding into four-fermion
 *       amplitude, because we use descending ordering for external fermions
 *          <out:q4,l3| T exp(i L_NUTNMSSM dx) |in:q2,l1> =: G4
 * @note For four fermion coefficient matching minus comes from the G4
 *       definition, because
 *          <out:q4,l3| T exp(i L_low dx) |in:q2,l1> = -i * C_XY [3X1]*[4Y2].
 * @tparam    Lepton   Type of a lepton fields.
 * @tparam    Up       Type of a up-quark fields.
 * @tparam    Down     Type of a down-quark fields.
 * @tparam    Photon   Type of a photon field.
 * @tparam    photon   Function name for photon form factors.
 * @tparam    npf_up   Function name for npf function for up-quark contribution.
 * @tparam    npf_down Function name for npf function for down-quark
 *                     contribution.
 * @param[in] in     Generation index for incoming lepton.
 * @param[in] out    Generation index for outgoing lepton.
 * @param[in] model  Mass eigenstates.
 * @param[in] parameters Parameters for the observable calculation.
 * @param[in] qedqcd Reference to low-energy data.
 * @return Set of (l to l conversion) process and Wilson coefficients.
 */
template<class Lepton, class Up, class Down, class Photon,
   ffv_function photon, npf_function npf_up, npf_function npf_down>
Eigen::Array<std::complex<double>, 13, 1> forge_conversion(int in, int out,
   const NUTNMSSM_l_to_l_conversion::Nucleus nucleus,
   const NUTNMSSM_mass_eigenstates& model_,
   const LToLConversion_settings& parameters,
   const softsusy::QedQcd& qedqcd) {

   NUTNMSSM_mass_eigenstates model(model_);
   try {
      run_to_MSUSY(model);
      model.get_physical().clear();
      model.solve_ewsb();
   } catch (const Error& e) {
      ERROR("NUTNMSSM_l_to_l_conversion:" << e.what_detailed());
      return std::numeric_limits<Eigen::Array<std::complex<double>,13,1>>::quiet_NaN();
   }

   context_base context {model};

   const auto form_factors = photon(in, out, model, false);
   const auto photon_u = F4F(Up,   0);
   const auto photon_d = F4F(Down, 0);
   const auto photon_s = F4F(Down, 1);

   const auto npf_u = NPF(npf_up,   0);
   const auto npf_d = NPF(npf_down, 0);
   const auto npf_s = NPF(npf_down, 1);

   // Matching
   const auto DL = - 0.5 * form_factors[2];
   const auto DR = - 0.5 * form_factors[3];
   const auto CXYu = - photon_u - npf_u;
   const auto CXYd = - photon_d - npf_d;
   const auto CXYs = - photon_s - npf_d;

   auto CSLu = ( CXYu[0] + CXYu[1] )/2.;
   auto CSRu = ( CXYu[2] + CXYu[3] )/2.;
   auto CSLd = ( CXYd[0] + CXYd[1] )/2.;
   auto CSRd = ( CXYd[2] + CXYd[3] )/2.;
   auto CSLs = ( CXYs[0] + CXYs[1] )/2.;
   auto CSRs = ( CXYs[2] + CXYs[3] )/2.;

   auto CVLu = ( CXYu[4] + CXYu[5] )/2.;
   auto CVRu = ( CXYu[6] + CXYu[7] )/2.;
   auto CVLd = ( CXYd[4] + CXYd[5] )/2.;
   auto CVRd = ( CXYd[6] + CXYd[7] )/2.;

   auto CTLu = CXYu[8], CTRu = CXYu[9];
   auto CTLd = CXYd[8], CTRd = CXYd[9];
   auto CTLs = CXYs[8], CTRs = CXYs[9];

   // Vector contribution.
   auto gpLV = ( GV(p,u)*CVLu + GV(p,d)*CVLd );
   auto gpRV = ( GV(p,u)*CVRu + GV(p,d)*CVRd );
   auto gnLV = ( GV(n,u)*CVLu + GV(n,d)*CVLd );
   auto gnRV = ( GV(n,u)*CVRu + GV(n,d)*CVRd );

   // Scalar contribution from scalar coefficients.
   Mass<Lepton, Up, Down> m(context, in);
   auto gpLS = ( GS(p,u)*CSLu + GS(p,d)*CSLd + GS(p,s)*CSLs );
   auto gpRS = ( GS(p,u)*CSRu + GS(p,d)*CSRd + GS(p,s)*CSRs );
   auto gnLS = ( GS(n,u)*CSLu + GS(n,d)*CSLd + GS(n,s)*CSLs );
   auto gnRS = ( GS(n,u)*CSRu + GS(n,d)*CSRd + GS(n,s)*CSRs );

   if (parameters.get(LToLConversion_settings::include_tensor_contribution)) {
      gpLS += (m.l/m.p)*(GT(p,u)*CTLu + GT(p,d)*CTLd + GT(p,s)*CTLs);
      gpRS += (m.l/m.p)*(GT(p,u)*CTRu + GT(p,d)*CTRd + GT(p,s)*CTRs);
      gnLS += (m.l/m.n)*(GT(n,u)*CTLu + GT(n,d)*CTLd + GT(n,s)*CTLs);
      gnRS += (m.l/m.n)*(GT(n,u)*CTRu + GT(n,d)*CTRd + GT(n,s)*CTRs);
   }

   if (parameters.get(LToLConversion_settings::include_gluonic_contribution)) {
      // Only scalar contributions is needed - no photon; minus from matching.
      const auto CXYc = - NPF(npf_up,   1);
      const auto CXYb = - NPF(npf_down, 2);

      auto CSLc = ( CXYc[0] + CXYc[1] )/2.;
      auto CSRc = ( CXYc[2] + CXYc[3] )/2.;
      auto CSLb = ( CXYb[0] + CXYb[1] )/2.;
      auto CSRb = ( CXYb[2] + CXYb[3] )/2.;

      const double Ggp = - 8.*Pi/9 * (1 - m.u/m.p * GS(p,u)
                                        - m.d/m.p * GS(p,d)
                                        - m.s/m.p * GS(p,s));
      const double Ggn = - 8.*Pi/9 * (1 - m.u/m.n * GS(n,u)
                                        - m.d/m.n * GS(n,d)
                                        - m.s/m.n * GS(n,s));

      // Scalar contribution from gluonic coefficients.
      gpLS += - 1./(12*Pi) * m.p*Ggp*(CSLc/m.c + CSLb/m.b);
      gpRS += - 1./(12*Pi) * m.p*Ggp*(CSRc/m.c + CSRb/m.b);
      gnLS += - 1./(12*Pi) * m.n*Ggn*(CSLc/m.c + CSLb/m.b);
      gnRS += - 1./(12*Pi) * m.n*Ggn*(CSRc/m.c + CSRb/m.b);
   }

   const auto ff = get_overlap_integrals(nucleus, qedqcd);
   const auto left  = DR*ff.D/4. - gpLV*ff.Vp - gnLV*ff.Vn - gpRS*ff.Sp - gnRS*ff.Sn;
   const auto right = DL*ff.D/4. - gpRV*ff.Vp - gnRV*ff.Vn - gpLS*ff.Sp - gnLS*ff.Sn;

   // eq. (3.31) of 1902.06650 (with some coefficient redefinitions)
   const double conversion_rate = 4*(std::norm(left) + std::norm(right));
   const double capture_rate = get_capture_rate(nucleus);
   Eigen::Array<std::complex<double>,13,1> res;
   res << conversion_rate/capture_rate,
          DL,      DR,
          CXYu[0], CXYu[1], CXYu[2], CXYu[3],
          CXYu[4], CXYu[5], CXYu[6], CXYu[7],
          CXYu[8], CXYu[9];
   return res;
}



} // namespace NUTNMSSM_l_to_l_conversion
} // namespace flexiblesusy
