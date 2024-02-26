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
 * @file lowMSSM_br_l_to_3l.cpp
 *
 * This file was generated at Sun 25 Feb 2024 21:03:41 with FlexibleSUSY
 * 2.8.0 and SARAH 4.15.1
 */

#include <valarray>
#include <complex>

#include "lowMSSM_mass_eigenstates.hpp"
#include "cxx_qft/lowMSSM_qft.hpp"

#include "lowMSSM_br_l_to_3l.hpp"
#include "lowMSSM_FFV_form_factors.hpp"

#define NPF(name) fix_tensors_sign<name>(model, nI, nO, nA);
#define F4F(lepton, g) embed_photon<lepton, Photon>(model, form_factors, g);

#include "loop_libraries/loop_library.hpp"
#include "cxx_qft/lowMSSM_npointfunctions_wilsoncoeffs.hpp"
#include "concatenate.hpp"
#include <limits>
#include <type_traits>
#include <boost/fusion/include/at_key.hpp>
#include "wrappers.hpp"

namespace flexiblesusy {

namespace lowMSSM_cxx_diagrams {
namespace npointfunctions {



} // namespace npointfunctions
} // namespace lowMSSM_cxx_diagrams

namespace {

std::valarray<std::complex<double>> zero(
   int generationIndex1,
   int generationIndex2,
   const lowMSSM_mass_eigenstates&,
   bool){
   std::valarray<std::complex<double>> res {0.0, 0.0, 0.0, 0.0};
   return res;
}

std::array<std::complex<double>,10> zero(
   const lowMSSM_mass_eigenstates&,
   const std::array<int,4>&,
   const std::array<Eigen::Vector4d,0>&){
   std::array<std::complex<double>,10> res{};
   return res;
}

double get_MSUSY(const lowMSSM_mass_eigenstates_interface& model)
{
   return Min(model.get_MSd().tail<6>().minCoeff(), model.get_MSu().tail<6>().
      minCoeff(), model.get_MSe().tail<6>().minCoeff(), model.get_MHpm().tail<1>()
      .minCoeff(), model.get_MCha().tail<2>().minCoeff());
}

void run_to_MSUSY(lowMSSM_mass_eigenstates& model)
{
   const double precision_goal = model.get_precision();
   const int Nmax = static_cast<int>(-std::log10(precision_goal) * 10);
   int it = 0;
   double precision = 0.;
   double MSUSY_old = 0., MSUSY_new = 0.;

   VERBOSE_MSG("lowMSSM_br_l_to_3l: iterate to run to MSUSY ...");

   do {
      MSUSY_new = get_MSUSY(model);

      if (MSUSY_new <= 0.)
         return;

      VERBOSE_MSG("lowMSSM_br_l_to_3l:    running to new MSUSY = " << MSUSY_new);

      model.run_to(MSUSY_new);
      model.calculate_DRbar_masses();

      precision = MaxRelDiff(MSUSY_old, MSUSY_new);
      MSUSY_old = MSUSY_new;
   } while (precision > precision_goal && ++it < Nmax);

   VERBOSE_MSG("lowMSSM_br_l_to_3l: iteration finished. MSUSY = " << MSUSY_new);

   if (precision > precision_goal) {
      WARNING("lowMSSM_br_l_to_3l: run_to_MSUSY(): "
              "Iteration did not converge after " << Nmax
              << " iterations (precision goal: " << precision_goal << ")");
   }
}

} // anonymous namespace

using namespace lowMSSM_cxx_diagrams;
using namespace lowMSSM_FFV_form_factors;

namespace lowMSSM_br_l_to_3l {

typedef std::valarray<std::complex<double>> (*ffv)(
   int, int, const lowMSSM_mass_eigenstates&, bool);

typedef std::array<std::complex<double>,10> (*npf)(
   const lowMSSM_mass_eigenstates&,
   const std::array<int,4>&,
   const std::array<Eigen::Vector4d,0>&);

/**
 * @param[in] nI     Generation index of incoming lepton.
 * @return Total decay width in GeV, according to PDG data.
 */
double get_total_width(int nI) {
   // https://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf
   constexpr double hbar = 6.582119569e-25, // [GeV*s]
   // https://pdg.lbl.gov/2020/listings/rpp2020-list-muon.pdf
                    muon = 2.1969811e-6,    // [s]
   // https://pdg.lbl.gov/2020/listings/rpp2020-list-tau.pdf
                    tau  = 290.3e-15;       // [s]
   switch (nI) {
      case 1: return hbar/muon;
      case 2: return hbar/tau;
      default: throw std::invalid_argument("Unrecognized lepton");
   }
}

/**
 * @tparam    Lepton Type of a lepton field.
 * @tparam    Photon Type of a photon field.
 * @param[in] model  Mass eigenstates.
 * @param[in] nA     Generation index.
 * @return Left part of ffv coupling, multiplied by -i.
 */
template <class Lepton, class Photon>
std::complex<double> left(const lowMSSM_mass_eigenstates& model, int nA) {
    context_base context {model};
    using vertex = Vertex<typename Lepton::lorentz_conjugate, Lepton, Photon>;
    std::array<int, 2> indices {nA, nA};
    const auto value =  vertex::evaluate(indices, context);
    return value.left();
}

/**
 * @tparam    Lepton Type of a lepton field.
 * @tparam    Photon Type of a photon field.
 * @param[in] model  Mass eigenstates.
 * @param[in] nA     Generation index.
 * @return Right part of ffv coupling, multiplied by -i.
 */
template <class Lepton, class Photon>
std::complex<double> right(const lowMSSM_mass_eigenstates& model, int nA) {
    context_base context {model};
    using vertex = Vertex<typename Lepton::lorentz_conjugate, Lepton, Photon>;
    std::array<int, 2> indices {nA, nA};
    const auto value =  vertex::evaluate(indices, context);
    return value.right();
}

/**
 * @tparam    L     Type of a lepton field.
 * @tparam    A     Type of a photon field.
 * @tparam    T     Type of a form factors.
 * @param[in] model Mass eigenstates.
 * @param[in] ff    Lepton-photon form factors.
 * @param[in] g     Generation index for leptons.
 * @return Set of four-fermion coefficients from photon penguin amplitudes
 *         (without overall i; with appropriate embedding).
 */
template <class L, class A, class T>
Eigen::Array<std::complex<double>,10,1> embed_photon(
   const lowMSSM_mass_eigenstates& model, const T& ff, int g) {
   // Get quark-photon couplings (without i, as everywhere):
   const auto lL =  left<L, A>(model, g);
   const auto lR = right<L, A>(model, g);
   Eigen::Array<std::complex<double>,10,1> res{};
   // Minus comes from the form_factors embedding into four-fermion amplitude:
   res[4] = - ff[0] * lL;
   res[5] = - ff[0] * lR;
   res[6] = - ff[1] * lL;
   res[7] = - ff[1] * lR;
   return res;
};

/**
 * @tparam    Name  Function name for npf function.
 * @param[in] model Mass eigenstates.
 * @param[in] nI    Generation index for incoming lepton.
 * @param[in] nO    Generation index for outgoing lepton.
 * @param[in] nA    Generation index for lepton pair.
 * @return Set of four-fermion coefficients for non-photonic amplitudes
 *         (without overall i; with fixed signs for tensor operators).
 */
template <npf Name>
Eigen::Array<std::complex<double>,10,1> fix_tensors_sign(
   const lowMSSM_mass_eigenstates& model, int nI, int nO, int nA) {
   const auto npf = Name(model,
         std::array<int,4>{nI, nA, nO, nA},
         std::array<Eigen::Vector4d, 0>{});
   Eigen::Array<std::complex<double>,10,1> res(npf.data());
   res[8] = - res[8];
   res[9] = - res[9];
   return res;
};

/**
 * @tparam    A      Type of input sets of coefficients.
 * @param[in] photon T-coefficients of all penguins.
 * @param[in] res    T-coefficients of all penguins.
 * @param[in] boxes  Coefficients of all boxes.
 * @return Wilson coefficients of four fermion operators.
 */
template<class A>
Eigen::Array<std::complex<double>,10,1> fierz(A photon, A rest, A box) {
   A t_channel = photon + rest;
   A tu_channels;
   tu_channels << t_channel[0]/2 - 6*t_channel[8],
                  t_channel[1] - 2*t_channel[5],
                  t_channel[2] - 2*t_channel[6],
                  t_channel[3]/2 - 6*t_channel[9],
                  2*t_channel[4],
                  t_channel[5] - t_channel[1]/2,
                  t_channel[6] - t_channel[2]/2,
                  2*t_channel[7],
                  1.5*t_channel[8] - t_channel[0]/8,
                  1.5*t_channel[9] - t_channel[3]/8;
   tu_channels = tu_channels + box;
   A res{};
   res << 2*tu_channels[0],
          0,
          0,
          2*tu_channels[3],
          tu_channels[4]/2,
          tu_channels[5],
          tu_channels[6],
          tu_channels[7]/2,
          0,
          0;
   return res;
}

struct LeptonOperators {
   std::complex<double> SLL, SLR, SRL, SRR,
                        VLL, VLR, VRL, VRR,
                        TLL, TRR;
   template<class T>
   LeptonOperators(T coeffs) {
      SLL = coeffs[0];
      SLR = coeffs[1];
      SRL = coeffs[2];
      SRR = coeffs[3];
      VLL = coeffs[4];
      VLR = coeffs[5];
      VRL = coeffs[6];
      VRR = coeffs[7];
      TLL = coeffs[8];
      TRR = coeffs[9];
   }
};

/**
 * @param[in] g Generation index for leptons.
 * @param[in] qedqcd Reference to low-energy data.
 * @return A mass of a lepton for a given generation number.
 */
double get_pole_mass(int g, const softsusy::QedQcd& qedqcd) {
   switch (g) {
      case 0: return qedqcd.displayMass(softsusy::mElectron);
      case 1: return qedqcd.displayMass(softsusy::mMuon);
      case 2: return qedqcd.displayMass(softsusy::mTau);
      default: throw std::invalid_argument("Unrecognized lepton");
   }
}

/**
 * @tparam    T1     Type of a dipole coefficient.
 * @tparam    T2     Type of a four-lepton coefficients.
 * @param[in] nI     Generation of decaying lepton.
 * @param[in] nA     Generation of any lepton in pair.
 * @param[in] model  Mass eigenstates.
 * @param[in] qedqcd Reference to low-energy data.
 * @return Width of decay with leptons of the same generation in the end.
 */
template<class T1, class T2>
double width_same(int nI, int nA, T1 DL, T1 DR, T2 coeffs,
   const lowMSSM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd) {

   context_base context {model};
   const double e = unit_charge(context);
   const auto mass_heavy = get_pole_mass(nI, qedqcd);
   const auto mass_light = get_pole_mass(nA, qedqcd);
   const auto r = mass_heavy / mass_light;
   const LeptonOperators C(coeffs);

   using std::norm;
   using std::real;
   using std::conj;
   using std::log;

   // 1802.06803: wrong eq. (3.12) for X term;
   // 1702.03020: wrong X term for their covariant derivative close to eq. (A.1).
   // hep-ph/0004025, hep-ph/9909265: same;
   double res = e*e * (norm(DL) + norm(DR)) * (-11. + 8.*log(r)) / 192.;
   res += e * (2.*real(C.VLL*conj(DR)) + real(C.VLR*conj(DR))) / 192.;
   res += e * (2.*real(C.VRR*conj(DL)) + real(C.VRL*conj(DL))) / 192.;
   res += (norm(C.SLL) + 16.*norm(C.VLL) + 8.*norm(C.VLR)) / 12288.;
   res += (norm(C.SRR) + 16.*norm(C.VRR) + 8.*norm(C.VRL)) / 12288.;
   return res * pow(mass_heavy, 5) / pow(Pi, 3);
}

/**
 * @tparam    T1     Type of a dipole coefficient.
 * @tparam    T2     Type of a four-lepton coefficients.
 * @param[in] nI     Generation of decaying lepton.
 * @param[in] nA     Generation of any lepton in pair.
 * @param[in] model  Mass eigenstates.
 * @param[in] qedqcd Reference to low-energy data.
 * @return Width of decay with leptons of different generations in the end.
 */
template<class T1, class T2>
double width_diff(int nI, int nA, T1 DL, T1 DR, T2 coeffs,
   const lowMSSM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd) {

   context_base context {model};
   const double e = unit_charge(context);
   const auto mass_heavy = get_pole_mass(nI, qedqcd);
   const auto mass_light = get_pole_mass(nA, qedqcd);
   const auto r = mass_heavy / mass_light;
   const LeptonOperators C(coeffs);

   using std::norm;
   using std::real;
   using std::conj;
   using std::log;

   // 1802.06803: wrong eq. (3.12) - they take heaviest mass of final states.
   //             correct eq. (3.11);
   // 1212.5939: wrong r in eq. (3.10) see eq. (4.9) of hep-ph/9403398, where
   //            text explanation is given;
   double res = e*e * (norm(DL) + norm(DR)) * (-3. + 2.*log(r)) / 48.;
   res += e * (real(C.VLL*conj(DR)) + real(C.VLR*conj(DR))) / 192.;
   res += e * (real(C.VRR*conj(DL)) + real(C.VRL*conj(DL))) / 192.;
   res += (norm(C.SLL) + norm(C.SLR) + norm(C.SRL) + norm(C.SRR)) / 6144.;
   res += (norm(C.VLL) + norm(C.VLR) + norm(C.VRL) + norm(C.VRR)) / 1536.;
   res += (norm(C.TLL) + norm(C.TRR)) / 128.;
   return res * pow(mass_heavy, 5) / pow(Pi, 3);
}

/**
 * @tparam    Lepton  Type of a lepton field.
 * @tparam    Photon  Type of a photon field.
 * @tparam    Factor  Function name for photon t-penguins.
 * @tparam    Scalars Function name for scalar t-penguins.
 * @tparam    Vectors Function name for vector t-penguins.
 * @tparam    Boxes   Function name for all box diagrams.
 * @param[in] nI      Generation of decaying lepton.
 * @param[in] nO      Generation of lepton, which can be separated from pair.
 * @param[in] nA      Generation of lepton and antilepton.
 * @param[in] model   Mass eigenstates.
 * @param[in] qedqcd  Reference to low-energy data.
 * @return Observable value and Wilson coefficients used to derive it.
 */
template<class Lepton, class Photon, ffv Factor, npf Scalars, npf Vectors, npf Boxes>
Eigen::Array<std::complex<double>,13,1> forge(int nI, int nO, int nA,
   const lowMSSM_mass_eigenstates& model_, const softsusy::QedQcd& qedqcd) {

   lowMSSM_mass_eigenstates model(model_);
   try {
      run_to_MSUSY(model);
      model.get_physical().clear();
      model.solve_ewsb();
   } catch (const Error& e) {
      ERROR("lowMSSM_br_l_to_3l:" << e.what_detailed());
      return std::numeric_limits<Eigen::Array<std::complex<double>,13,1>>::quiet_NaN();
   }

   context_base context {model};

   const auto form_factors = Factor(nI, nO, model, false);
   const auto photon_amp = F4F(Lepton, nA);

   // Full amplitude is calculated with the following convention (for the case
   // 4!=3, otherwise Fierz transformations should be used hep-ph/0412245, note
   // the half for every sigma-sigma summation in eqs. 40 and 42):
   //    <out:4,3| T exp(i L_full dx) |in:2,1> = i * npf,
   // which means, that it has to be matched with -1 * C_XY, because
   //    <out:4,3| T exp(i L_low dx) |in:2,1> = -i * C_XY [3 X 1]*[4 Y 2].
   auto t_amp = NPF(Scalars);
   t_amp = t_amp + NPF(Vectors);
   const auto box_amp = NPF(Boxes);

   // If final leptons are the same, then we need to add u-penguins as well. If
   // not, then we can just sum them (we always calculate all boxes in npf).
   auto CXYl = (nO == nA) ? fierz(photon_amp, t_amp, box_amp)
                          : photon_amp + t_amp + box_amp;
   // Matching
   const auto DL = - 0.5 * form_factors[2];
   const auto DR = - 0.5 * form_factors[3];
   CXYl = - CXYl;

   // @todo Running of coefficients

   double partial_width;
   if (nO == nA) partial_width = width_same(nI, nA, DL, DR, CXYl, model, qedqcd);
   else          partial_width = width_diff(nI, nA, DL, DR, CXYl, model, qedqcd);

   const double total_width = get_total_width(nI);
   Eigen::Array<std::complex<double>,13,1> res;

   res << partial_width/total_width,
          DL,
          DR,
          CXYl[0],
          CXYl[1],
          CXYl[2],
          CXYl[3],
          CXYl[4],
          CXYl[5],
          CXYl[6],
          CXYl[7],
          CXYl[8],
          CXYl[9];

   return res;
}



} // namespace lowMSSM_br_l_to_3l
} // namespace flexiblesusy
