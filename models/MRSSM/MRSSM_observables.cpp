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

// File generated at Sun 26 Aug 2018 15:32:52

#include "MRSSM_observables.hpp"
#include "MRSSM_mass_eigenstates.hpp"
#include "MRSSM_a_muon.hpp"
#include "MRSSM_edm.hpp"
#include "MRSSM_effective_couplings.hpp"
#include "config.h"
#include "eigen_utils.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"
#include "lowe.h"
#include "physical_input.hpp"

#ifdef ENABLE_GM2Calc
#include "gm2calc_interface.hpp"
#endif

#define MODEL model
#define AMU a_muon
#define AMUUNCERTAINTY a_muon_uncertainty
#define AMUGM2CALC a_muon_gm2calc
#define AMUGM2CALCUNCERTAINTY a_muon_gm2calc_uncertainty
#define EDM0(p) edm_ ## p
#define EDM1(p,idx) edm_ ## p ## _ ## idx
#define EFFCPHIGGSPHOTONPHOTON eff_cp_higgs_photon_photon
#define EFFCPHIGGSGLUONGLUON eff_cp_higgs_gluon_gluon
#define EFFCPPSEUDOSCALARPHOTONPHOTON eff_cp_pseudoscalar_photon_photon
#define EFFCPPSEUDOSCALARGLUONGLUON eff_cp_pseudoscalar_gluon_gluon

#define ALPHA_S_MZ qedqcd.displayAlpha(softsusy::ALPHAS)
#define MWPole qedqcd.displayPoleMW()
#define MZPole qedqcd.displayPoleMZ()
#define MTPole qedqcd.displayPoleMt()
#define MBMB qedqcd.displayMbMb()
#define MTauPole qedqcd.displayPoleMtau()
#define MMPole qedqcd.displayPoleMmuon()

namespace flexiblesusy {

const int MRSSM_observables::NUMBER_OF_OBSERVABLES;

MRSSM_observables::MRSSM_observables()
   : a_muon(0)
   , eff_cp_higgs_photon_photon(Eigen::Array<std::complex<double>,4,1>::Zero())
   , eff_cp_higgs_gluon_gluon(Eigen::Array<std::complex<double>,4,1>::Zero())
   , eff_cp_pseudoscalar_photon_photon(Eigen::Array<std::complex<double>,3,1>::
      Zero())
   , eff_cp_pseudoscalar_gluon_gluon(Eigen::Array<std::complex<double>,3,1>::Zero(
      ))

{
}

Eigen::ArrayXd MRSSM_observables::get() const
{
   Eigen::ArrayXd vec(MRSSM_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = a_muon;
   vec(1) = Re(eff_cp_higgs_photon_photon(0));
   vec(2) = Im(eff_cp_higgs_photon_photon(0));
   vec(3) = Re(eff_cp_higgs_photon_photon(1));
   vec(4) = Im(eff_cp_higgs_photon_photon(1));
   vec(5) = Re(eff_cp_higgs_photon_photon(2));
   vec(6) = Im(eff_cp_higgs_photon_photon(2));
   vec(7) = Re(eff_cp_higgs_photon_photon(3));
   vec(8) = Im(eff_cp_higgs_photon_photon(3));
   vec(9) = Re(eff_cp_higgs_gluon_gluon(0));
   vec(10) = Im(eff_cp_higgs_gluon_gluon(0));
   vec(11) = Re(eff_cp_higgs_gluon_gluon(1));
   vec(12) = Im(eff_cp_higgs_gluon_gluon(1));
   vec(13) = Re(eff_cp_higgs_gluon_gluon(2));
   vec(14) = Im(eff_cp_higgs_gluon_gluon(2));
   vec(15) = Re(eff_cp_higgs_gluon_gluon(3));
   vec(16) = Im(eff_cp_higgs_gluon_gluon(3));
   vec(17) = Re(eff_cp_pseudoscalar_photon_photon(0));
   vec(18) = Im(eff_cp_pseudoscalar_photon_photon(0));
   vec(19) = Re(eff_cp_pseudoscalar_photon_photon(1));
   vec(20) = Im(eff_cp_pseudoscalar_photon_photon(1));
   vec(21) = Re(eff_cp_pseudoscalar_photon_photon(2));
   vec(22) = Im(eff_cp_pseudoscalar_photon_photon(2));
   vec(23) = Re(eff_cp_pseudoscalar_gluon_gluon(0));
   vec(24) = Im(eff_cp_pseudoscalar_gluon_gluon(0));
   vec(25) = Re(eff_cp_pseudoscalar_gluon_gluon(1));
   vec(26) = Im(eff_cp_pseudoscalar_gluon_gluon(1));
   vec(27) = Re(eff_cp_pseudoscalar_gluon_gluon(2));
   vec(28) = Im(eff_cp_pseudoscalar_gluon_gluon(2));

   return vec;
}

std::vector<std::string> MRSSM_observables::get_names()
{
   std::vector<std::string> names(MRSSM_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "a_muon";
   names[1] = "Re(eff_cp_higgs_photon_photon(0))";
   names[2] = "Im(eff_cp_higgs_photon_photon(0))";
   names[3] = "Re(eff_cp_higgs_photon_photon(1))";
   names[4] = "Im(eff_cp_higgs_photon_photon(1))";
   names[5] = "Re(eff_cp_higgs_photon_photon(2))";
   names[6] = "Im(eff_cp_higgs_photon_photon(2))";
   names[7] = "Re(eff_cp_higgs_photon_photon(3))";
   names[8] = "Im(eff_cp_higgs_photon_photon(3))";
   names[9] = "Re(eff_cp_higgs_gluon_gluon(0))";
   names[10] = "Im(eff_cp_higgs_gluon_gluon(0))";
   names[11] = "Re(eff_cp_higgs_gluon_gluon(1))";
   names[12] = "Im(eff_cp_higgs_gluon_gluon(1))";
   names[13] = "Re(eff_cp_higgs_gluon_gluon(2))";
   names[14] = "Im(eff_cp_higgs_gluon_gluon(2))";
   names[15] = "Re(eff_cp_higgs_gluon_gluon(3))";
   names[16] = "Im(eff_cp_higgs_gluon_gluon(3))";
   names[17] = "Re(eff_cp_pseudoscalar_photon_photon(0))";
   names[18] = "Im(eff_cp_pseudoscalar_photon_photon(0))";
   names[19] = "Re(eff_cp_pseudoscalar_photon_photon(1))";
   names[20] = "Im(eff_cp_pseudoscalar_photon_photon(1))";
   names[21] = "Re(eff_cp_pseudoscalar_photon_photon(2))";
   names[22] = "Im(eff_cp_pseudoscalar_photon_photon(2))";
   names[23] = "Re(eff_cp_pseudoscalar_gluon_gluon(0))";
   names[24] = "Im(eff_cp_pseudoscalar_gluon_gluon(0))";
   names[25] = "Re(eff_cp_pseudoscalar_gluon_gluon(1))";
   names[26] = "Im(eff_cp_pseudoscalar_gluon_gluon(1))";
   names[27] = "Re(eff_cp_pseudoscalar_gluon_gluon(2))";
   names[28] = "Im(eff_cp_pseudoscalar_gluon_gluon(2))";

   return names;
}

void MRSSM_observables::clear()
{
   a_muon = 0.;
   eff_cp_higgs_photon_photon = Eigen::Array<std::complex<double>,4,1>::Zero();
   eff_cp_higgs_gluon_gluon = Eigen::Array<std::complex<double>,4,1>::Zero();
   eff_cp_pseudoscalar_photon_photon = Eigen::Array<std::complex<double>,3,1>::Zero();
   eff_cp_pseudoscalar_gluon_gluon = Eigen::Array<std::complex<double>,3,1>::Zero();

}

void MRSSM_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == MRSSM_observables::NUMBER_OF_OBSERVABLES);

   a_muon = vec(0);
   eff_cp_higgs_photon_photon(0) = std::complex<double>(vec(1), vec(2));
   eff_cp_higgs_photon_photon(1) = std::complex<double>(vec(3), vec(4));
   eff_cp_higgs_photon_photon(2) = std::complex<double>(vec(5), vec(6));
   eff_cp_higgs_photon_photon(3) = std::complex<double>(vec(7), vec(8));
   eff_cp_higgs_gluon_gluon(0) = std::complex<double>(vec(9), vec(10));
   eff_cp_higgs_gluon_gluon(1) = std::complex<double>(vec(11), vec(12));
   eff_cp_higgs_gluon_gluon(2) = std::complex<double>(vec(13), vec(14));
   eff_cp_higgs_gluon_gluon(3) = std::complex<double>(vec(15), vec(16));
   eff_cp_pseudoscalar_photon_photon(0) = std::complex<double>(vec(17), vec(18));
   eff_cp_pseudoscalar_photon_photon(1) = std::complex<double>(vec(19), vec(20));
   eff_cp_pseudoscalar_photon_photon(2) = std::complex<double>(vec(21), vec(22));
   eff_cp_pseudoscalar_gluon_gluon(0) = std::complex<double>(vec(23), vec(24));
   eff_cp_pseudoscalar_gluon_gluon(1) = std::complex<double>(vec(25), vec(26));
   eff_cp_pseudoscalar_gluon_gluon(2) = std::complex<double>(vec(27), vec(28));

}

MRSSM_observables calculate_observables(MRSSM_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input,
                                              double scale)
{
   auto model_at_scale = model;

   if (scale > 0.) {
      try {
         model_at_scale.run_to(scale);
      } catch (const Error& e) {
         model.get_problems().flag_thrown(e.what());
         return MRSSM_observables();
      }
   }

   return calculate_observables(model_at_scale, qedqcd, physical_input);
}

MRSSM_observables calculate_observables(MRSSM_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input)
{
   MRSSM_observables observables;

   try {
      MRSSM_effective_couplings effective_couplings(model, qedqcd, physical_input);
      effective_couplings.calculate_effective_couplings();

      observables.AMU = MRSSM_a_muon::calculate_a_muon(MODEL);
      observables.EFFCPHIGGSPHOTONPHOTON(0) = effective_couplings.get_eff_CphhVPVP(0);
      observables.EFFCPHIGGSPHOTONPHOTON(1) = effective_couplings.get_eff_CphhVPVP(1);
      observables.EFFCPHIGGSPHOTONPHOTON(2) = effective_couplings.get_eff_CphhVPVP(2);
      observables.EFFCPHIGGSPHOTONPHOTON(3) = effective_couplings.get_eff_CphhVPVP(3);
      observables.EFFCPHIGGSGLUONGLUON(0) = effective_couplings.get_eff_CphhVGVG(0);
      observables.EFFCPHIGGSGLUONGLUON(1) = effective_couplings.get_eff_CphhVGVG(1);
      observables.EFFCPHIGGSGLUONGLUON(2) = effective_couplings.get_eff_CphhVGVG(2);
      observables.EFFCPHIGGSGLUONGLUON(3) = effective_couplings.get_eff_CphhVGVG(3);
      observables.EFFCPPSEUDOSCALARPHOTONPHOTON(0) = effective_couplings.get_eff_CpAhVPVP(1);
      observables.EFFCPPSEUDOSCALARPHOTONPHOTON(1) = effective_couplings.get_eff_CpAhVPVP(2);
      observables.EFFCPPSEUDOSCALARPHOTONPHOTON(2) = effective_couplings.get_eff_CpAhVPVP(3);
      observables.EFFCPPSEUDOSCALARGLUONGLUON(0) = effective_couplings.get_eff_CpAhVGVG(1);
      observables.EFFCPPSEUDOSCALARGLUONGLUON(1) = effective_couplings.get_eff_CpAhVGVG(2);
      observables.EFFCPPSEUDOSCALARGLUONGLUON(2) = effective_couplings.get_eff_CpAhVGVG(3);
   } catch (const Error& e) {
      model.get_problems().flag_thrown(e.what());
   }

   return observables;
}

} // namespace flexiblesusy
