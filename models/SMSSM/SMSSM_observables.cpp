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

// File generated at Mon 9 May 2016 13:03:57

#include "SMSSM_observables.hpp"
#include "SMSSM_mass_eigenstates.hpp"
#include "SMSSM_effective_couplings.hpp"
#include "gm2calc_interface.hpp"
#include "eigen_utils.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"
#include "lowe.h"
#include "physical_input.hpp"

#define MODEL model
#define AMUGM2CALC a_muon_gm2calc
#define AMUGM2CALCUNCERTAINTY a_muon_gm2calc_uncertainty
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

const unsigned SMSSM_observables::NUMBER_OF_OBSERVABLES;

SMSSM_observables::SMSSM_observables()
   : eff_cp_higgs_photon_photon(Eigen::Array<std::complex<double>,3,1>::Zero())
   , eff_cp_higgs_gluon_gluon(Eigen::Array<std::complex<double>,3,1>::Zero())
   , eff_cp_pseudoscalar_photon_photon(Eigen::Array<std::complex<double>,2,1>
      ::Zero())
   , eff_cp_pseudoscalar_gluon_gluon(Eigen::Array<std::complex<double>,2,1>
      ::Zero())

{
}

Eigen::ArrayXd SMSSM_observables::get() const
{
   Eigen::ArrayXd vec(SMSSM_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = Re(eff_cp_higgs_photon_photon(0));
   vec(1) = Im(eff_cp_higgs_photon_photon(0));
   vec(2) = Re(eff_cp_higgs_photon_photon(1));
   vec(3) = Im(eff_cp_higgs_photon_photon(1));
   vec(4) = Re(eff_cp_higgs_photon_photon(2));
   vec(5) = Im(eff_cp_higgs_photon_photon(2));
   vec(6) = Re(eff_cp_higgs_gluon_gluon(0));
   vec(7) = Im(eff_cp_higgs_gluon_gluon(0));
   vec(8) = Re(eff_cp_higgs_gluon_gluon(1));
   vec(9) = Im(eff_cp_higgs_gluon_gluon(1));
   vec(10) = Re(eff_cp_higgs_gluon_gluon(2));
   vec(11) = Im(eff_cp_higgs_gluon_gluon(2));
   vec(12) = Re(eff_cp_pseudoscalar_photon_photon(0));
   vec(13) = Im(eff_cp_pseudoscalar_photon_photon(0));
   vec(14) = Re(eff_cp_pseudoscalar_photon_photon(1));
   vec(15) = Im(eff_cp_pseudoscalar_photon_photon(1));
   vec(16) = Re(eff_cp_pseudoscalar_gluon_gluon(0));
   vec(17) = Im(eff_cp_pseudoscalar_gluon_gluon(0));
   vec(18) = Re(eff_cp_pseudoscalar_gluon_gluon(1));
   vec(19) = Im(eff_cp_pseudoscalar_gluon_gluon(1));

   return vec;
}

std::vector<std::string> SMSSM_observables::get_names()
{
   std::vector<std::string> names(SMSSM_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "Re(eff_cp_higgs_photon_photon(0))";
   names[1] = "Im(eff_cp_higgs_photon_photon(0))";
   names[2] = "Re(eff_cp_higgs_photon_photon(1))";
   names[3] = "Im(eff_cp_higgs_photon_photon(1))";
   names[4] = "Re(eff_cp_higgs_photon_photon(2))";
   names[5] = "Im(eff_cp_higgs_photon_photon(2))";
   names[6] = "Re(eff_cp_higgs_gluon_gluon(0))";
   names[7] = "Im(eff_cp_higgs_gluon_gluon(0))";
   names[8] = "Re(eff_cp_higgs_gluon_gluon(1))";
   names[9] = "Im(eff_cp_higgs_gluon_gluon(1))";
   names[10] = "Re(eff_cp_higgs_gluon_gluon(2))";
   names[11] = "Im(eff_cp_higgs_gluon_gluon(2))";
   names[12] = "Re(eff_cp_pseudoscalar_photon_photon(0))";
   names[13] = "Im(eff_cp_pseudoscalar_photon_photon(0))";
   names[14] = "Re(eff_cp_pseudoscalar_photon_photon(1))";
   names[15] = "Im(eff_cp_pseudoscalar_photon_photon(1))";
   names[16] = "Re(eff_cp_pseudoscalar_gluon_gluon(0))";
   names[17] = "Im(eff_cp_pseudoscalar_gluon_gluon(0))";
   names[18] = "Re(eff_cp_pseudoscalar_gluon_gluon(1))";
   names[19] = "Im(eff_cp_pseudoscalar_gluon_gluon(1))";

   return names;
}

void SMSSM_observables::clear()
{
   eff_cp_higgs_photon_photon = Eigen::Array<std::complex<double>,3,1>::Zero();
   eff_cp_higgs_gluon_gluon = Eigen::Array<std::complex<double>,3,1>::Zero();
   eff_cp_pseudoscalar_photon_photon = Eigen::Array<std::complex<double>,2,1>::Zero();
   eff_cp_pseudoscalar_gluon_gluon = Eigen::Array<std::complex<double>,2,1>::Zero();

}

void SMSSM_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == SMSSM_observables::NUMBER_OF_OBSERVABLES);

   eff_cp_higgs_photon_photon(0) = std::complex<double>(vec(0), vec(1));
   eff_cp_higgs_photon_photon(1) = std::complex<double>(vec(2), vec(3));
   eff_cp_higgs_photon_photon(2) = std::complex<double>(vec(4), vec(5));
   eff_cp_higgs_gluon_gluon(0) = std::complex<double>(vec(6), vec(7));
   eff_cp_higgs_gluon_gluon(1) = std::complex<double>(vec(8), vec(9));
   eff_cp_higgs_gluon_gluon(2) = std::complex<double>(vec(10), vec(11));
   eff_cp_pseudoscalar_photon_photon(0) = std::complex<double>(vec(12), vec(13));
   eff_cp_pseudoscalar_photon_photon(1) = std::complex<double>(vec(14), vec(15));
   eff_cp_pseudoscalar_gluon_gluon(0) = std::complex<double>(vec(16), vec(17));
   eff_cp_pseudoscalar_gluon_gluon(1) = std::complex<double>(vec(18), vec(19));

}

SMSSM_observables calculate_observables(const SMSSM_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input)
{
   SMSSM_observables observables;

   SMSSM_effective_couplings effective_couplings(model, qedqcd, physical_input);
   effective_couplings.calculate_effective_couplings();

   observables.EFFCPHIGGSPHOTONPHOTON(0) = effective_couplings.get_eff_CphhVPVP(0);
   observables.EFFCPHIGGSPHOTONPHOTON(1) = effective_couplings.get_eff_CphhVPVP(1);
   observables.EFFCPHIGGSPHOTONPHOTON(2) = effective_couplings.get_eff_CphhVPVP(2);
   observables.EFFCPHIGGSGLUONGLUON(0) = effective_couplings.get_eff_CphhVGVG(0);
   observables.EFFCPHIGGSGLUONGLUON(1) = effective_couplings.get_eff_CphhVGVG(1);
   observables.EFFCPHIGGSGLUONGLUON(2) = effective_couplings.get_eff_CphhVGVG(2);
   observables.EFFCPPSEUDOSCALARPHOTONPHOTON(0) = effective_couplings.get_eff_CpAhVPVP(1);
   observables.EFFCPPSEUDOSCALARPHOTONPHOTON(1) = effective_couplings.get_eff_CpAhVPVP(2);
   observables.EFFCPPSEUDOSCALARGLUONGLUON(0) = effective_couplings.get_eff_CpAhVGVG(1);
   observables.EFFCPPSEUDOSCALARGLUONGLUON(1) = effective_couplings.get_eff_CpAhVGVG(2);

   return observables;
}

} // namespace flexiblesusy
