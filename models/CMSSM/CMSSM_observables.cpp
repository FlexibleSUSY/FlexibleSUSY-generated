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

// File generated at Tue 8 Mar 2016 18:45:02

#include "CMSSM_observables.hpp"
#include "CMSSM_mass_eigenstates.hpp"
#include "CMSSM_effective_couplings.hpp"
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

const unsigned CMSSM_observables::NUMBER_OF_OBSERVABLES;

CMSSM_observables::CMSSM_observables()
   : eff_cp_higgs_photon_photon(Eigen::Array<std::complex<double>,2,1>::Zero())
   , eff_cp_higgs_gluon_gluon(Eigen::Array<std::complex<double>,2,1>::Zero())
   , eff_cp_pseudoscalar_photon_photon(0)
   , eff_cp_pseudoscalar_gluon_gluon(0)

{
}

Eigen::ArrayXd CMSSM_observables::get() const
{
   Eigen::ArrayXd vec(CMSSM_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = Re(eff_cp_higgs_photon_photon(0));
   vec(1) = Im(eff_cp_higgs_photon_photon(0));
   vec(2) = Re(eff_cp_higgs_photon_photon(1));
   vec(3) = Im(eff_cp_higgs_photon_photon(1));
   vec(4) = Re(eff_cp_higgs_gluon_gluon(0));
   vec(5) = Im(eff_cp_higgs_gluon_gluon(0));
   vec(6) = Re(eff_cp_higgs_gluon_gluon(1));
   vec(7) = Im(eff_cp_higgs_gluon_gluon(1));
   vec(8) = Re(eff_cp_pseudoscalar_photon_photon);
   vec(9) = Im(eff_cp_pseudoscalar_photon_photon);
   vec(10) = Re(eff_cp_pseudoscalar_gluon_gluon);
   vec(11) = Im(eff_cp_pseudoscalar_gluon_gluon);

   return vec;
}

std::vector<std::string> CMSSM_observables::get_names()
{
   std::vector<std::string> names(CMSSM_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "Re(eff_cp_higgs_photon_photon(0))";
   names[1] = "Im(eff_cp_higgs_photon_photon(0))";
   names[2] = "Re(eff_cp_higgs_photon_photon(1))";
   names[3] = "Im(eff_cp_higgs_photon_photon(1))";
   names[4] = "Re(eff_cp_higgs_gluon_gluon(0))";
   names[5] = "Im(eff_cp_higgs_gluon_gluon(0))";
   names[6] = "Re(eff_cp_higgs_gluon_gluon(1))";
   names[7] = "Im(eff_cp_higgs_gluon_gluon(1))";
   names[8] = "Re(eff_cp_pseudoscalar_photon_photon)";
   names[9] = "Im(eff_cp_pseudoscalar_photon_photon)";
   names[10] = "Re(eff_cp_pseudoscalar_gluon_gluon)";
   names[11] = "Im(eff_cp_pseudoscalar_gluon_gluon)";

   return names;
}

void CMSSM_observables::clear()
{
   eff_cp_higgs_photon_photon = Eigen::Array<std::complex<double>,2,1>::Zero();
   eff_cp_higgs_gluon_gluon = Eigen::Array<std::complex<double>,2,1>::Zero();
   eff_cp_pseudoscalar_photon_photon = std::complex<double>(0.,0.);
   eff_cp_pseudoscalar_gluon_gluon = std::complex<double>(0.,0.);

}

void CMSSM_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == CMSSM_observables::NUMBER_OF_OBSERVABLES);

   eff_cp_higgs_photon_photon(0) = std::complex<double>(vec(0), vec(1));
   eff_cp_higgs_photon_photon(1) = std::complex<double>(vec(2), vec(3));
   eff_cp_higgs_gluon_gluon(0) = std::complex<double>(vec(4), vec(5));
   eff_cp_higgs_gluon_gluon(1) = std::complex<double>(vec(6), vec(7));
   eff_cp_pseudoscalar_photon_photon = std::complex<double>(vec(8), vec(9));
   eff_cp_pseudoscalar_gluon_gluon = std::complex<double>(vec(10), vec(11));

}

CMSSM_observables calculate_observables(const CMSSM_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input)
{
   CMSSM_observables observables;

   CMSSM_effective_couplings effective_couplings(model, qedqcd, physical_input);
   effective_couplings.calculate_effective_couplings();

   observables.EFFCPHIGGSPHOTONPHOTON(0) = effective_couplings.get_eff_CphhVPVP(0);
   observables.EFFCPHIGGSPHOTONPHOTON(1) = effective_couplings.get_eff_CphhVPVP(1);
   observables.EFFCPHIGGSGLUONGLUON(0) = effective_couplings.get_eff_CphhVGVG(0);
   observables.EFFCPHIGGSGLUONGLUON(1) = effective_couplings.get_eff_CphhVGVG(1);
   observables.EFFCPPSEUDOSCALARPHOTONPHOTON = effective_couplings.get_eff_CpAhVPVP(1);
   observables.EFFCPPSEUDOSCALARGLUONGLUON = effective_couplings.get_eff_CpAhVGVG(1);

   return observables;
}

} // namespace flexiblesusy
