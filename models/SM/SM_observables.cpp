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

// File generated at Tue 5 Sep 2017 10:42:15

#include "SM_observables.hpp"
#include "SM_mass_eigenstates.hpp"
#include "SM_effective_couplings.hpp"
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

const unsigned SM_observables::NUMBER_OF_OBSERVABLES;

SM_observables::SM_observables()
   : eff_cp_higgs_photon_photon(0)
   , eff_cp_higgs_gluon_gluon(0)

{
}

Eigen::ArrayXd SM_observables::get() const
{
   Eigen::ArrayXd vec(SM_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = Re(eff_cp_higgs_photon_photon);
   vec(1) = Im(eff_cp_higgs_photon_photon);
   vec(2) = Re(eff_cp_higgs_gluon_gluon);
   vec(3) = Im(eff_cp_higgs_gluon_gluon);

   return vec;
}

std::vector<std::string> SM_observables::get_names()
{
   std::vector<std::string> names(SM_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "Re(eff_cp_higgs_photon_photon)";
   names[1] = "Im(eff_cp_higgs_photon_photon)";
   names[2] = "Re(eff_cp_higgs_gluon_gluon)";
   names[3] = "Im(eff_cp_higgs_gluon_gluon)";

   return names;
}

void SM_observables::clear()
{
   eff_cp_higgs_photon_photon = std::complex<double>(0.,0.);
   eff_cp_higgs_gluon_gluon = std::complex<double>(0.,0.);

}

void SM_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == SM_observables::NUMBER_OF_OBSERVABLES);

   eff_cp_higgs_photon_photon = std::complex<double>(vec(0), vec(1));
   eff_cp_higgs_gluon_gluon = std::complex<double>(vec(2), vec(3));

}

SM_observables calculate_observables(const SM_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input)
{
   SM_observables observables;

   SM_effective_couplings effective_couplings(model, qedqcd, physical_input);
   effective_couplings.calculate_effective_couplings();

   observables.EFFCPHIGGSPHOTONPHOTON = effective_couplings.get_eff_CphhVPVP();
   observables.EFFCPHIGGSGLUONGLUON = effective_couplings.get_eff_CphhVGVG();

   return observables;
}

} // namespace flexiblesusy
