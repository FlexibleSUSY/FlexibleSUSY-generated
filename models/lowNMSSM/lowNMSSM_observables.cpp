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

// File generated at Mon 5 Mar 2018 18:23:07

#include "lowNMSSM_observables.hpp"
#include "lowNMSSM_mass_eigenstates.hpp"
#include "lowNMSSM_a_muon.hpp"
#include "lowNMSSM_edm.hpp"
#include "lowNMSSM_effective_couplings.hpp"
#include "gm2calc_interface.hpp"
#include "eigen_utils.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"
#include "lowe.h"
#include "physical_input.hpp"

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

const int lowNMSSM_observables::NUMBER_OF_OBSERVABLES;

lowNMSSM_observables::lowNMSSM_observables()

{
}

Eigen::ArrayXd lowNMSSM_observables::get() const
{
   Eigen::ArrayXd vec(1);

   vec(0) = 0.;

   return vec;
}

std::vector<std::string> lowNMSSM_observables::get_names()
{
   std::vector<std::string> names(1);

   names[0] = "no observables defined";

   return names;
}

void lowNMSSM_observables::clear()
{

}

void lowNMSSM_observables::set(const Eigen::ArrayXd& vec)
{

}

lowNMSSM_observables calculate_observables(lowNMSSM_mass_eigenstates& model,
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
         return lowNMSSM_observables();
      }
   }

   return calculate_observables(model_at_scale, qedqcd, physical_input);
}

lowNMSSM_observables calculate_observables(lowNMSSM_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input)
{
   lowNMSSM_observables observables;

   try {
      

   } catch (const Error& e) {
      model.get_problems().flag_thrown(e.what());
   }

   return observables;
}

} // namespace flexiblesusy
