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

// File generated at Sun 26 Aug 2018 14:44:47

#include "NUTNMSSM_observables.hpp"
#include "NUTNMSSM_mass_eigenstates.hpp"
#include "NUTNMSSM_a_muon.hpp"
#include "NUTNMSSM_edm.hpp"
#include "NUTNMSSM_effective_couplings.hpp"
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

const int NUTNMSSM_observables::NUMBER_OF_OBSERVABLES;

NUTNMSSM_observables::NUTNMSSM_observables()
   : a_muon(0)

{
}

Eigen::ArrayXd NUTNMSSM_observables::get() const
{
   Eigen::ArrayXd vec(NUTNMSSM_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = a_muon;

   return vec;
}

std::vector<std::string> NUTNMSSM_observables::get_names()
{
   std::vector<std::string> names(NUTNMSSM_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "a_muon";

   return names;
}

void NUTNMSSM_observables::clear()
{
   a_muon = 0.;

}

void NUTNMSSM_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == NUTNMSSM_observables::NUMBER_OF_OBSERVABLES);

   a_muon = vec(0);

}

NUTNMSSM_observables calculate_observables(NUTNMSSM_mass_eigenstates& model,
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
         return NUTNMSSM_observables();
      }
   }

   return calculate_observables(model_at_scale, qedqcd, physical_input);
}

NUTNMSSM_observables calculate_observables(NUTNMSSM_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input)
{
   NUTNMSSM_observables observables;

   try {
      
      observables.AMU = NUTNMSSM_a_muon::calculate_a_muon(MODEL);
   } catch (const Error& e) {
      model.get_problems().flag_thrown(e.what());
   }

   return observables;
}

} // namespace flexiblesusy
