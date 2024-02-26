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


#include "CMSSMSemiAnalytic_observables.hpp"
#include "CMSSMSemiAnalytic_mass_eigenstates.hpp"
#include "CMSSMSemiAnalytic_amm.hpp"
#include "CMSSMSemiAnalytic_edm.hpp"
#include "CMSSMSemiAnalytic_b_to_s_gamma.hpp"
#include "observables/l_to_l_conversion/settings.hpp"
#include "observables/CMSSMSemiAnalytic_br_l_to_3l.hpp"
#include "observables/CMSSMSemiAnalytic_br_l_to_l_gamma.hpp"
#include "observables/CMSSMSemiAnalytic_l_to_l_conversion.hpp"
#include "cxx_qft/CMSSMSemiAnalytic_qft.hpp"
#include "config.h"
#include "eigen_utils.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"
#include "lowe.h"
#include "physical_input.hpp"

#ifdef ENABLE_GM2CALC
#include "gm2calc_interface.hpp"
#endif

#define MODEL model
#define AMM0(p) amm_ ## p
#define AMM1(p,idx) amm_ ## p ## _ ## idx
#define AMMUNCERTAINTY0(p) amm_uncertainty_ ## p
#define AMMUNCERTAINTY1(p,idx) amm_uncertainty_ ## p ## _ ## idx
#define AMUGM2CALC a_muon_gm2calc
#define AMUGM2CALCUNCERTAINTY a_muon_gm2calc_uncertainty
#define DERIVEDPARAMETER(p) model.p()
#define EXTRAPARAMETER(p) model.get_##p()
#define INPUTPARAMETER(p) model.get_input().p
#define MODELPARAMETER(p) model.get_##p()
#define PHASE(p) model.get_##p()
#define LowEnergyConstant(p) Electroweak_constants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model.get_physical().p
#define EDM0(p) edm_ ## p
#define EDM1(p,idx) edm_ ## p ## _ ## idx
#define BSGAMMA b_to_s_gamma

#define ALPHA_EM_MZ qedqcd.displayAlpha(softsusy::ALPHA)
#define ALPHA_EM_0 physical_input.get(Physical_input::alpha_em_0)
#define ALPHA_S_MZ qedqcd.displayAlpha(softsusy::ALPHAS)
#define MHPole physical_input.get(Physical_input::mh_pole)
#define MWPole qedqcd.displayPoleMW()
#define MZPole qedqcd.displayPoleMZ()
#define MU2GeV qedqcd.displayMu2GeV()
#define MS2GeV qedqcd.displayMs2GeV()
#define MTPole qedqcd.displayPoleMt()
#define MD2GeV qedqcd.displayMd2GeV()
#define MCMC qedqcd.displayMcMc()
#define MBMB qedqcd.displayMbMb()
#define Mv1Pole qedqcd.displayNeutrinoPoleMass(1)
#define Mv2Pole qedqcd.displayNeutrinoPoleMass(2)
#define Mv3Pole qedqcd.displayNeutrinoPoleMass(3)
#define MEPole qedqcd.displayPoleMel()
#define MMPole qedqcd.displayPoleMmuon()
#define MTauPole qedqcd.displayPoleMtau()
#define CKMInput qedqcd.get_complex_ckm()

namespace flexiblesusy {

const int CMSSMSemiAnalytic_observables::NUMBER_OF_OBSERVABLES;

CMSSMSemiAnalytic_observables::CMSSMSemiAnalytic_observables()
   : amm_Fe_1(0)

{
}

Eigen::ArrayXd CMSSMSemiAnalytic_observables::get() const
{
   Eigen::ArrayXd vec(CMSSMSemiAnalytic_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = amm_Fe_1;

   return vec;
}

std::vector<std::string> CMSSMSemiAnalytic_observables::get_names()
{
   std::vector<std::string> names(CMSSMSemiAnalytic_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "amm_Fe_1";

   return names;
}

void CMSSMSemiAnalytic_observables::clear()
{
   amm_Fe_1 = 0.;

}

void CMSSMSemiAnalytic_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == CMSSMSemiAnalytic_observables::NUMBER_OF_OBSERVABLES);

   amm_Fe_1 = vec(0);

}

CMSSMSemiAnalytic_observables calculate_observables(const CMSSMSemiAnalytic_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              
                                              const Physical_input& physical_input,
                                              const Spectrum_generator_settings& settings,
                                              double scale)
{
   auto model_at_scale = model;

   if (scale > 0.) {
      try {
         model_at_scale.run_to(scale);
      } catch (const NonPerturbativeRunningError& e) {
         CMSSMSemiAnalytic_observables observables;
         observables.problems.general.flag_non_perturbative_running(scale);
         return observables;
      } catch (const Error& e) {
         CMSSMSemiAnalytic_observables observables;
         observables.problems.general.flag_thrown(e.what());
         return observables;
      } catch (const std::exception& e) {
         CMSSMSemiAnalytic_observables observables;
         observables.problems.general.flag_thrown(e.what());
         return observables;
      }
   }

   return calculate_observables(model_at_scale,
                                qedqcd,
                                
                                physical_input,
                                settings);
}

CMSSMSemiAnalytic_observables calculate_observables(const CMSSMSemiAnalytic_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              
                                              const Physical_input& physical_input,
                                              const Spectrum_generator_settings& settings)
{
   CMSSMSemiAnalytic_observables observables;

   try {
      
      observables.AMM1(Fe, 1) = CMSSMSemiAnalytic_amm::calculate_amm<CMSSMSemiAnalytic_cxx_diagrams::fields::Fe>(MODEL, qedqcd, settings,1);
   } catch (const NonPerturbativeRunningError& e) {
      observables.problems.general.flag_non_perturbative_running(e.get_scale());
   } catch (const Error& e) {
      observables.problems.general.flag_thrown(e.what());
   } catch (const std::exception& e) {
      observables.problems.general.flag_thrown(e.what());
   }

   return observables;
}

} // namespace flexiblesusy
