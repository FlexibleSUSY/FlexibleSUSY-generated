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


#include "THDMII_observables.hpp"
#include "THDMII_mass_eigenstates.hpp"
#include "THDMII_amm.hpp"
#include "THDMII_edm.hpp"
#include "THDMII_b_to_s_gamma.hpp"
#include "observables/l_to_l_conversion/settings.hpp"
#include "observables/THDMII_br_l_to_3l.hpp"
#include "observables/THDMII_br_l_to_l_gamma.hpp"
#include "observables/THDMII_l_to_l_conversion.hpp"
#include "cxx_qft/THDMII_qft.hpp"
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

const int THDMII_observables::NUMBER_OF_OBSERVABLES;

THDMII_observables::THDMII_observables()
   : amm_Fe_1(0)
   , amm_uncertainty_Fe_1(0)
   , a_muon_gm2calc(0)
   , a_muon_gm2calc_uncertainty(0)

{
}

Eigen::ArrayXd THDMII_observables::get() const
{
   Eigen::ArrayXd vec(THDMII_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = amm_Fe_1;
   vec(1) = amm_uncertainty_Fe_1;
   vec(2) = a_muon_gm2calc;
   vec(3) = a_muon_gm2calc_uncertainty;

   return vec;
}

std::vector<std::string> THDMII_observables::get_names()
{
   std::vector<std::string> names(THDMII_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "amm_Fe_1";
   names[1] = "amm_uncertainty_Fe_1";
   names[2] = "a_muon_gm2calc";
   names[3] = "a_muon_gm2calc_uncertainty";

   return names;
}

void THDMII_observables::clear()
{
   amm_Fe_1 = 0.;
   amm_uncertainty_Fe_1 = 0.;
   a_muon_gm2calc = 0.;
   a_muon_gm2calc_uncertainty = 0.;

}

void THDMII_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == THDMII_observables::NUMBER_OF_OBSERVABLES);

   amm_Fe_1 = vec(0);
   amm_uncertainty_Fe_1 = vec(1);
   a_muon_gm2calc = vec(2);
   a_muon_gm2calc_uncertainty = vec(3);

}

THDMII_observables calculate_observables(const THDMII_mass_eigenstates& model,
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
         THDMII_observables observables;
         observables.problems.general.flag_non_perturbative_running(scale);
         return observables;
      } catch (const Error& e) {
         THDMII_observables observables;
         observables.problems.general.flag_thrown(e.what());
         return observables;
      } catch (const std::exception& e) {
         THDMII_observables observables;
         observables.problems.general.flag_thrown(e.what());
         return observables;
      }
   }

   return calculate_observables(model_at_scale,
                                qedqcd,
                                
                                physical_input,
                                settings);
}

THDMII_observables calculate_observables(const THDMII_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              
                                              const Physical_input& physical_input,
                                              const Spectrum_generator_settings& settings)
{
   THDMII_observables observables;

   try {
      #ifdef ENABLE_GM2CALC
      GM2Calc_THDM_data gm2calc_data;
      {
         const auto Lambda1 = MODELPARAMETER(Lambda1);
         const auto Lambda2 = MODELPARAMETER(Lambda2);
         const auto Lambda3 = MODELPARAMETER(Lambda3);
         const auto Lambda4 = MODELPARAMETER(Lambda4);
         const auto Lambda5 = MODELPARAMETER(Lambda5);
         const auto Lambda6 = MODELPARAMETER(Lambda6);
         const auto Lambda7 = MODELPARAMETER(Lambda7);
         const auto v1 = MODELPARAMETER(v1);
         const auto v2 = MODELPARAMETER(v2);
         const auto M122 = MODELPARAMETER(M122);

         gm2calc_data.alpha_em_mz = ALPHA_EM_MZ;
         gm2calc_data.alpha_em_0 = ALPHA_EM_0;
         gm2calc_data.alpha_s_mz = ALPHA_S_MZ;
         if (!is_zero(MODEL.get_physical().Mhh(0))) {
            gm2calc_data.mh = MODEL.get_physical().Mhh(0);
         } else if (!is_zero(MHPole)) {
            gm2calc_data.mh = MHPole;
         }
         if (!is_zero(MODEL.get_physical().MVWm)) {
            gm2calc_data.mw = MODEL.get_physical().MVWm;
         } else if (!is_zero(MWPole)) {
            gm2calc_data.mw = MWPole;
         }
         gm2calc_data.mz = MZPole;
         gm2calc_data.mu(0) = MU2GeV;
         gm2calc_data.mu(1) = MCMC;
         gm2calc_data.mu(2) = MTPole;
         gm2calc_data.md(0) = MD2GeV;
         gm2calc_data.md(1) = MS2GeV;
         gm2calc_data.md(2) = MBMB;
         gm2calc_data.mv(0) = Mv1Pole;
         gm2calc_data.mv(1) = Mv2Pole;
         gm2calc_data.mv(2) = Mv3Pole;
         gm2calc_data.ml(0) = MEPole;
         gm2calc_data.ml(1) = MMPole;
         gm2calc_data.ml(2) = MTauPole;
         gm2calc_data.ckm = CKMInput;
         gm2calc_data.yukawa_type = 2;
         gm2calc_data.lambda(0) = 2*Lambda1;
         gm2calc_data.lambda(1) = 2*Lambda2;
         gm2calc_data.lambda(2) = Lambda3;
         gm2calc_data.lambda(3) = Lambda4;
         gm2calc_data.lambda(4) = Lambda5;
         gm2calc_data.lambda(5) = Lambda6;
         gm2calc_data.lambda(6) = Lambda7;
         gm2calc_data.tan_beta = v2/v1;
         gm2calc_data.m122 = M122;
         gm2calc_data.zeta_u = 0;
         gm2calc_data.zeta_d = 0;
         gm2calc_data.zeta_l = 0;
         gm2calc_data.delta_u = ZEROMATRIX(3,3);
         gm2calc_data.delta_d = ZEROMATRIX(3,3);
         gm2calc_data.delta_l = ZEROMATRIX(3,3);
         gm2calc_data.pi_u = ZEROMATRIX(3,3);
         gm2calc_data.pi_d = ZEROMATRIX(3,3);
         gm2calc_data.pi_l = ZEROMATRIX(3,3);
      }
      #endif


      observables.AMM1(Fe, 1) = THDMII_amm::calculate_amm<THDMII_cxx_diagrams::fields::Fe>(MODEL, qedqcd, settings,1);
      observables.AMMUNCERTAINTY1(Fe, 1) = THDMII_amm::calculate_amm_uncertainty<THDMII_cxx_diagrams::fields::Fe>(MODEL, qedqcd, settings, 1);
      #ifdef ENABLE_GM2CALC
      observables.AMUGM2CALC = gm2calc_calculate_amu(gm2calc_data);
      #endif
      #ifdef ENABLE_GM2CALC
      observables.AMUGM2CALCUNCERTAINTY = gm2calc_calculate_amu_uncertainty(gm2calc_data);
      #endif
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
