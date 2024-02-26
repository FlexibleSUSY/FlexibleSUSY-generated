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


#include "CMSSMNoFV_observables.hpp"
#include "CMSSMNoFV_mass_eigenstates.hpp"
#include "CMSSMNoFV_amm.hpp"
#include "CMSSMNoFV_edm.hpp"
#include "CMSSMNoFV_b_to_s_gamma.hpp"
#include "observables/l_to_l_conversion/settings.hpp"
#include "observables/CMSSMNoFV_br_l_to_3l.hpp"
#include "observables/CMSSMNoFV_br_l_to_l_gamma.hpp"
#include "observables/CMSSMNoFV_l_to_l_conversion.hpp"
#include "cxx_qft/CMSSMNoFV_qft.hpp"
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

const int CMSSMNoFV_observables::NUMBER_OF_OBSERVABLES;

CMSSMNoFV_observables::CMSSMNoFV_observables()
   : amm_Fm(0)
   , amm_uncertainty_Fm(0)
   , a_muon_gm2calc(0)
   , a_muon_gm2calc_uncertainty(0)

{
}

Eigen::ArrayXd CMSSMNoFV_observables::get() const
{
   Eigen::ArrayXd vec(CMSSMNoFV_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = amm_Fm;
   vec(1) = amm_uncertainty_Fm;
   vec(2) = a_muon_gm2calc;
   vec(3) = a_muon_gm2calc_uncertainty;

   return vec;
}

std::vector<std::string> CMSSMNoFV_observables::get_names()
{
   std::vector<std::string> names(CMSSMNoFV_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "amm_Fm";
   names[1] = "amm_uncertainty_Fm";
   names[2] = "a_muon_gm2calc";
   names[3] = "a_muon_gm2calc_uncertainty";

   return names;
}

void CMSSMNoFV_observables::clear()
{
   amm_Fm = 0.;
   amm_uncertainty_Fm = 0.;
   a_muon_gm2calc = 0.;
   a_muon_gm2calc_uncertainty = 0.;

}

void CMSSMNoFV_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == CMSSMNoFV_observables::NUMBER_OF_OBSERVABLES);

   amm_Fm = vec(0);
   amm_uncertainty_Fm = vec(1);
   a_muon_gm2calc = vec(2);
   a_muon_gm2calc_uncertainty = vec(3);

}

CMSSMNoFV_observables calculate_observables(const CMSSMNoFV_mass_eigenstates& model,
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
         CMSSMNoFV_observables observables;
         observables.problems.general.flag_non_perturbative_running(scale);
         return observables;
      } catch (const Error& e) {
         CMSSMNoFV_observables observables;
         observables.problems.general.flag_thrown(e.what());
         return observables;
      } catch (const std::exception& e) {
         CMSSMNoFV_observables observables;
         observables.problems.general.flag_thrown(e.what());
         return observables;
      }
   }

   return calculate_observables(model_at_scale,
                                qedqcd,
                                
                                physical_input,
                                settings);
}

CMSSMNoFV_observables calculate_observables(const CMSSMNoFV_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              
                                              const Physical_input& physical_input,
                                              const Spectrum_generator_settings& settings)
{
   CMSSMNoFV_observables observables;

   try {
      #ifdef ENABLE_GM2CALC
      GM2Calc_MSSMNoFV_data gm2calc_data;
      gm2calc_data.scale = MODEL.get_scale();
      gm2calc_data.alpha_em_MZ = ALPHA_EM_MZ;
      gm2calc_data.alpha_em_0 = ALPHA_EM_0;
      gm2calc_data.alpha_s_MZ = ALPHA_S_MZ;
      gm2calc_data.MZ    = MZPole;
      if (!is_zero(MODEL.get_physical().MVWm)) {
         gm2calc_data.MW = MODEL.get_physical().MVWm;
      } else if (!is_zero(MWPole)) {
         gm2calc_data.MW = MWPole;
      }
      gm2calc_data.mb_mb = MBMB;
      gm2calc_data.MT    = MTPole;
      gm2calc_data.MTau  = MTauPole;
      gm2calc_data.MM    = MMPole;
      gm2calc_data.MA0   = MODEL.get_physical().MAh(1);
      gm2calc_data.MSvm  = MODEL.get_physical().MSvmL;
      gm2calc_data.TB    = MODEL.get_vu() / MODEL.get_vd();
      gm2calc_data.Mu    = MODEL.get_Mu();
      gm2calc_data.M1    = MODEL.get_MassB();
      gm2calc_data.M2    = MODEL.get_MassWB();
      gm2calc_data.M3    = MODEL.get_MassG();
      gm2calc_data.MSm   = MODEL.get_physical().MSm;
      gm2calc_data.MCha  = MODEL.get_physical().MCha;
      gm2calc_data.MChi  = MODEL.get_physical().MChi;
      gm2calc_data.mq2   = MODEL.get_mq2();
      gm2calc_data.mu2   = MODEL.get_mu2();
      gm2calc_data.md2   = MODEL.get_md2();
      gm2calc_data.ml2   = MODEL.get_ml2();
      gm2calc_data.me2   = MODEL.get_me2();
      gm2calc_data.Au    = div_safe(MODEL.get_TYu(), MODEL.get_Yu());
      gm2calc_data.Ad    = div_safe(MODEL.get_TYd(), MODEL.get_Yd());
      gm2calc_data.Ae    = div_safe(MODEL.get_TYe(), MODEL.get_Ye());
      #endif


      observables.AMM0(Fm) = CMSSMNoFV_amm::calculate_amm<CMSSMNoFV_cxx_diagrams::fields::Fm>(MODEL, qedqcd, settings);
      observables.AMMUNCERTAINTY0(Fm) = CMSSMNoFV_amm::calculate_amm_uncertainty<CMSSMNoFV_cxx_diagrams::fields::Fm>(MODEL, qedqcd, settings);
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
