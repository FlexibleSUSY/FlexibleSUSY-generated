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


#include "NUHMSSMNoFVHimalaya_observables.hpp"
#include "NUHMSSMNoFVHimalaya_mass_eigenstates.hpp"
#include "NUHMSSMNoFVHimalaya_a_muon.hpp"
#include "NUHMSSMNoFVHimalaya_edm.hpp"
#include "NUHMSSMNoFVHimalaya_l_to_lgamma.hpp"
#include "NUHMSSMNoFVHimalaya_b_to_s_gamma.hpp"
#include "NUHMSSMNoFVHimalaya_f_to_f_conversion.hpp"
#include "NUHMSSMNoFVHimalaya_effective_couplings.hpp"
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
#define AMU a_muon
#define AMUUNCERTAINTY a_muon_uncertainty
#define AMUGM2CALC a_muon_gm2calc
#define AMUGM2CALCUNCERTAINTY a_muon_gm2calc_uncertainty
#define EDM0(p) edm_ ## p
#define EDM1(p,idx) edm_ ## p ## _ ## idx
#define LToLGamma0(pIn, pOut, spec) pIn ## _to_ ## pOut ## _ ## spec
#define LToLGamma1(pIn,idxIn,pOut,idxOut,spec) pIn ## _to_ ## pOut ## _ ## spec
#define FToFConversion1(pIn,idxIn,pOut,idxOut,nuclei,qedqcd) pIn ## _to_ ## pOut ## _in_ ## nuclei
#define BSGAMMA b_to_s_gamma
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

const int NUHMSSMNoFVHimalaya_observables::NUMBER_OF_OBSERVABLES;

NUHMSSMNoFVHimalaya_observables::NUHMSSMNoFVHimalaya_observables()
   : a_muon(0)
   , a_muon_uncertainty(0)
   , a_muon_gm2calc(0)
   , a_muon_gm2calc_uncertainty(0)

{
}

Eigen::ArrayXd NUHMSSMNoFVHimalaya_observables::get() const
{
   Eigen::ArrayXd vec(NUHMSSMNoFVHimalaya_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = a_muon;
   vec(1) = a_muon_uncertainty;
   vec(2) = a_muon_gm2calc;
   vec(3) = a_muon_gm2calc_uncertainty;

   return vec;
}

std::vector<std::string> NUHMSSMNoFVHimalaya_observables::get_names()
{
   std::vector<std::string> names(NUHMSSMNoFVHimalaya_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "a_muon";
   names[1] = "a_muon_uncertainty";
   names[2] = "a_muon_gm2calc";
   names[3] = "a_muon_gm2calc_uncertainty";

   return names;
}

void NUHMSSMNoFVHimalaya_observables::clear()
{
   a_muon = 0.;
   a_muon_uncertainty = 0.;
   a_muon_gm2calc = 0.;
   a_muon_gm2calc_uncertainty = 0.;

}

void NUHMSSMNoFVHimalaya_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == NUHMSSMNoFVHimalaya_observables::NUMBER_OF_OBSERVABLES);

   a_muon = vec(0);
   a_muon_uncertainty = vec(1);
   a_muon_gm2calc = vec(2);
   a_muon_gm2calc_uncertainty = vec(3);

}

NUHMSSMNoFVHimalaya_observables calculate_observables(const NUHMSSMNoFVHimalaya_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input,
                                              double scale)
{
   auto model_at_scale = model;

   if (scale > 0.) {
      try {
         model_at_scale.run_to(scale);
      } catch (const NonPerturbativeRunningError& e) {
         NUHMSSMNoFVHimalaya_observables observables;
         observables.problems.general.flag_non_perturbative_running(scale);
         return observables;
      } catch (const Error& e) {
         NUHMSSMNoFVHimalaya_observables observables;
         observables.problems.general.flag_thrown(e.what());
         return observables;
      } catch (const std::exception& e) {
         NUHMSSMNoFVHimalaya_observables observables;
         observables.problems.general.flag_thrown(e.what());
         return observables;
      }
   }

   return calculate_observables(model_at_scale, qedqcd, physical_input);
}

NUHMSSMNoFVHimalaya_observables calculate_observables(const NUHMSSMNoFVHimalaya_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input)
{
   NUHMSSMNoFVHimalaya_observables observables;

   try {
      #ifdef ENABLE_GM2CALC
      GM2Calc_data gm2calc_data;
      gm2calc_data.alpha_s_MZ = ALPHA_S_MZ;
      gm2calc_data.MZ    = MZPole;
      if (!is_zero(MODEL.get_physical().MVWm))
         gm2calc_data.MW = MODEL.get_physical().MVWm;
      else if (!is_zero(MWPole))
         gm2calc_data.MW = MWPole;
      gm2calc_data.mb_mb = MBMB;
      gm2calc_data.MT    = MTPole;
      gm2calc_data.MTau  = MTauPole;
      gm2calc_data.MM    = MMPole;
      gm2calc_data.MA0   = MODEL.get_physical().MAh(1);
      gm2calc_data.MSvm  = MODEL.get_physical().MSvmL;
      gm2calc_data.MSm   = MODEL.get_physical().MSm;
      gm2calc_data.MCha  = MODEL.get_physical().MCha;
      gm2calc_data.MChi  = MODEL.get_physical().MChi;
      gm2calc_data.scale = MODEL.get_scale();
      gm2calc_data.TB    = MODEL.get_vu() / MODEL.get_vd();
      gm2calc_data.Mu    = MODEL.get_Mu();
      gm2calc_data.M1    = MODEL.get_MassB();
      gm2calc_data.M2    = MODEL.get_MassWB();
      gm2calc_data.M3    = MODEL.get_MassG();
      gm2calc_data.mq2   = MODEL.get_mq2();
      gm2calc_data.mu2   = MODEL.get_mu2();
      gm2calc_data.md2   = MODEL.get_md2();
      gm2calc_data.ml2   = MODEL.get_ml2();
      gm2calc_data.me2   = MODEL.get_me2();
      gm2calc_data.Au    = div_safe(MODEL.get_TYu(), MODEL.get_Yu());
      gm2calc_data.Ad    = div_safe(MODEL.get_TYd(), MODEL.get_Yd());
      gm2calc_data.Ae    = div_safe(MODEL.get_TYe(), MODEL.get_Ye());
      #endif


      observables.AMU = NUHMSSMNoFVHimalaya_a_muon::calculate_a_muon(MODEL, qedqcd);
      observables.AMUUNCERTAINTY = NUHMSSMNoFVHimalaya_a_muon::calculate_a_muon_uncertainty(MODEL, qedqcd);
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
