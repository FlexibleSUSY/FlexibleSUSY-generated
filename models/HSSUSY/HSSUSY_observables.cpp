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

// File generated at Wed 16 Oct 2019 21:37:02

#include "HSSUSY_observables.hpp"
#include "HSSUSY_mass_eigenstates.hpp"
#include "HSSUSY_a_muon.hpp"
#include "HSSUSY_edm.hpp"
#include "HSSUSY_l_to_lgamma.hpp"
//#include "HSSUSY_f_to_f_conversion.hpp"
#include "HSSUSY_effective_couplings.hpp"
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
#define LToLGamma0(pIn, pOut, spec) pIn ## _to_ ## pOut ## _ ## spec
#define LToLGamma1(pIn,idxIn,pOut,idxOut,spec) pIn ## _to_ ## pOut ## _ ## spec
#define FToFConversion1(pIn,idxIn,pOut,idxOut,nuclei) pIn ## _to_ ## pOut ## _in_ ## nuclei
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

const int HSSUSY_observables::NUMBER_OF_OBSERVABLES;

HSSUSY_observables::HSSUSY_observables()
   : eff_cp_higgs_photon_photon(0)
   , eff_cp_higgs_gluon_gluon(0)

{
}

Eigen::ArrayXd HSSUSY_observables::get() const
{
   Eigen::ArrayXd vec(HSSUSY_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = Re(eff_cp_higgs_photon_photon);
   vec(1) = Im(eff_cp_higgs_photon_photon);
   vec(2) = Re(eff_cp_higgs_gluon_gluon);
   vec(3) = Im(eff_cp_higgs_gluon_gluon);

   return vec;
}

std::vector<std::string> HSSUSY_observables::get_names()
{
   std::vector<std::string> names(HSSUSY_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "Re(eff_cp_higgs_photon_photon)";
   names[1] = "Im(eff_cp_higgs_photon_photon)";
   names[2] = "Re(eff_cp_higgs_gluon_gluon)";
   names[3] = "Im(eff_cp_higgs_gluon_gluon)";

   return names;
}

void HSSUSY_observables::clear()
{
   eff_cp_higgs_photon_photon = std::complex<double>(0.,0.);
   eff_cp_higgs_gluon_gluon = std::complex<double>(0.,0.);

}

void HSSUSY_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == HSSUSY_observables::NUMBER_OF_OBSERVABLES);

   eff_cp_higgs_photon_photon = std::complex<double>(vec(0), vec(1));
   eff_cp_higgs_gluon_gluon = std::complex<double>(vec(2), vec(3));

}

HSSUSY_observables calculate_observables(HSSUSY_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input,
                                              double scale)
{
   auto model_at_scale = model;

   if (scale > 0.) {
      try {
         model_at_scale.run_to(scale);
      } catch (const Error& e) {
         model.get_problems().flag_thrown(e.what_detailed());
         return HSSUSY_observables();
      }
   }

   return calculate_observables(model_at_scale, qedqcd, physical_input);
}

HSSUSY_observables calculate_observables(HSSUSY_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input)
{
   HSSUSY_observables observables;

   try {
      HSSUSY_effective_couplings effective_couplings(model, qedqcd, physical_input);
      effective_couplings.calculate_effective_couplings();

      observables.EFFCPHIGGSPHOTONPHOTON = effective_couplings.get_eff_CphhVPVP();
      observables.EFFCPHIGGSGLUONGLUON = effective_couplings.get_eff_CphhVGVG();
   } catch (const Error& e) {
      model.get_problems().flag_thrown(e.what_detailed());
   }

   return observables;
}

} // namespace flexiblesusy
