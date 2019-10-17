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

// File generated at Wed 16 Oct 2019 21:28:52

#include "THDMII_observables.hpp"
#include "THDMII_mass_eigenstates.hpp"
#include "THDMII_a_muon.hpp"
#include "THDMII_edm.hpp"
#include "THDMII_l_to_lgamma.hpp"
//#include "THDMII_f_to_f_conversion.hpp"
#include "THDMII_effective_couplings.hpp"
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

const int THDMII_observables::NUMBER_OF_OBSERVABLES;

THDMII_observables::THDMII_observables()
   : eff_cp_higgs_photon_photon(Eigen::Array<std::complex<double>,2,1>::Zero())
   , eff_cp_higgs_gluon_gluon(Eigen::Array<std::complex<double>,2,1>::Zero())
   , eff_cp_pseudoscalar_photon_photon(0)
   , eff_cp_pseudoscalar_gluon_gluon(0)

{
}

Eigen::ArrayXd THDMII_observables::get() const
{
   Eigen::ArrayXd vec(THDMII_observables::NUMBER_OF_OBSERVABLES);

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

std::vector<std::string> THDMII_observables::get_names()
{
   std::vector<std::string> names(THDMII_observables::NUMBER_OF_OBSERVABLES);

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

void THDMII_observables::clear()
{
   eff_cp_higgs_photon_photon = Eigen::Array<std::complex<double>,2,1>::Zero();
   eff_cp_higgs_gluon_gluon = Eigen::Array<std::complex<double>,2,1>::Zero();
   eff_cp_pseudoscalar_photon_photon = std::complex<double>(0.,0.);
   eff_cp_pseudoscalar_gluon_gluon = std::complex<double>(0.,0.);

}

void THDMII_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == THDMII_observables::NUMBER_OF_OBSERVABLES);

   eff_cp_higgs_photon_photon(0) = std::complex<double>(vec(0), vec(1));
   eff_cp_higgs_photon_photon(1) = std::complex<double>(vec(2), vec(3));
   eff_cp_higgs_gluon_gluon(0) = std::complex<double>(vec(4), vec(5));
   eff_cp_higgs_gluon_gluon(1) = std::complex<double>(vec(6), vec(7));
   eff_cp_pseudoscalar_photon_photon = std::complex<double>(vec(8), vec(9));
   eff_cp_pseudoscalar_gluon_gluon = std::complex<double>(vec(10), vec(11));

}

THDMII_observables calculate_observables(THDMII_mass_eigenstates& model,
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
         return THDMII_observables();
      }
   }

   return calculate_observables(model_at_scale, qedqcd, physical_input);
}

THDMII_observables calculate_observables(THDMII_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input)
{
   THDMII_observables observables;

   try {
      THDMII_effective_couplings effective_couplings(model, qedqcd, physical_input);
      effective_couplings.calculate_effective_couplings();

      observables.EFFCPHIGGSPHOTONPHOTON(0) = effective_couplings.get_eff_CphhVPVP(0);
      observables.EFFCPHIGGSPHOTONPHOTON(1) = effective_couplings.get_eff_CphhVPVP(1);
      observables.EFFCPHIGGSGLUONGLUON(0) = effective_couplings.get_eff_CphhVGVG(0);
      observables.EFFCPHIGGSGLUONGLUON(1) = effective_couplings.get_eff_CphhVGVG(1);
      observables.EFFCPPSEUDOSCALARPHOTONPHOTON = effective_couplings.get_eff_CpAhVPVP(1);
      observables.EFFCPPSEUDOSCALARGLUONGLUON = effective_couplings.get_eff_CpAhVGVG(1);
   } catch (const Error& e) {
      model.get_problems().flag_thrown(e.what_detailed());
   }

   return observables;
}

} // namespace flexiblesusy
