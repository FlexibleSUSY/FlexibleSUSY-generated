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


#include "MSSM_observables.hpp"
#include "MSSM_mass_eigenstates.hpp"
#include "MSSM_a_muon.hpp"
#include "MSSM_edm.hpp"
#include "MSSM_l_to_lgamma.hpp"
#include "MSSM_b_to_s_gamma.hpp"
#include "MSSM_f_to_f_conversion.hpp"
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

#define ALPHA_S_MZ qedqcd.displayAlpha(softsusy::ALPHAS)
#define MWPole qedqcd.displayPoleMW()
#define MZPole qedqcd.displayPoleMZ()
#define MTPole qedqcd.displayPoleMt()
#define MBMB qedqcd.displayMbMb()
#define MTauPole qedqcd.displayPoleMtau()
#define MMPole qedqcd.displayPoleMmuon()

namespace flexiblesusy {

const int MSSM_observables::NUMBER_OF_OBSERVABLES;

MSSM_observables::MSSM_observables()
   : a_muon(0)
   , edm_Fe_0(0)
   , edm_Fe_1(0)
   , edm_Fe_2(0)
   , Fe_to_Fe_VP(0)

{
}

Eigen::ArrayXd MSSM_observables::get() const
{
   Eigen::ArrayXd vec(MSSM_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = a_muon;
   vec(1) = edm_Fe_0;
   vec(2) = edm_Fe_1;
   vec(3) = edm_Fe_2;
   vec(4) = Fe_to_Fe_VP;

   return vec;
}

std::vector<std::string> MSSM_observables::get_names()
{
   std::vector<std::string> names(MSSM_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "a_muon";
   names[1] = "edm_Fe_0";
   names[2] = "edm_Fe_1";
   names[3] = "edm_Fe_2";
   names[4] = "Fe_to_Fe_VP";

   return names;
}

void MSSM_observables::clear()
{
   a_muon = 0.;
   edm_Fe_0 = 0.;
   edm_Fe_1 = 0.;
   edm_Fe_2 = 0.;
   Fe_to_Fe_VP = 0.;

}

void MSSM_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == MSSM_observables::NUMBER_OF_OBSERVABLES);

   a_muon = vec(0);
   edm_Fe_0 = vec(1);
   edm_Fe_1 = vec(2);
   edm_Fe_2 = vec(3);
   Fe_to_Fe_VP = vec(4);

}

MSSM_observables calculate_observables(const MSSM_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input,
                                              double scale)
{
   auto model_at_scale = model;

   if (scale > 0.) {
      try {
         model_at_scale.run_to(scale);
      } catch (const NonPerturbativeRunningError& e) {
         MSSM_observables observables;
         observables.problems.general.flag_non_perturbative_running(scale);
         return observables;
      } catch (const Error& e) {
         MSSM_observables observables;
         observables.problems.general.flag_thrown(e.what());
         return observables;
      } catch (const std::exception& e) {
         MSSM_observables observables;
         observables.problems.general.flag_thrown(e.what());
         return observables;
      }
   }

   return calculate_observables(model_at_scale, qedqcd, physical_input);
}

MSSM_observables calculate_observables(const MSSM_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input)
{
   MSSM_observables observables;

   try {
      
      observables.AMU = MSSM_a_muon::calculate_a_muon(MODEL, qedqcd);
      observables.EDM1(Fe, 0) = MSSM_edm::calculate_edm_Fe(0, MODEL);
      observables.EDM1(Fe, 1) = MSSM_edm::calculate_edm_Fe(1, MODEL);
      observables.EDM1(Fe, 2) = MSSM_edm::calculate_edm_Fe(2, MODEL);
      observables.LToLGamma1(Fe, 1, Fe, 0, VP) = MSSM_l_to_lgamma::calculate_Fe_to_Fe_VP(1, 0, MODEL, qedqcd, physical_input);
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
