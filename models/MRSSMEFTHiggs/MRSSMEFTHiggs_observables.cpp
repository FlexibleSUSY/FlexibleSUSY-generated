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


#include "MRSSMEFTHiggs_observables.hpp"
#include "MRSSMEFTHiggs_mass_eigenstates.hpp"
#include "MRSSMEFTHiggs_a_muon.hpp"
#include "MRSSMEFTHiggs_edm.hpp"
#include "MRSSMEFTHiggs_l_to_lgamma.hpp"
#include "MRSSMEFTHiggs_b_to_s_gamma.hpp"
#include "MRSSMEFTHiggs_f_to_f_conversion.hpp"
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

const int MRSSMEFTHiggs_observables::NUMBER_OF_OBSERVABLES;

MRSSMEFTHiggs_observables::MRSSMEFTHiggs_observables()
   : a_muon(0)

{
}

Eigen::ArrayXd MRSSMEFTHiggs_observables::get() const
{
   Eigen::ArrayXd vec(MRSSMEFTHiggs_observables::NUMBER_OF_OBSERVABLES);

   vec(0) = a_muon;

   return vec;
}

std::vector<std::string> MRSSMEFTHiggs_observables::get_names()
{
   std::vector<std::string> names(MRSSMEFTHiggs_observables::NUMBER_OF_OBSERVABLES);

   names[0] = "a_muon";

   return names;
}

void MRSSMEFTHiggs_observables::clear()
{
   a_muon = 0.;

}

void MRSSMEFTHiggs_observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == MRSSMEFTHiggs_observables::NUMBER_OF_OBSERVABLES);

   a_muon = vec(0);

}

MRSSMEFTHiggs_observables calculate_observables(const MRSSMEFTHiggs_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input,
                                              double scale)
{
   auto model_at_scale = model;

   if (scale > 0.) {
      try {
         model_at_scale.run_to(scale);
      } catch (const NonPerturbativeRunningError& e) {
         MRSSMEFTHiggs_observables observables;
         observables.problems.general.flag_non_perturbative_running(scale);
         return observables;
      } catch (const Error& e) {
         MRSSMEFTHiggs_observables observables;
         observables.problems.general.flag_thrown(e.what());
         return observables;
      } catch (const std::exception& e) {
         MRSSMEFTHiggs_observables observables;
         observables.problems.general.flag_thrown(e.what());
         return observables;
      }
   }

   return calculate_observables(model_at_scale, qedqcd, physical_input);
}

MRSSMEFTHiggs_observables calculate_observables(const MRSSMEFTHiggs_mass_eigenstates& model,
                                              const softsusy::QedQcd& qedqcd,
                                              const Physical_input& physical_input)
{
   MRSSMEFTHiggs_observables observables;

   try {
      
      observables.AMU = MRSSMEFTHiggs_a_muon::calculate_a_muon(MODEL, qedqcd);
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
