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


/**
 * @file MSSM_l_to_lgamma.cpp
 *
 * This file was generated with FlexibleSUSY 2.5.0 and SARAH 4.14.4 .
 */

#include <valarray>
#include <complex>

#include "MSSM_l_to_lgamma.hpp"
#include "MSSM_mass_eigenstates.hpp"

#include "cxx_qft/MSSM_qft.hpp"
#include "MSSM_FFV_form_factors.hpp"

#include "lowe.h"
#include "physical_input.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

using namespace MSSM_cxx_diagrams;
using namespace MSSM_cxx_diagrams::fields;
using namespace MSSM_FFV_form_factors;

namespace MSSM_l_to_lgamma {

template <typename FIn, typename FOut, typename T1, typename T2>
double lepton_total_decay_width(
      T1 const& indices1, T2 const& indices2, 
      const MSSM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd) {

   context_base context {model};

   const auto leptonInMass = context.mass<FIn>(indices1);
   const auto leptonOutMass = context.mass<FOut>(indices2);
   const auto r = pow(leptonOutMass/leptonInMass, 2);
   const auto GF2 = pow(qedqcd.displayFermiConstant(), 2);
   const double Alpha = Sqr(FIn::electric_charge * unit_charge(context))/(4*Pi);

   // eq. 10.6 of http://pdg.lbl.gov/2018/reviews/rpp2018-rev-standard-model.pdf
   return 
      // tree-level (with massless outgoing lepton)
      GF2*pow(leptonInMass,5)/(192.0*pow(Pi,3)) 
      // outgoing lepton mass correction
      * (1 - 8*r + 8*pow(r,3) - pow(r,4) - 12*r*r*std::log(r));
      // higher order correction
      //* (1.0 + 0.5 * Alpha * (6.25-Sqr(Pi))/Pi)
}

double calculate_Fe_to_Fe_VP (
   int generationIndex1, int generationIndex2, 
   const MSSM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Physical_input& physical_input) {

      context_base context {model};
   std::array<int, 1> indices1 = { generationIndex1 };
   std::array<int, 1> indices2 = { generationIndex2 };

   const auto form_factors = calculate_Fe_Fe_VP_form_factors (generationIndex1, generationIndex2, model, false);
   double leptonInMassOS;
   switch (generationIndex1) {
      case 0: leptonInMassOS = qedqcd.displayMass(softsusy::mElectron); break;
      case 1: leptonInMassOS = qedqcd.displayMass(softsusy::mMuon);     break;
      case 2: leptonInMassOS = qedqcd.displayMass(softsusy::mTau);      break;
      default: throw std::invalid_argument("Unrecognized lepton");
   }

   // eq. 51 of arXiv:hep-ph/9510309 (note that we include 'e' in the definition of form_factor)
   const double partial_width = pow(leptonInMassOS,5)/(16.0*Pi) * (std::norm(form_factors[2]) + std::norm(form_factors[3]));
   const double total_width = lepton_total_decay_width<Fe, Fe>(indices1, indices2, model, qedqcd);

   return partial_width/total_width;
}

}
} // namespace flexiblesusy
