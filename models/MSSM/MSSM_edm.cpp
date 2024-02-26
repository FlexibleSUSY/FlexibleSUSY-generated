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
 * @file MSSM_edm.cpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "MSSM_edm.hpp"
#include "cxx_qft/MSSM_qft.hpp"
#include "MSSM_FFV_form_factors.hpp"

namespace flexiblesusy {

using namespace MSSM_cxx_diagrams;

namespace MSSM_edm {

// spliting templated function into header and cpp
// https://isocpp.org/wiki/faq/templates#separate-template-fn-defn-from-decl
template double MSSM_edm::calculate_edm<MSSM_cxx_diagrams::fields::Fe>(const MSSM_mass_eigenstates&, const softsusy::QedQcd&, int);

template <typename Lepton>
double calculate_edm(const MSSM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, int idx)
{
   context_base context{ model };

   using namespace MSSM_cxx_diagrams::fields;

   const auto form_factors = MSSM_FFV_form_factors::calculate_form_factors<Lepton,Lepton,VP>(idx, idx, model, false);

   const double val =
      // imaginary part of the axial-vector form factor
      0.5*(form_factors[3] - form_factors[2]).imag()
      // reinstate the mass that was factored out in definition of form factor
      * context.mass<Lepton>({idx});

   return val;
}

} // MSSM_edm
} // namespace flexiblesusy
