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
 * @file MRSSM2_edm.cpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "MRSSM2_edm.hpp"
#include "cxx_qft/MRSSM2_qft.hpp"
#include "MRSSM2_FFV_form_factors.hpp"

namespace flexiblesusy {

using namespace MRSSM2_cxx_diagrams;

namespace MRSSM2_edm {

// spliting templated function into header and cpp
// https://isocpp.org/wiki/faq/templates#separate-template-fn-defn-from-decl
template double MRSSM2_edm::calculate_edm<MRSSM2_cxx_diagrams::fields::Fe>(const MRSSM2_mass_eigenstates&, const softsusy::QedQcd&, int);

template <typename Lepton>
double calculate_edm(const MRSSM2_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, int idx)
{
   context_base context{ model };

   using namespace MRSSM2_cxx_diagrams::fields;

   const auto form_factors = MRSSM2_FFV_form_factors::calculate_form_factors<Lepton,Lepton,VP>(idx, idx, model, false);

   const double val =
      // imaginary part of the axial-vector form factor
      0.5*(form_factors[3] - form_factors[2]).imag()
      // reinstate the mass that was factored out in definition of form factor
      * context.mass<Lepton>({idx});

   return val;
}

} // MRSSM2_edm
} // namespace flexiblesusy
