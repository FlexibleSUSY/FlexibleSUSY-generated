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
 * @file MSSMNoFVHimalaya_edm.cpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "MSSMNoFVHimalaya_edm.hpp"
#include "cxx_qft/MSSMNoFVHimalaya_qft.hpp"
#include "MSSMNoFVHimalaya_FFV_form_factors.hpp"

namespace flexiblesusy {

using namespace MSSMNoFVHimalaya_cxx_diagrams;

namespace MSSMNoFVHimalaya_edm {

// spliting templated function into header and cpp
// https://isocpp.org/wiki/faq/templates#separate-template-fn-defn-from-decl


template <typename Lepton>
double calculate_edm(const MSSMNoFVHimalaya_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)
{
   context_base context{ model };

   using namespace MSSMNoFVHimalaya_cxx_diagrams::fields;

   const std::valarray<std::complex<double>> form_factors {0., 0., 0., 0.};

   const double val =
      // imaginary part of the axial-vector form factor
      0.5*(form_factors[3] - form_factors[2]).imag()
      // reinstate the mass that was factored out in definition of form factor
      * context.mass<Lepton>({});

   return val;
}

} // MSSMNoFVHimalaya_edm
} // namespace flexiblesusy
