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
 * @file cxx_qft/NUHMSSMNoFVHimalaya_vertices.cpp
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#include "NUHMSSMNoFVHimalaya_context_base.hpp"
#include "NUHMSSMNoFVHimalaya_input_parameters.hpp"
#include "NUHMSSMNoFVHimalaya_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input_parameters().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy {
namespace NUHMSSMNoFVHimalaya_cxx_diagrams {
namespace detail {

ChiralVertex unit_charge(const context_base& context)
{
   std::array<int, 0> electron_indices = {};
   std::array<int, 0> photon_indices = {};
   std::array<int, 0> indices = concatenate(photon_indices, electron_indices, electron_indices);

      const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   const std::complex<double> right = 0.7745966692414834*g1*Cos(ThetaW);

   return {left, right};
}

} // namespace detail
} // namespace NUHMSSMNoFVHimalaya_cxx_diagrams
} // namespace flexiblesusy
