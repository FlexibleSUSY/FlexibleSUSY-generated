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

// File generated at Wed 12 Apr 2017 11:25:50

#include "E6SSMtower_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of MuPr.
 *
 * @return one-loop beta function
 */
double E6SSMtower_susy_parameters::calc_beta_MuPr_one_loop(const Susy_traces& susy_traces) const
{


   double beta_MuPr;

   beta_MuPr = Re(-0.2*oneOver16PiSqr*MuPr*(3*Sqr(g1) + 15*Sqr(g2) + 2*
      Sqr(gN)));


   return beta_MuPr;
}

/**
 * Calculates the two-loop beta function of MuPr.
 *
 * @return two-loop beta function
 */
double E6SSMtower_susy_parameters::calc_beta_MuPr_two_loop(const Susy_traces& susy_traces) const
{


   double beta_MuPr;

   beta_MuPr = Re(0.06*twoLoop*MuPr*(99*Power(g1,4) + 275*Power(g2,4) +
      64*Power(gN,4) + 20*Sqr(g2)*Sqr(gN) + 6*Sqr(g1)*(5*Sqr(g2) + 2*Sqr(gN))))
      ;


   return beta_MuPr;
}

/**
 * Calculates the three-loop beta function of MuPr.
 *
 * @return three-loop beta function
 */
double E6SSMtower_susy_parameters::calc_beta_MuPr_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MuPr;

   beta_MuPr = 0;


   return beta_MuPr;
}

} // namespace flexiblesusy
