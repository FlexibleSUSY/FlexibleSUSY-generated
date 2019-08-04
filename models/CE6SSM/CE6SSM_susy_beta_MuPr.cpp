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

// File generated at Sun 4 Aug 2019 17:10:43

#include "CE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of MuPr.
 *
 * @return 1-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_MuPr_1_loop(const Susy_traces& susy_traces) const
{


   double beta_MuPr;

   beta_MuPr = Re(-0.2*oneOver16PiSqr*MuPr*(3*Sqr(g1) + 15*Sqr(g2) + 2*Sqr(gN))
      );


   return beta_MuPr;
}

/**
 * Calculates the 2-loop beta function of MuPr.
 *
 * @return 2-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_MuPr_2_loop(const Susy_traces& susy_traces) const
{


   double beta_MuPr;

   beta_MuPr = Re(0.06*twoLoop*MuPr*(99*Quad(g1) + 275*Quad(g2) + 64*Quad(gN) +
      30*Sqr(g1)*Sqr(g2) + 12*Sqr(g1)*Sqr(gN) + 20*Sqr(g2)*Sqr(gN)));


   return beta_MuPr;
}

/**
 * Calculates the 3-loop beta function of MuPr.
 *
 * @return 3-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_MuPr_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MuPr;

   beta_MuPr = 0;


   return beta_MuPr;
}

/**
 * Calculates the 4-loop beta function of MuPr.
 *
 * @return 4-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_MuPr_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MuPr;

   beta_MuPr = 0;


   return beta_MuPr;
}

/**
 * Calculates the 5-loop beta function of MuPr.
 *
 * @return 5-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_MuPr_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MuPr;

   beta_MuPr = 0;


   return beta_MuPr;
}

} // namespace flexiblesusy
