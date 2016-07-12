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

// File generated at Tue 12 Jul 2016 10:31:02

#include "THDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Ye.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> THDMIIMSSMBC_susy_parameters::calc_beta_Ye_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (0.25*oneOver16PiSqr*(Ye*(12*traceYdAdjYd + 4*traceYeAdjYe -
      9*(Sqr(g1) + Sqr(g2))) + 6*(Ye*Ye.adjoint()*Ye))).real();


   return beta_Ye;
}

/**
 * Calculates the two-loop beta function of Ye.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> THDMIIMSSMBC_susy_parameters::calc_beta_Ye_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (twoLoop*(Ye*(7.245*Power(g1,4) - 5.25*Power(g2,4) + Lambda3
      *Lambda4 - 6.75*traceYdAdjYdYdAdjYd - 2.25*traceYdAdjYuYuAdjYd - 2.25*
      traceYeAdjYeYeAdjYe + 1.35*Sqr(g1)*Sqr(g2) + 1.875*traceYeAdjYe*(Sqr(g1)
      + Sqr(g2)) + 0.625*traceYdAdjYd*(Sqr(g1) + 9*Sqr(g2) + 32*Sqr(g3)) + 6*
      Sqr(Lambda1) + Sqr(Lambda3) + Sqr(Lambda4) + 1.5*Sqr(Lambda5) + 4.5*Sqr(
      Lambda6) + 1.5*Sqr(Lambda7)) + 0.0375*((-320*Lambda1 - 180*traceYdAdjYd -
      60*traceYeAdjYe + 129*Sqr(g1) + 225*Sqr(g2))*(Ye*Ye.adjoint()*Ye) + 40*(
      Ye*Ye.adjoint()*Ye*Ye.adjoint()*Ye)))).real();


   return beta_Ye;
}

/**
 * Calculates the three-loop beta function of Ye.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> THDMIIMSSMBC_susy_parameters::calc_beta_Ye_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

} // namespace flexiblesusy
