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

// File generated at Mon 9 May 2016 12:04:17

#include "HGTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Ye.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Ye_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (0.25*oneOver16PiSqr*(Ye*(12*traceYdAdjYd + 4*traceYeAdjYe -
      9*Sqr(g1) + 6*Sqr(g1d) + 2*Sqr(g1dp) - 9*Sqr(g2)) + 6*(Ye*Ye.adjoint()*
      Ye))).real();


   return beta_Ye;
}

/**
 * Calculates the two-loop beta function of Ye.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Ye_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (0.0025*twoLoop*(Ye*(3162*Power(g1,4) - 1125*Power(g1d,4) -
      225*Power(g1dp,4) - 1500*Power(g2,4) + 400*Lambda3*Lambda4 - 2700*
      traceYdAdjYdYdAdjYd - 900*traceYdAdjYuYuAdjYd - 900*traceYeAdjYeYeAdjYe +
      225*Sqr(g1)*Sqr(g1d) + 75*Sqr(g1)*Sqr(g1dp) - 450*Sqr(g1d)*Sqr(g1dp) +
      540*Sqr(g1)*Sqr(g2) + 4125*Sqr(g1d)*Sqr(g2) + 375*Sqr(g1dp)*Sqr(g2) + 750
      *traceYeAdjYe*(Sqr(g1) + Sqr(g2)) - 450*Sqr(g1d)*Sqr(g2u) - 150*Sqr(g1dp)
      *Sqr(g2up) + 250*traceYdAdjYd*(Sqr(g1) + 9*Sqr(g2) + 32*Sqr(g3)) + 2400*
      Sqr(Lambda1) + 400*Sqr(Lambda3) + 400*Sqr(Lambda4) + 600*Sqr(Lambda5) +
      1800*Sqr(Lambda6) + 600*Sqr(Lambda7)) + 15*((-320*Lambda1 - 180*
      traceYdAdjYd - 60*traceYeAdjYe + 129*Sqr(g1) - 90*Sqr(g1d) - 30*Sqr(g1dp)
      + 225*Sqr(g2))*(Ye*Ye.adjoint()*Ye) + 40*(Ye*Ye.adjoint()*Ye*Ye.adjoint(
      )*Ye)))).real();


   return beta_Ye;
}

/**
 * Calculates the three-loop beta function of Ye.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Ye_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

} // namespace flexiblesusy
