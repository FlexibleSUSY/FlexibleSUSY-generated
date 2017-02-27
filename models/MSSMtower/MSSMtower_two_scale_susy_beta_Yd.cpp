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

// File generated at Mon 27 Feb 2017 13:20:15

#include "MSSMtower_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMtower_susy_parameters::calc_beta_Yd_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (oneOver16PiSqr*(Yd*(3*traceYdAdjYd + traceYeAdjYe -
      0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3)) + 3*(
      Yd*Yd.adjoint()*Yd) + Yd*Yu.adjoint()*Yu)).real();


   return beta_Yd;
}

/**
 * Calculates the two-loop beta function of Yd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMtower_susy_parameters::calc_beta_Yd_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (twoLoop*(0.011111111111111112*Yd*(287*Power(g1,4) + 2*Sqr(
      g1)*(-18*traceYdAdjYd + 54*traceYeAdjYe + 45*Sqr(g2) + 40*Sqr(g3)) + 5*(
      135*Power(g2,4) + 144*Sqr(g2)*Sqr(g3) - 2*(16*Power(g3,4) + 27*(3*
      traceYdAdjYdYdAdjYd + traceYdAdjYuYuAdjYd + traceYeAdjYeYeAdjYe) - 144*
      traceYdAdjYd*Sqr(g3)))) + (-3*(3*traceYdAdjYd + traceYeAdjYe) + 0.8*Sqr(
      g1) + 6*Sqr(g2))*(Yd*Yd.adjoint()*Yd) + (-3*traceYuAdjYu + 0.8*Sqr(g1))*(
      Yd*Yu.adjoint()*Yu) - 4*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yd*
      Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu)
      )).real();


   return beta_Yd;
}

/**
 * Calculates the three-loop beta function of Yd.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMtower_susy_parameters::calc_beta_Yd_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

} // namespace flexiblesusy
