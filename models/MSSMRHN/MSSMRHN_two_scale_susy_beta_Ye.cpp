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

// File generated at Mon 9 May 2016 13:20:47

#include "MSSMRHN_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Ye.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Ye_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (oneOver16PiSqr*(Ye*(3*traceYdAdjYd + traceYeAdjYe - 1.8*Sqr
      (g1) - 3*Sqr(g2)) + 3*(Ye*Ye.adjoint()*Ye) + Ye*Yv.adjoint()*Yv)).real();


   return beta_Ye;
}

/**
 * Calculates the two-loop beta function of Ye.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Ye_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYvYvAdjYe = TRACE_STRUCT.traceYeAdjYvYvAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (twoLoop*(0.1*Ye*(135*Power(g1,4) + 75*Power(g2,4) - 90*
      traceYdAdjYdYdAdjYd - 30*traceYdAdjYuYuAdjYd - 30*traceYeAdjYeYeAdjYe -
      10*traceYeAdjYvYvAdjYe + 12*traceYeAdjYe*Sqr(g1) + 18*Sqr(g1)*Sqr(g2) - 4
      *traceYdAdjYd*(Sqr(g1) - 40*Sqr(g3))) + (-9*traceYdAdjYd - 3*traceYeAdjYe
      + 6*Sqr(g2))*(Ye*Ye.adjoint()*Ye) - 3*traceYuAdjYu*(Ye*Yv.adjoint()*Yv)
      - traceYvAdjYv*(Ye*Yv.adjoint()*Yv) - 4*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*
      Ye) - 2*(Ye*Yv.adjoint()*Yv*Ye.adjoint()*Ye) - 2*(Ye*Yv.adjoint()*Yv*
      Yv.adjoint()*Yv))).real();


   return beta_Ye;
}

/**
 * Calculates the three-loop beta function of Ye.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Ye_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

} // namespace flexiblesusy
