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

// File generated at Mon 27 Feb 2017 13:42:43

#include "MSSMRHN_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yv.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Yv_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = (oneOver16PiSqr*(Yv*(3*traceYuAdjYu + traceYvAdjYv - 0.6*Sqr
      (g1) - 3*Sqr(g2)) + Yv*Ye.adjoint()*Ye + 3*(Yv*Yv.adjoint()*Yv))).real();


   return beta_Yv;
}

/**
 * Calculates the two-loop beta function of Yv.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Yv_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYvYvAdjYe = TRACE_STRUCT.traceYeAdjYvYvAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = (twoLoop*(0.02*Yv*(207*Power(g1,4) + 10*Sqr(g1)*(4*
      traceYuAdjYu + 9*Sqr(g2)) + 25*(15*Power(g2,4) - 2*(3*traceYdAdjYuYuAdjYd
      + traceYeAdjYvYvAdjYe + 9*traceYuAdjYuYuAdjYu + 3*traceYvAdjYvYvAdjYv -
      16*traceYuAdjYu*Sqr(g3)))) + (-3*traceYdAdjYd - traceYeAdjYe + 1.2*Sqr(g1
      ))*(Yv*Ye.adjoint()*Ye) + (-3*(3*traceYuAdjYu + traceYvAdjYv) + 1.2*Sqr(
      g1) + 6*Sqr(g2))*(Yv*Yv.adjoint()*Yv) - 2*(Yv*Ye.adjoint()*Ye*Ye.adjoint(
      )*Ye) - 2*(Yv*Ye.adjoint()*Ye*Yv.adjoint()*Yv) - 4*(Yv*Yv.adjoint()*Yv*
      Yv.adjoint()*Yv))).real();


   return beta_Yv;
}

/**
 * Calculates the three-loop beta function of Yv.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Yv_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = ZEROMATRIX(3,3);


   return beta_Yv;
}

} // namespace flexiblesusy
