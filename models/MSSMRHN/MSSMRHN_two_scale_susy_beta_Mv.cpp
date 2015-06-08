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

// File generated at Mon 8 Jun 2015 18:00:20

#include "MSSMRHN_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Mv.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Mv_one_loop(const Susy_traces& susy_traces) const
{


   Eigen::Matrix<double,3,3> beta_Mv;

   beta_Mv = (2*oneOver16PiSqr*(Mv*Yv.conjugate()*Yv.transpose() + Yv*
      Yv.adjoint()*Mv)).real();


   return beta_Mv;
}

/**
 * Calculates the two-loop beta function of Mv.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Mv_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_Mv;

   beta_Mv = (twoLoop*((-6*traceYuAdjYu - 2*traceYvAdjYv + 1.2*Sqr(g1) +
      6*Sqr(g2))*(Mv*Yv.conjugate()*Yv.transpose()) + (-6*traceYuAdjYu - 2*
      traceYvAdjYv + 1.2*Sqr(g1) + 6*Sqr(g2))*(Yv*Yv.adjoint()*Mv) - 2*(Mv*
      Yv.conjugate()*Ye.transpose()*Ye.conjugate()*Yv.transpose() + Mv*
      Yv.conjugate()*Yv.transpose()*Yv.conjugate()*Yv.transpose() + Yv*
      Ye.adjoint()*Ye*Yv.adjoint()*Mv + Yv*Yv.adjoint()*Yv*Yv.adjoint()*Mv)))
      .real();


   return beta_Mv;
}

} // namespace flexiblesusy
