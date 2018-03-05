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

// File generated at Mon 5 Mar 2018 18:45:56

#include "MSSMRHN_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Mv.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Mv_1_loop(const Susy_traces& susy_traces) const
{


   Eigen::Matrix<double,3,3> beta_Mv;

   beta_Mv = (oneOver16PiSqr*(2*(Mv*Yv.conjugate()*Yv.transpose()) + 2*(
      Yv*Yv.adjoint()*Mv))).real();


   return beta_Mv;
}

/**
 * Calculates the 2-loop beta function of Mv.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Mv_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_Mv;

   beta_Mv = (twoLoop*((-2*(3*traceYuAdjYu + traceYvAdjYv) + 1.2*Sqr(g1)
      + 6*Sqr(g2))*(Mv*Yv.conjugate()*Yv.transpose()) + (-2*(3*traceYuAdjYu +
      traceYvAdjYv) + 1.2*Sqr(g1) + 6*Sqr(g2))*(Yv*Yv.adjoint()*Mv) - 2*(Mv*
      Yv.conjugate()*Ye.transpose()*Ye.conjugate()*Yv.transpose()) - 2*(Mv*
      Yv.conjugate()*Yv.transpose()*Yv.conjugate()*Yv.transpose()) - 2*(Yv*
      Ye.adjoint()*Ye*Yv.adjoint()*Mv) - 2*(Yv*Yv.adjoint()*Yv*Yv.adjoint()*Mv)
      )).real();


   return beta_Mv;
}

/**
 * Calculates the 3-loop beta function of Mv.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Mv_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Mv;

   beta_Mv = ZEROMATRIX(3,3);


   return beta_Mv;
}

/**
 * Calculates the 4-loop beta function of Mv.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Mv_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Mv;

   beta_Mv = ZEROMATRIX(3,3);


   return beta_Mv;
}

} // namespace flexiblesusy
