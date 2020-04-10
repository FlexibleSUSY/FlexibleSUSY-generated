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

// File generated at Fri 10 Apr 2020 20:34:55

#include "MSSMRHN_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yv.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Yv_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = (oneOver16PiSqr*(-0.2*Yv*(-15*traceYuAdjYu - 5*traceYvAdjYv + 3*
      Sqr(g1) + 15*Sqr(g2)) + Yv*Ye.adjoint()*Ye + 3*(Yv*Yv.adjoint()*Yv))).
      real();


   return beta_Yv;
}

/**
 * Calculates the 2-loop beta function of Yv.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Yv_2_loop(const Susy_traces& susy_traces) const
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

   beta_Yv = (twoLoop*(0.02*Yv*(-150*traceYdAdjYuYuAdjYd - 50*
      traceYeAdjYvYvAdjYe - 450*traceYuAdjYuYuAdjYu - 150*traceYvAdjYvYvAdjYv +
      207*Quad(g1) + 375*Quad(g2) + 40*traceYuAdjYu*Sqr(g1) + 90*Sqr(g1)*Sqr(g2
      ) + 800*traceYuAdjYu*Sqr(g3)) + 0.2*(-15*traceYdAdjYd - 5*traceYeAdjYe +
      6*Sqr(g1))*(Yv*Ye.adjoint()*Ye) + 0.6*(-15*traceYuAdjYu - 5*traceYvAdjYv
      + 2*Sqr(g1) + 10*Sqr(g2))*(Yv*Yv.adjoint()*Yv) - 2*(Yv*Ye.adjoint()*Ye*Ye
      .adjoint()*Ye) - 2*(Yv*Ye.adjoint()*Ye*Yv.adjoint()*Yv) - 4*(Yv*Yv.
      adjoint()*Yv*Yv.adjoint()*Yv))).real();


   return beta_Yv;
}

/**
 * Calculates the 3-loop beta function of Yv.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Yv_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = ZEROMATRIX(3,3);


   return beta_Yv;
}

/**
 * Calculates the 4-loop beta function of Yv.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Yv_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = ZEROMATRIX(3,3);


   return beta_Yv;
}

/**
 * Calculates the 5-loop beta function of Yv.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMRHN_susy_parameters::calc_beta_Yv_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = ZEROMATRIX(3,3);


   return beta_Yv;
}

} // namespace flexiblesusy
