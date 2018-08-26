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

// File generated at Sun 26 Aug 2018 14:09:41

#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Ye.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Ye_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (oneOver16PiSqr*(-0.25*Ye*(9*Sqr(g1) + 9*Sqr(g2) - 2*(6*
      traceYdAdjYd + 2*traceYeAdjYe + 6*traceYuAdjYu + 3*Sqr(g2d) + 3*Sqr(g2u)
      + Sqr(gYd) + Sqr(gYu))) + 1.5*(Ye*Ye.adjoint()*Ye))).real();


   return beta_Ye;
}

/**
 * Calculates the 2-loop beta function of Ye.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Ye_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (twoLoop*(0.0025*Ye*(3006*Quad(g1) + 5*Sqr(g1)*(108*Sqr(g2) + 5*(
      10*traceYdAdjYd + 30*traceYeAdjYe + 34*traceYuAdjYu + 9*Sqr(g2d) + 9*Sqr(
      g2u) + 3*Sqr(gYd) + 3*Sqr(gYu))) - 25*(48*g2d*g2u*gYd*gYu + 108*
      traceYdAdjYdYdAdjYd - 24*traceYdAdjYuYuAdjYd + 36*traceYeAdjYeYeAdjYe +
      108*traceYuAdjYuYuAdjYu + 68*Quad(g2) + 45*Quad(g2d) + 45*Quad(g2u) + 9*
      Quad(gYd) + 9*Quad(gYu) - 320*traceYdAdjYd*Sqr(g3) - 320*traceYuAdjYu*Sqr
      (g3) + 6*Sqr(g2d)*(2*Sqr(g2u) + 3*Sqr(gYd)) + 18*Sqr(g2u)*Sqr(gYu) + 20*
      Sqr(gYd)*Sqr(gYu) - 15*Sqr(g2)*(6*traceYdAdjYd + 2*traceYeAdjYe + 6*
      traceYuAdjYu + 11*Sqr(g2d) + 11*Sqr(g2u) + Sqr(gYd) + Sqr(gYu)) - 24*Sqr(
      Lambdax))) + 0.0375*(129*Sqr(g1) + 5*(45*Sqr(g2) - 2*(18*traceYdAdjYd + 6
      *traceYeAdjYe + 18*traceYuAdjYu + 16*Lambdax + 9*Sqr(g2d) + 9*Sqr(g2u) +
      3*Sqr(gYd) + 3*Sqr(gYu))))*(Ye*Ye.adjoint()*Ye) + 1.5*(Ye*Ye.adjoint()*Ye
      *Ye.adjoint()*Ye))).real();


   return beta_Ye;
}

/**
 * Calculates the 3-loop beta function of Ye.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Ye_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

/**
 * Calculates the 4-loop beta function of Ye.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Ye_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

} // namespace flexiblesusy
