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

// File generated at Wed 12 Apr 2017 11:14:17

#include "SplitMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yd_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (oneOver16PiSqr*(-0.25*Yd*(Sqr(g1) + 9*Sqr(g2) - 2*(6*
      traceYdAdjYd + 2*traceYeAdjYe + 6*traceYuAdjYu + 3*Sqr(g2d) + 3*Sqr(g2u)
      - 16*Sqr(g3) + Sqr(gYd) + Sqr(gYu))) + 1.5*(Yd*Yd.adjoint()*Yd) - 1.5*(Yd
      *Yu.adjoint()*Yu))).real();


   return beta_Yd;
}

/**
 * Calculates the two-loop beta function of Yd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yd_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (twoLoop*(-0.0008333333333333334*Yd*(262*Power(g1,4) + 5*Sqr
      (g1)*(-150*traceYdAdjYd - 450*traceYeAdjYe - 510*traceYuAdjYu + 324*Sqr(
      g2) - 135*Sqr(g2d) - 135*Sqr(g2u) - 496*Sqr(g3) - 45*Sqr(gYd) - 45*Sqr(
      gYu)) + 25*(204*Power(g2,4) + 135*Power(g2d,4) + 135*Power(g2u,4) + 4544*
      Power(g3,4) + 27*Power(gYd,4) + 144*g2d*g2u*gYd*gYu + 27*Power(gYu,4) +
      324*traceYdAdjYdYdAdjYd - 72*traceYdAdjYuYuAdjYd + 108*
      traceYeAdjYeYeAdjYe + 324*traceYuAdjYuYuAdjYu - 960*traceYdAdjYd*Sqr(g3)
      - 960*traceYuAdjYu*Sqr(g3) + 18*Sqr(g2d)*(2*Sqr(g2u) + 3*Sqr(gYd)) + 54*
      Sqr(g2u)*Sqr(gYu) + 60*Sqr(gYd)*Sqr(gYu) - 9*Sqr(g2)*(30*traceYdAdjYd +
      10*traceYeAdjYe + 30*traceYuAdjYu + 55*Sqr(g2d) + 55*Sqr(g2u) + 48*Sqr(g3
      ) + 5*Sqr(gYd) + 5*Sqr(gYu)) - 72*Sqr(Lambdax))) + 0.0125*(187*Sqr(g1) +
      5*(135*Sqr(g2) - 2*(54*traceYdAdjYd + 18*traceYeAdjYe + 54*traceYuAdjYu +
      48*Lambdax + 27*Sqr(g2d) + 27*Sqr(g2u) - 128*Sqr(g3) + 9*Sqr(gYd) + 9*
      Sqr(gYu))))*(Yd*Yd.adjoint()*Yd) + 0.0125*(-79*Sqr(g1) + 5*(60*
      traceYdAdjYd + 20*traceYeAdjYe + 60*traceYuAdjYu + 9*Sqr(g2) + 30*Sqr(g2d
      ) + 30*Sqr(g2u) - 256*Sqr(g3) + 10*Sqr(gYd) + 10*Sqr(gYu)))*(Yd*
      Yu.adjoint()*Yu) + 1.5*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - Yd*
      Yd.adjoint()*Yd*Yu.adjoint()*Yu - 0.25*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*
      Yd) + 2.75*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yd;
}

/**
 * Calculates the three-loop beta function of Yd.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yd_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

} // namespace flexiblesusy
