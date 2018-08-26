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

// File generated at Sun 26 Aug 2018 14:09:40

#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yu.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yu_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(-0.05*Yu*(17*Sqr(g1) + 5*(9*Sqr(g2) - 2*(6*
      traceYdAdjYd + 2*traceYeAdjYe + 6*traceYuAdjYu + 3*Sqr(g2d) + 3*Sqr(g2u)
      - 16*Sqr(g3) + Sqr(gYd) + Sqr(gYu)))) - 1.5*(Yu*Yd.adjoint()*Yd) + 1.5*(
      Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the 2-loop beta function of Yu.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yu_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(0.0008333333333333334*Yu*(2606*Quad(g1) + 5*Sqr(g1)*(150
      *traceYdAdjYd + 450*traceYeAdjYe + 510*traceYuAdjYu - 108*Sqr(g2) + 135*
      Sqr(g2d) + 135*Sqr(g2u) + 304*Sqr(g3) + 45*Sqr(gYd) + 45*Sqr(gYu)) - 25*(
      144*g2d*g2u*gYd*gYu + 324*traceYdAdjYdYdAdjYd - 72*traceYdAdjYuYuAdjYd +
      108*traceYeAdjYeYeAdjYe + 324*traceYuAdjYuYuAdjYu + 204*Quad(g2) + 135*
      Quad(g2d) + 135*Quad(g2u) + 4544*Quad(g3) + 27*Quad(gYd) + 27*Quad(gYu) -
      960*traceYdAdjYd*Sqr(g3) - 960*traceYuAdjYu*Sqr(g3) + 18*Sqr(g2d)*(2*Sqr(
      g2u) + 3*Sqr(gYd)) + 54*Sqr(g2u)*Sqr(gYu) + 60*Sqr(gYd)*Sqr(gYu) - 9*Sqr(
      g2)*(30*traceYdAdjYd + 10*traceYeAdjYe + 30*traceYuAdjYu + 55*Sqr(g2d) +
      55*Sqr(g2u) + 48*Sqr(g3) + 5*Sqr(gYd) + 5*Sqr(gYu)) - 72*Sqr(Lambdax))) +
      0.0125*(-43*Sqr(g1) + 5*(60*traceYdAdjYd + 20*traceYeAdjYe + 60*
      traceYuAdjYu + 9*Sqr(g2) + 30*Sqr(g2d) + 30*Sqr(g2u) - 256*Sqr(g3) + 10*
      Sqr(gYd) + 10*Sqr(gYu)))*(Yu*Yd.adjoint()*Yd) + 0.0125*(223*Sqr(g1) + 5*(
      135*Sqr(g2) - 2*(54*traceYdAdjYd + 18*traceYeAdjYe + 54*traceYuAdjYu + 48
      *Lambdax + 27*Sqr(g2d) + 27*Sqr(g2u) - 128*Sqr(g3) + 9*Sqr(gYd) + 9*Sqr(
      gYu))))*(Yu*Yu.adjoint()*Yu) + 2.75*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd)
      - 0.25*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) - Yu*Yu.adjoint()*Yu*Yd.
      adjoint()*Yd + 1.5*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the 3-loop beta function of Yu.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

/**
 * Calculates the 4-loop beta function of Yu.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

} // namespace flexiblesusy
