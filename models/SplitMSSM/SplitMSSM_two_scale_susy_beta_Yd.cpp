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

// File generated at Sun 10 Jan 2016 15:29:40

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

   beta_Yd = (0.25*oneOver16PiSqr*(Yd*(12*traceYdAdjYd + 4*traceYeAdjYe +
      12*traceYuAdjYu - Sqr(g1) - 9*Sqr(g2) + 6*Sqr(g2d) + 6*Sqr(g2u) - 32*Sqr
      (g3) + 2*Sqr(gYd) + 2*Sqr(gYu)) + 6*(Yd*Yd.adjoint()*Yd - Yd*Yu.adjoint()
      *Yu))).real();


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

   beta_Yd = (twoLoop*(-0.0008333333333333334*Yd*(262*Power(g1,4) + 5100*
      Power(g2,4) + 3375*Power(g2d,4) + 3375*Power(g2u,4) + 113600*Power(g3,4)
      + 675*Power(gYd,4) + 3600*g2d*g2u*gYd*gYu + 675*Power(gYu,4) + 8100*
      traceYdAdjYdYdAdjYd - 1800*traceYdAdjYuYuAdjYd + 2700*traceYeAdjYeYeAdjYe
      + 8100*traceYuAdjYuYuAdjYu - 2550*traceYuAdjYu*Sqr(g1) - 6750*
      traceYuAdjYu*Sqr(g2) + 1620*Sqr(g1)*Sqr(g2) - 2250*traceYeAdjYe*(Sqr(g1)
      + Sqr(g2)) - 675*Sqr(g1)*Sqr(g2d) - 12375*Sqr(g2)*Sqr(g2d) - 675*Sqr(g1)*
      Sqr(g2u) - 12375*Sqr(g2)*Sqr(g2u) + 900*Sqr(g2d)*Sqr(g2u) - 24000*
      traceYuAdjYu*Sqr(g3) - 2480*Sqr(g1)*Sqr(g3) - 10800*Sqr(g2)*Sqr(g3) - 750
      *traceYdAdjYd*(Sqr(g1) + 9*Sqr(g2) + 32*Sqr(g3)) - 225*Sqr(g1)*Sqr(gYd) -
      1125*Sqr(g2)*Sqr(gYd) + 1350*Sqr(g2d)*Sqr(gYd) - 225*Sqr(g1)*Sqr(gYu) -
      1125*Sqr(g2)*Sqr(gYu) + 1350*Sqr(g2u)*Sqr(gYu) + 1500*Sqr(gYd)*Sqr(gYu) -
      1800*Sqr(Lambdax)) + 0.0125*((-540*traceYdAdjYd - 180*traceYeAdjYe - 540
      *traceYuAdjYu - 480*Lambdax + 187*Sqr(g1) + 675*Sqr(g2) - 270*Sqr(g2d) -
      270*Sqr(g2u) + 1280*Sqr(g3) - 90*Sqr(gYd) - 90*Sqr(gYu))*(Yd*Yd.adjoint()
      *Yd) + (300*traceYdAdjYd + 100*traceYeAdjYe + 300*traceYuAdjYu - 79*Sqr(
      g1) + 45*Sqr(g2) + 150*Sqr(g2d) + 150*Sqr(g2u) - 1280*Sqr(g3) + 50*Sqr(
      gYd) + 50*Sqr(gYu))*(Yd*Yu.adjoint()*Yu) + 20*(6*(Yd*Yd.adjoint()*Yd*
      Yd.adjoint()*Yd) - 4*(Yd*Yd.adjoint()*Yd*Yu.adjoint()*Yu) - Yd*Yu.adjoint
      ()*Yu*Yd.adjoint()*Yd + 11*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))))).real(
      );


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
