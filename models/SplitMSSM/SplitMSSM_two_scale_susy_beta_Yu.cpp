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

// File generated at Fri 8 Jan 2016 11:56:19

#include "SplitMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(Yu*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu - 0.85*Sqr(g1) - 2.25*Sqr(g2) + 1.5*Sqr(g2d) + 1.5*Sqr(g2u)
      - 8*Sqr(g3) + 0.5*Sqr(gYd) + 0.5*Sqr(gYu)) - 1.5*(Yu*Yd.adjoint()*Yd - Yu
      *Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the two-loop beta function of Yu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (0.0008333333333333334*twoLoop*(Yu*(2606*Power(g1,4) - 5100*
      Power(g2,4) - 3375*Power(g2d,4) - 3375*Power(g2u,4) - 113600*Power(g3,4)
      - 675*Power(gYd,4) - 3600*g2d*g2u*gYd*gYu - 675*Power(gYu,4) - 8100*
      traceYdAdjYdYdAdjYd + 1800*traceYdAdjYuYuAdjYd - 2700*traceYeAdjYeYeAdjYe
      - 8100*traceYuAdjYuYuAdjYu + 2550*traceYuAdjYu*Sqr(g1) + 6750*
      traceYuAdjYu*Sqr(g2) - 540*Sqr(g1)*Sqr(g2) + 2250*traceYeAdjYe*(Sqr(g1) +
      Sqr(g2)) + 675*Sqr(g1)*Sqr(g2d) + 12375*Sqr(g2)*Sqr(g2d) + 675*Sqr(g1)*
      Sqr(g2u) + 12375*Sqr(g2)*Sqr(g2u) - 900*Sqr(g2d)*Sqr(g2u) + 24000*
      traceYuAdjYu*Sqr(g3) + 1520*Sqr(g1)*Sqr(g3) + 10800*Sqr(g2)*Sqr(g3) + 750
      *traceYdAdjYd*(Sqr(g1) + 9*Sqr(g2) + 32*Sqr(g3)) + 225*Sqr(g1)*Sqr(gYd) +
      1125*Sqr(g2)*Sqr(gYd) - 1350*Sqr(g2d)*Sqr(gYd) + 225*Sqr(g1)*Sqr(gYu) +
      1125*Sqr(g2)*Sqr(gYu) - 1350*Sqr(g2u)*Sqr(gYu) - 1500*Sqr(gYd)*Sqr(gYu) +
      1800*Sqr(Lambdax)) - 15*((-300*traceYdAdjYd - 100*traceYeAdjYe - 300*
      traceYuAdjYu + 43*Sqr(g1) - 45*Sqr(g2) - 150*Sqr(g2d) - 150*Sqr(g2u) +
      1280*Sqr(g3) - 50*Sqr(gYd) - 50*Sqr(gYu))*(Yu*Yd.adjoint()*Yd) + (540*
      traceYdAdjYd + 180*traceYeAdjYe + 540*traceYuAdjYu + 480*Lambdax - 223*
      Sqr(g1) - 675*Sqr(g2) + 270*Sqr(g2d) + 270*Sqr(g2u) - 1280*Sqr(g3) + 90*
      Sqr(gYd) + 90*Sqr(gYu))*(Yu*Yu.adjoint()*Yu) - 20*(11*(Yu*Yd.adjoint()*Yd
      *Yd.adjoint()*Yd) - Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu - 4*(Yu*Yu.adjoint
      ()*Yu*Yd.adjoint()*Yd) + 6*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu))))).real(
      );


   return beta_Yu;
}

/**
 * Calculates the three-loop beta function of Yu.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> SplitMSSM_susy_parameters::calc_beta_Yu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

} // namespace flexiblesusy
