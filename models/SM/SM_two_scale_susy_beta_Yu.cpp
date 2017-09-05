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

// File generated at Tue 5 Sep 2017 10:40:00

#include "SM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> SM_susy_parameters::calc_beta_Yu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(Yu*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu - 0.85*Sqr(g1) - 2.25*Sqr(g2) - 8*Sqr(g3)) - 1.5*(Yu*
      Yd.adjoint()*Yd) + 1.5*(Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the two-loop beta function of Yu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> SM_susy_parameters::calc_beta_Yu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(0.0016666666666666668*Yu*(1187*Power(g1,4) + 5*Sqr
      (g1)*(75*traceYdAdjYd + 225*traceYeAdjYe + 255*traceYuAdjYu - 54*Sqr(g2)
      + 152*Sqr(g3)) - 75*(46*Power(g2,4) - 3*Sqr(g2)*(5*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu) + 24*Sqr(g3)) + 2*(432*Power(g3,4) - 80*(
      traceYdAdjYd + traceYuAdjYu)*Sqr(g3) + 3*(9*traceYdAdjYdYdAdjYd - 2*
      traceYdAdjYuYuAdjYd + 3*traceYeAdjYeYeAdjYe + 9*traceYuAdjYuYuAdjYu - 2*
      Sqr(Lambdax))))) + 0.0125*(-43*Sqr(g1) + 5*(20*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu) + 9*Sqr(g2) - 256*Sqr(g3)))*(Yu*Yd.adjoint
      ()*Yd) + 0.0125*(223*Sqr(g1) + 675*Sqr(g2) + 20*(-3*(9*traceYdAdjYd + 3*
      traceYeAdjYe + 9*traceYuAdjYu + 8*Lambdax) + 64*Sqr(g3)))*(Yu*Yu.adjoint(
      )*Yu) + 2.75*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 0.25*(Yu*Yd.adjoint()
      *Yd*Yu.adjoint()*Yu) - Yu*Yu.adjoint()*Yu*Yd.adjoint()*Yd + 1.5*(Yu*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the three-loop beta function of Yu.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> SM_susy_parameters::calc_beta_Yu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (0.00005*PROJECTOR*threeLoop*Yu(2,2)*(321980*Power(g1,6) +
      3396580*Power(g2,6) - 5*Power(g2,4)*(21375*Lambdax + 84288*Sqr(g3) -
      67960*Sqr(Yu(2,2))) - 5*Power(g1,4)*(5445*Lambdax + 17768*Sqr(g2) + 89276
      *Sqr(g3) + 97688*Sqr(Yu(2,2))) - 10*Sqr(g1)*(9486*Power(g2,4) + 30192*
      Power(g3,4) - 4500*Sqr(Lambdax) + Sqr(g2)*(-2925*Lambdax + 32100*Sqr(g3)
      - 69658*Sqr(Yu(2,2))) + 12700*Lambdax*Sqr(Yu(2,2)) - 36148*Sqr(g3)*Sqr(Yu
      (2,2)) + 60925*Power(Yu(2,2),4)) + 10*Sqr(g2)*(147308*Power(g3,4) + 96740
      *Sqr(g3)*Sqr(Yu(2,2)) + 1125*(20*Sqr(Lambdax) - 60*Lambdax*Sqr(Yu(2,2)) -
      177*Power(Yu(2,2),4))) + 2*(-6193500*Power(g3,6) - 45000*Power(Lambdax,3
      ) + 3637640*Power(g3,4)*Sqr(Yu(2,2)) + 9375*Sqr(Lambdax)*Sqr(Yu(2,2)) +
      990000*Lambdax*Power(Yu(2,2),4) + 586028*Power(Yu(2,2),6) + 10000*Sqr(g3)
      *(8*Lambdax*Sqr(Yu(2,2)) - 157*Power(Yu(2,2),4))))).real();


   return beta_Yu;
}

} // namespace flexiblesusy
