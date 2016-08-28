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

// File generated at Sun 28 Aug 2016 15:02:36

#include "HGTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Yu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(0.05*Yu*(60*traceYuAdjYu - 17*Sqr(g1) - 45*
      Sqr(g2) + 30*Sqr(g2u) + 10*Sqr(g2up) - 160*Sqr(g3)) + 0.5*(Yu*Yd.adjoint(
      )*Yd) + 1.5*(Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the two-loop beta function of Yu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Yu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(Yu*(2.305*Power(g1,4) - 3.75*Power(g2,4) - 2.8125*
      Power(g2u,4) - 0.5625*Power(g2up,4) - 94.66666666666667*Power(g3,4) +
      Lambda3*Lambda4 - 2.25*traceYdAdjYuYuAdjYd - 6.75*traceYuAdjYuYuAdjYu -
      0.45*Sqr(g1)*Sqr(g2) + 0.5625*Sqr(g1)*Sqr(g2u) - 1.125*Sqr(g1d)*Sqr(g2u)
      + 10.3125*Sqr(g2)*Sqr(g2u) + 0.1875*Sqr(g1)*Sqr(g2up) - 0.375*Sqr(g1dp)*
      Sqr(g2up) + 0.9375*Sqr(g2)*Sqr(g2up) - 1.125*Sqr(g2u)*Sqr(g2up) +
      1.2666666666666666*Sqr(g1)*Sqr(g3) + 9*Sqr(g2)*Sqr(g3) + 0.125*
      traceYuAdjYu*(17*Sqr(g1) + 45*Sqr(g2) + 160*Sqr(g3)) + 6*Sqr(Lambda2) +
      Sqr(Lambda3) + Sqr(Lambda4) + 1.5*Sqr(Lambda5) + 1.5*Sqr(Lambda6) + 4.5*
      Sqr(Lambda7)) + 0.004166666666666667*(-480*Lambda3 + 480*Lambda4 - 540*
      traceYdAdjYd - 180*traceYeAdjYe - 41*Sqr(g1) - 270*Sqr(g1d) - 90*Sqr(g1dp
      ) + 495*Sqr(g2) + 1280*Sqr(g3))*(Yu*Yd.adjoint()*Yd) + 0.0125*(-960*
      Lambda2 - 540*traceYuAdjYu + 223*Sqr(g1) + 675*Sqr(g2) - 270*Sqr(g2u) -
      90*Sqr(g2up) + 1280*Sqr(g3))*(Yu*Yu.adjoint()*Yu) - 0.25*(Yu*Yd.adjoint()
      *Yd*Yd.adjoint()*Yd) - 0.25*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) + 1.5*(
      Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the three-loop beta function of Yu.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> HGTHDMIIMSSMBC_susy_parameters::calc_beta_Yu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

} // namespace flexiblesusy
