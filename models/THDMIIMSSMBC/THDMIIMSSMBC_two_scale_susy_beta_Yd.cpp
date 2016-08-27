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

// File generated at Sat 27 Aug 2016 11:41:41

#include "THDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> THDMIIMSSMBC_susy_parameters::calc_beta_Yd_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (oneOver16PiSqr*(-0.25*Yd*(-12*traceYdAdjYd - 4*traceYeAdjYe
      + Sqr(g1) + 9*Sqr(g2) + 32*Sqr(g3)) + 1.5*(Yd*Yd.adjoint()*Yd) + 0.5*(Yd
      *Yu.adjoint()*Yu))).real();


   return beta_Yd;
}

/**
 * Calculates the two-loop beta function of Yd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> THDMIIMSSMBC_susy_parameters::calc_beta_Yd_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (twoLoop*(Yd*(-0.18833333333333332*Power(g1,4) - 5.25*Power(
      g2,4) - 108*Power(g3,4) + Lambda3*Lambda4 - 6.75*traceYdAdjYdYdAdjYd -
      2.25*traceYdAdjYuYuAdjYd - 2.25*traceYeAdjYeYeAdjYe - 1.35*Sqr(g1)*Sqr(g2
      ) + 1.875*traceYeAdjYe*(Sqr(g1) + Sqr(g2)) + 2.066666666666667*Sqr(g1)*
      Sqr(g3) + 9*Sqr(g2)*Sqr(g3) + 0.625*traceYdAdjYd*(Sqr(g1) + 9*Sqr(g2) +
      32*Sqr(g3)) + 6*Sqr(Lambda1) + Sqr(Lambda3) + Sqr(Lambda4) + 1.5*Sqr(
      Lambda5) + 4.5*Sqr(Lambda6) + 1.5*Sqr(Lambda7)) + 0.0125*(-960*Lambda1 -
      540*traceYdAdjYd - 180*traceYeAdjYe + 187*Sqr(g1) + 675*Sqr(g2) + 1280*
      Sqr(g3))*(Yd*Yd.adjoint()*Yd) + 0.004166666666666667*(-540*traceYuAdjYu -
      53*Sqr(g1) + 5*(-96*Lambda3 + 96*Lambda4 + 99*Sqr(g2) + 256*Sqr(g3)))*(
      Yd*Yu.adjoint()*Yu) + 1.5*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 0.25*(Yd
      *Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 0.25*(Yd*Yu.adjoint()*Yu*Yu.adjoint()
      *Yu))).real();


   return beta_Yd;
}

/**
 * Calculates the three-loop beta function of Yd.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> THDMIIMSSMBC_susy_parameters::calc_beta_Yd_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

} // namespace flexiblesusy
