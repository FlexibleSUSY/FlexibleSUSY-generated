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

// File generated at Thu 15 Dec 2016 12:40:44

#include "HTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> HTHDMIIMSSMBC_susy_parameters::calc_beta_Yu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(Yu*(3*traceYuAdjYu - 0.85*Sqr(g1) - 2.25*
      Sqr(g2) - 8*Sqr(g3)) + 0.5*(Yu*Yd.adjoint()*Yd) + 1.5*(Yu*Yu.adjoint()*Yu
      ))).real();


   return beta_Yu;
}

/**
 * Calculates the two-loop beta function of Yu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> HTHDMIIMSSMBC_susy_parameters::calc_beta_Yu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(Yu*(2.305*Power(g1,4) - 4.75*Power(g2,4) - 108*
      Power(g3,4) + Lambda3*Lambda4 - 2.25*traceYdAdjYuYuAdjYd - 6.75*
      traceYuAdjYuYuAdjYu + 20*traceYuAdjYu*Sqr(g3) + Sqr(g1)*(2.125*
      traceYuAdjYu - 0.45*Sqr(g2) + 1.2666666666666666*Sqr(g3)) + 1.125*Sqr(g2)
      *(5*traceYuAdjYu + 8*Sqr(g3)) + 6*Sqr(Lambda2) + Sqr(Lambda3) + Sqr(
      Lambda4) + 1.5*Sqr(Lambda5) + 1.5*Sqr(Lambda6) + 4.5*Sqr(Lambda7)) + (-2*
      Lambda3 + 2*Lambda4 - 2.25*traceYdAdjYd - 0.75*traceYeAdjYe -
      0.17083333333333334*Sqr(g1) + 2.0625*Sqr(g2) + 5.333333333333333*Sqr(g3))
      *(Yu*Yd.adjoint()*Yd) + (-12*Lambda2 - 6.75*traceYuAdjYu + 2.7875*Sqr(g1)
      + 8.4375*Sqr(g2) + 16*Sqr(g3))*(Yu*Yu.adjoint()*Yu) - 0.25*(Yu*
      Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 0.25*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*
      Yu) + 1.5*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the three-loop beta function of Yu.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> HTHDMIIMSSMBC_susy_parameters::calc_beta_Yu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

} // namespace flexiblesusy
